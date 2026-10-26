#!/usr/bin/env bash
# Dry run of tools/kestrel/nightly_gpu.sh with stub git, module, sbatch,
# srun, cmake and ctest.
#
#   nightly_gpu_dryrun.sh SRC WORK
#
# run: refused outside an allocation; inside, a passing run configures one
#   rank per GPU, runs `ctest -L smoke` with a JUnit file and writes "pass";
#   a smoke run of gust without the GPU flags stops the job before the tests
#   ("fail: arena flags"); failing tests give "fail: tests" and exit 1; old
#   run logs are pruned to GUST_NIGHTLY_KEEP.
# submit: needs GUST_ACCOUNT; the job script requests the GPUs and calls run.
# A bad GUST_NIGHTLY_GPUS is refused. A failed fetch writes "fail: fetch"; a
# second run while one holds the work directory's lock is refused.
set -euo pipefail
src=$1; work=$2
script="$src/tools/kestrel/nightly_gpu.sh"

rm -rf "$work"
mkdir -p "$work/bin"
bin="$work/bin"
for stub in module srun; do
    printf '#!/bin/sh\necho "%s stub: $*"\n' "$stub" > "$bin/$stub"
done
cat > "$bin/sbatch" <<'EOF'
#!/bin/sh
echo "sbatch stub: $*"
cp "$1" "$STUB_LOG/submitted_job.sh"
EOF
# git clone makes a checkout whose build scripts record their calls.
cat > "$bin/git" <<'EOF'
#!/bin/bash
echo "git $*" >> "$STUB_LOG/git.log"
case " $* " in *" ${STUB_GIT_FAIL:-none} "*) exit 128 ;; esac
if [ "$1" = clone ]; then
    mkdir -p "$3/.git" "$3/tools/kestrel"
    printf '#!/bin/sh\necho "build_kestrel $* jobs=$GUST_BUILD_JOBS" >> "$STUB_LOG/build.log"\nmkdir -p "$(dirname "$0")/../../build-kestrel-gpu"\n' \
        > "$3/tools/kestrel/build_kestrel.sh"
    chmod +x "$3/tools/kestrel/build_kestrel.sh"
    echo 'echo "kestrel_env $*" >> "$STUB_LOG/build.log"' > "$3/tools/kestrel/kestrel_env.sh"
    exit 0
fi
[ "$1" = -C ] && shift 2
case "$1" in
    rev-parse) if [ "$2" = HEAD ]; then echo 0123456789abcdef0123456789abcdef01234567; else echo 012345678; fi ;;
esac
EOF
cat > "$bin/cmake" <<'EOF'
#!/bin/sh
echo "cmake $*" >> "$STUB_LOG/cmake.log"
EOF
cat > "$bin/ctest" <<'EOF'
#!/bin/sh
case "$*" in
    *--show-only=json-v1*) cat "$STUB_TESTS_JSON"; exit 0 ;;
esac
echo "ctest $*" >> "$STUB_LOG/ctest.log"
exit "${STUB_CTEST_STATUS:-0}"
EOF
chmod +x "$bin"/*

flags='"amrex.the_arena_init_size=0", "amrex.init_snan=1", "amrex.fpe_trap_invalid=1"'
cat > "$work/tests_good.json" <<EOF
{"tests": [
 {"name": "smoke_run", "command": ["srun", "-n", "2", "--gpus-per-task=1", "/b/gust", "inputs_smoke", $flags]},
 {"name": "rt_coriolis_repeat_plt_matches", "command": ["srun", "-n", "1", "/b/fcompare", "a", "b"]}]}
EOF
cat > "$work/tests_bad.json" <<EOF
{"tests": [
 {"name": "smoke_run", "command": ["srun", "-n", "2", "/b/gust", "inputs_smoke", $flags]},
 {"name": "rt_no_flags", "command": ["srun", "-n", "2", "/b/gust", "inputs_x", "amrex.init_snan=1"]}]}
EOF

# On Kestrel module (and ml) are exported shell functions, which a script
# finds before the stub on PATH, and every bash script defines them again
# from BASH_ENV (Lmod's init): a real `module reset` would unload the loaded
# modules of the caller (ctest on a GPU node).
unset -f module ml 2>/dev/null || true
unset BASH_ENV
# The job's own inputs from the caller's environment would change what it does.
unset GUST_BUILD_JOBS GUST_NIGHTLY_REPO GUST_NIGHTLY_REF GUST_NIGHTLY_GPUS GUST_NIGHTLY_KEEP GUST_NIGHTLY_MAIL GUST_ACCOUNT
export PATH="$bin:$PATH" STUB_LOG="$work/stublog" GUST_NIGHTLY_WORK="$work/nightly" SLURM_CPUS_ON_NODE=32
mkdir -p "$STUB_LOG"
fail () { echo "nightly dry run: $*"; exit 1; }

# Outside an allocation: refused.
set +e; env -u SLURM_JOB_ID "$script" run > "$work/out0" 2>&1; rc=$?; set -e
[ "$rc" -eq 2 ] || fail "run outside an allocation exited $rc, expected 2"
grep -q "run inside a GPU allocation" "$work/out0" || fail "no allocation message"

# Old run logs, to be pruned to GUST_NIGHTLY_KEEP = 2 by the passing run.
mkdir -p "$work/nightly/logs/2000-01-01_000000_aaa" "$work/nightly/logs/2000-01-02_000000_bbb" \
         "$work/nightly/logs/2000-01-03_000000_ccc"

# A passing run.
export SLURM_JOB_ID=42 STUB_TESTS_JSON="$work/tests_good.json"
GUST_NIGHTLY_KEEP=2 "$script" run > "$work/out1" 2>&1 || { cat "$work/out1"; fail "passing run exited non-zero"; }
latest="$work/nightly/logs/latest"
[ "$(cat "$latest/status")" = pass ] || fail "status is not pass"
grep -q "commit     0123456789abcdef" "$latest/info" || fail "info lacks the commit"
grep -q "git clone https://github.com/erf-model/gust-abl.git $work/nightly/gust-abl" "$STUB_LOG/git.log" \
    || fail "did not clone erf-model/gust-abl"
grep -q "checkout --detach -q origin/main" "$STUB_LOG/git.log" || fail "did not check out origin/main"
grep -q "build_kestrel gpu jobs=32" "$STUB_LOG/build.log" || fail "build_kestrel.sh gpu not called with 32 jobs"
grep -q -- "-DGUST_MPIEXEC_PREFLAGS=--gpus-per-task=1;--gpu-bind=none" "$STUB_LOG/cmake.log" \
    || fail "no one-rank-per-GPU flag"
grep -q "git -C $work/nightly/gust-abl remote set-url origin https://github.com/erf-model/gust-abl.git" \
    "$STUB_LOG/git.log" || fail "origin not set to GUST_NIGHTLY_REPO"
grep -q -- "-L smoke --output-on-failure --output-junit $latest/junit.xml" "$STUB_LOG/ctest.log" \
    || grep -q -- "-L smoke --output-on-failure --output-junit $(readlink "$latest")/junit.xml" "$STUB_LOG/ctest.log" \
    || fail "ctest -L smoke with a JUnit file not run"
grep -qF -- "-E ^kestrel_dryrun$" "$STUB_LOG/ctest.log" || fail "kestrel_dryrun not left out"
grep -q "1 smoke runs of gust, 0 without the GPU flags" "$latest/flags_check.log" || fail "flag check did not count the run"
n=$(cd "$work/nightly/logs" && ls -1d 20*_* | wc -l)
[ "$n" -eq 2 ] || fail "$n run logs kept, expected 2"
[ ! -e "$work/nightly/logs/2000-01-01_000000_aaa" ] || fail "the oldest log was kept"

# A smoke run of gust without the GPU flags stops the job before the tests.
sleep 1
: > "$STUB_LOG/ctest.log"
set +e; STUB_TESTS_JSON="$work/tests_bad.json" "$script" run > "$work/out2" 2>&1; rc=$?; set -e
[ "$rc" -eq 1 ] || fail "flag check exited $rc, expected 1"
[ "$(cat "$latest/status")" = "fail: arena flags" ] || fail "status is not 'fail: arena flags'"
grep -q "missing: rt_no_flags" "$latest/flags_check.log" || fail "the run without flags is not named"
[ ! -s "$STUB_LOG/ctest.log" ] || fail "the tests ran after the flag check failed"

# Failing tests.
sleep 1
set +e; STUB_CTEST_STATUS=8 "$script" run > "$work/out3" 2>&1; rc=$?; set -e
[ "$rc" -eq 1 ] || fail "failing tests exited $rc, expected 1"
[ "$(cat "$latest/status")" = "fail: tests" ] || fail "status is not 'fail: tests'"
set +e; "$script" status > "$work/out4" 2>&1; rc=$?; set -e
[ "$rc" -eq 1 ] && grep -q "status     fail: tests" "$work/out4" || fail "status did not report the failure"

# A fetch that fails: a log with "fail: fetch", not the previous night's status;
# the logs are pruned after a failure too.
sleep 1
set +e; STUB_GIT_FAIL=fetch GUST_NIGHTLY_KEEP=2 "$script" run > "$work/out8" 2>&1; rc=$?; set -e
[ "$rc" -eq 1 ] || fail "failed fetch exited $rc, expected 1"
[ "$(cat "$latest/status")" = "fail: fetch" ] || fail "status is not 'fail: fetch'"
case "$(readlink "$latest")" in *_nofetch) ;; *) fail "latest does not point to the _nofetch log" ;; esac
n=$(cd "$work/nightly/logs" && ls -1d 20*_* | wc -l)
[ "$n" -eq 2 ] || fail "$n run logs kept after a failure, expected 2"

# A second run while one holds the lock is refused.
set +e; (flock 9; "$script" run > "$work/out9" 2>&1) 9> "$work/nightly/.lock"; rc=$?; set -e
[ "$rc" -eq 3 ] && grep -q "another run holds" "$work/out9" || fail "a concurrent run was not refused (exit $rc)"

# submit.
set +e; env -u GUST_ACCOUNT "$script" submit > "$work/out5" 2>&1; rc=$?; set -e
[ "$rc" -ne 0 ] || fail "submit without GUST_ACCOUNT succeeded"
GUST_ACCOUNT=test GUST_NIGHTLY_MAIL=a@b.c "$script" submit > "$work/out6" 2>&1 || fail "submit failed"
job="$STUB_LOG/submitted_job.sh"
for line in "#!/bin/bash -l" "#SBATCH --dependency=singleton" "#SBATCH --account=test" "#SBATCH --gpus=4" "#SBATCH --ntasks=4" "#SBATCH --mail-type=FAIL" \
            "#SBATCH --mail-user=a@b.c"; do
    grep -qx -- "$line" "$job" || fail "job script lacks '$line'"
done
grep -q "nightly_gpu.sh\" run$" "$job" || fail "job script does not call run"

# A bad GPU count.
set +e; GUST_NIGHTLY_GPUS=two "$script" run > "$work/out7" 2>&1; rc=$?; set -e
[ "$rc" -eq 2 ] && grep -q "GUST_NIGHTLY_GPUS must be a positive integer (got two)" "$work/out7" \
    || fail "a bad GUST_NIGHTLY_GPUS was not refused"
set +e; GUST_NIGHTLY_GPUS=2 "$script" run > "$work/out10" 2>&1; rc=$?; set -e
[ "$rc" -eq 2 ] && grep -q "GUST_NIGHTLY_GPUS must be 3 or 4" "$work/out10" || fail "GUST_NIGHTLY_GPUS=2 was not refused"

echo "nightly gpu dry run passed"
