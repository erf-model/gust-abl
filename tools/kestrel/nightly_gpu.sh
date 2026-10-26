#!/usr/bin/env bash
# Nightly GPU job on NLR Kestrel (H100): build main for CUDA and run every
# smoke test on the GPUs. The smoke tests include the 2-D hill restart chain
# (rt_steady_hill2d_restart_*, audit A1) and the Coriolis case run twice
# (rt_coriolis_repeat_*, audit A4). Every smoke run of gust carries
# amrex.the_arena_init_size=0 amrex.init_snan=1 amrex.fpe_trap_invalid=1
# (GUST_TEST_AMREX_FLAGS in Tests/CMakeLists.txt); the job checks that before
# it runs them.
#
#   nightly_gpu.sh run              build and test here (inside a GPU allocation)
#   nightly_gpu.sh submit           sbatch the same as a batch job
#   nightly_gpu.sh status           print the result of the latest run
#
# Environment:
#   GUST_ACCOUNT        Kestrel allocation (required for submit)
#   GUST_NIGHTLY_WORK   work directory (default /scratch/$USER/gust-nightly)
#   GUST_NIGHTLY_REPO   repository to test (default https://github.com/erf-model/gust-abl.git)
#   GUST_NIGHTLY_REF    branch or commit (default main)
#   GUST_NIGHTLY_GPUS   GPUs to request (default 4; the smoke tests use up to 3 ranks)
#   GUST_NIGHTLY_MAIL   mail address for a failed job (submit only; optional)
#   GUST_NIGHTLY_KEEP   number of run logs kept (default 30)
#
# The job clones the repository into its work directory and builds with the
# build scripts of that checkout; the copy this script is started from is not
# changed. The build
# directory is kept between nights (an incremental build after the first).
# Each run writes $work/logs/<date>_<sha>/ with build.log, ctest.log,
# junit.xml and status ("pass" or "fail: <step>"); $work/logs/latest points
# to it. Exit status 0 only when every smoke test passed.
#
# One rank per GPU: the tests launch srun with --gpus-per-task=1
# (GUST_MPIEXEC_PREFLAGS), one test at a time.
#
# Scheduling: scrontab (where Kestrel enables it), or cron on a login node:
#   0 1 * * * GUST_ACCOUNT=<account> /path/to/nightly_gpu.sh submit
#
# Partition and GPU request syntax follow Kestrel's documentation as of
# October 2026 and may need adjusting after the system update.
set -euo pipefail

work=${GUST_NIGHTLY_WORK:-/scratch/$USER/gust-nightly}
repo=${GUST_NIGHTLY_REPO:-https://github.com/erf-model/gust-abl.git}
ref=${GUST_NIGHTLY_REF:-main}
gpus=${GUST_NIGHTLY_GPUS:-4}
keep=${GUST_NIGHTLY_KEEP:-30}
src="$work/gust-abl"
self=$(cd "$(dirname "$0")" && pwd)/$(basename "$0")

usage () { sed -n '2,38p' "$0" >&2; exit 2; }

case "$gpus" in ''|*[!0-9]*|0) echo "GUST_NIGHTLY_GPUS must be a positive integer (got $gpus)" >&2; exit 2 ;; esac
case "$keep" in ''|*[!0-9]*|0) echo "GUST_NIGHTLY_KEEP must be a positive integer (got $keep)" >&2; exit 2 ;; esac

stamp () { date '+%Y-%m-%dT%H:%M:%S%z'; }

# Writes the status file and exits: fail <step> (exit 1) or pass (exit 0).
finish () {
    echo "$1" > "$log/status"
    echo "nightly_gpu: $1 ($log)"
    if [ "$1" = pass ]; then exit 0; fi
    exit 1
}

cmd=${1:-}
case "$cmd" in

run)
    if [ -z "${SLURM_JOB_ID:-}" ]; then
        echo "run inside a GPU allocation, e.g. salloc -A \$GUST_ACCOUNT --gpus=$gpus -t 03:00:00" >&2
        exit 2
    fi
    mkdir -p "$work/logs"
    if [ ! -d "$src/.git" ]; then
        git clone "$repo" "$src"
    fi
    git -C "$src" fetch --prune origin
    git -C "$src" checkout --detach -q "origin/$ref" 2>/dev/null || git -C "$src" checkout --detach -q "$ref"
    sha=$(git -C "$src" rev-parse --short=9 HEAD)
    log="$work/logs/$(date +%Y-%m-%d_%H%M%S)_$sha"
    mkdir -p "$log"
    ln -sfn "$log" "$work/logs/latest"
    {
        echo "repository $repo"
        echo "ref        $ref"
        echo "commit     $(git -C "$src" rev-parse HEAD)"
        echo "job        $SLURM_JOB_ID on $(hostname)"
        echo "started    $(stamp)"
    } > "$log/info"

    build="$src/build-kestrel-gpu"
    GUST_BUILD_JOBS=${GUST_BUILD_JOBS:-${SLURM_CPUS_ON_NODE:-16}} \
        "$src/tools/kestrel/build_kestrel.sh" gpu > "$log/build.log" 2>&1 || finish "fail: build"
    # One rank per GPU for every test.
    # shellcheck source=kestrel_env.sh
    source "$src/tools/kestrel/kestrel_env.sh" gpu >> "$log/build.log" 2>&1
    cmake "$build" -DGUST_MPIEXEC_PREFLAGS="--gpus-per-task=1" >> "$log/build.log" 2>&1 \
        || finish "fail: configure"
    cmake --build "$build" -j "${GUST_BUILD_JOBS:-${SLURM_CPUS_ON_NODE:-16}}" >> "$log/build.log" 2>&1 \
        || finish "fail: build"

    # Every smoke run of gust must carry the GPU runtime flags.
    ctest --test-dir "$build" -L smoke --show-only=json-v1 > "$log/tests.json"
    python3 -I - "$log/tests.json" > "$log/flags_check.log" 2>&1 <<'PY' || finish "fail: arena flags"
import json, os, sys
tests = json.load(open(sys.argv[1]))["tests"]
need = ("amrex.the_arena_init_size=0", "amrex.init_snan=1", "amrex.fpe_trap_invalid=1")
runs = [t for t in tests if any(os.path.basename(a) == "gust" for a in t.get("command", []))]
bad = [t["name"] for t in runs if not all(f in t["command"] for f in need)]
print(f"{len(runs)} smoke runs of gust, {len(bad)} without the GPU flags")
for name in bad:
    print("  missing:", name)
sys.exit(1 if bad or not runs else 0)
PY

    echo "ctest      $(stamp)" >> "$log/info"
    ctest --test-dir "$build" -L smoke --output-on-failure --output-junit "$log/junit.xml" \
        > "$log/ctest.log" 2>&1 || { tail -n 40 "$log/ctest.log"; finish "fail: tests"; }
    echo "finished   $(stamp)" >> "$log/info"

    # Keep the newest $keep run logs.
    (cd "$work/logs" && ls -1d 20*_* 2>/dev/null | sort |
        awk -v k="$keep" '{ d[NR] = $0 } END { for (i = 1; i <= NR - k; i++) print d[i] }' |
        while read -r d; do rm -rf "./$d"; done)
    finish pass
    ;;

submit)
    : "${GUST_ACCOUNT:?set GUST_ACCOUNT to the Kestrel allocation}"
    mkdir -p "$work/jobs"
    job="$work/jobs/nightly_gpu_$(date +%Y-%m-%d_%H%M).sh"
    mail=""
    if [ -n "${GUST_NIGHTLY_MAIL:-}" ]; then
        mail="#SBATCH --mail-type=FAIL
#SBATCH --mail-user=$GUST_NIGHTLY_MAIL"
    fi
    cat > "$job" <<JOB
#!/usr/bin/env bash
#SBATCH --account=$GUST_ACCOUNT
#SBATCH --job-name=gust-nightly-gpu
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --gpus=$gpus
#SBATCH --ntasks=$gpus
#SBATCH --cpus-per-task=8
#SBATCH --mem=240G
#SBATCH --output=$work/jobs/%x-%j.log
$mail
set -euo pipefail
export GUST_NIGHTLY_WORK="$work" GUST_NIGHTLY_REPO="$repo" GUST_NIGHTLY_REF="$ref"
export GUST_NIGHTLY_GPUS="$gpus" GUST_NIGHTLY_KEEP="$keep"
"$self" run
JOB
    (cd "$work/jobs" && sbatch "$job")
    ;;

status)
    if [ ! -e "$work/logs/latest/status" ]; then
        echo "no finished run in $work/logs" >&2
        exit 1
    fi
    cat "$work/logs/latest/info"
    echo "status     $(cat "$work/logs/latest/status")"
    [ "$(cat "$work/logs/latest/status")" = pass ]
    ;;

*) usage ;;
esac
