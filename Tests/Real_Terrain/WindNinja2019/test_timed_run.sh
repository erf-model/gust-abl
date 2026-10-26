#!/bin/bash
# CTest wn2019_timed_run: timed_run.sh on this shell (bash 3.2 on macOS, 5 on Linux): the wall time line, the
# exit status, done only on success, arguments with spaces kept whole; and neither driver uses /usr/bin/time.
set -u
D=$(cd "$(dirname "$0")" && pwd)
. "$D/timed_run.sh"
tmp=$(mktemp -d "${TMPDIR:-/tmp}/wn2019 timed.XXXXXX")
trap 'rm -rf "$tmp"' EXIT
fail () { echo "wn2019_timed_run: $*"; exit 1; }

ok="$tmp/run ok"; mkdir -p "$ok"
timed_run "$ok" sh -c 'echo out; echo err >&2; printf "%s\n" "$1" > arg.txt' _ "terrain.file=/a b/c d.txt" \
  || fail "a successful command returned $?"
[ -f "$ok/done" ] || fail "done missing after success"
grep -qx out "$ok/log.txt" || fail "stdout not in log.txt"
grep -qx err "$ok/err.txt" || fail "stderr not in err.txt"
grep -qE '^real [0-9]+\.[0-9]+$' "$ok/time.txt" || fail "time.txt is not 'real <s>': $(cat "$ok/time.txt")"
[ "$(wc -l < "$ok/time.txt" | tr -d ' ')" = 1 ] || fail "time.txt holds more than the time line"
[ "$(cat "$ok/arg.txt")" = "terrain.file=/a b/c d.txt" ] || fail "an argument with spaces was split: $(cat "$ok/arg.txt")"

bad="$tmp/run bad"; mkdir -p "$bad"
timed_run "$bad" sh -c 'exit 3' 2> "$tmp/msg"
rc=$?
[ $rc = 3 ] || fail "a failing command returned $rc, not 3"
[ ! -f "$bad/done" ] || fail "done written after a failure"
grep -q "exit status 3" "$tmp/msg" || fail "the failure was not reported"

# sweep.sh end to end with stubs on PATH (mpirun runs its command, gust records its arguments, python3 does
# nothing), on a run directory whose path holds a space: the terrain argument reaches gust whole, the run is
# timed, and the summary line carries the wall time.
stub="$tmp/stub bin"; mkdir -p "$stub" "$tmp/build"
printf '#!/bin/sh\nshift 2\nexec "$@"\n' > "$stub/mpirun"
printf '#!/bin/sh\nprintf "%%s\\n" "$@" > args.txt\necho "gust: steady after 5 iterations"\n' > "$tmp/build/gust"
printf '#!/bin/sh\nexit 0\n' > "$stub/python3"
chmod +x "$stub/mpirun" "$stub/python3" "$tmp/build/gust"
run="$tmp/run dir"; mkdir -p "$run"; run=$(cd "$run" && pwd)
echo "v1 steady.max_iterations=5" > "$tmp/variants.txt"
out=$(PATH="$stub:$PATH" bash "$D/sweep.sh" "$tmp/build" "$run" bolund_270 "$tmp/variants.txt" 2 2>&1) \
  || fail "sweep.sh failed: $out"
v="$run/bolund_270_rans__v1"
grep -qx "terrain.file=$run/bolund.xyz" "$v/args.txt" || fail "the terrain argument was split: $(cat "$v/args.txt")"
[ -f "$v/done" ] || fail "sweep.sh: done missing"
echo "$out" | grep -qE "^v1: gust: steady after 5 iterations [0-9]+\.[0-9]+ s$" || fail "sweep.sh summary: $out"

# run_all.sh end to end with the same stubs: an mc run (no extra argument: the empty array under set -u, which
# bash 3.2 calls unbound unless spelled ${extra[@]+...}), a rans run (convection.limiter added) and an unknown
# case, skipped.
out=$(PATH="$stub:$PATH" CASES="nope bolund_270" MODES="mc rans" bash "$D/run_all.sh" "$tmp/build" "$tmp/wntools" "$run" 2 2>&1) \
  || fail "run_all.sh failed: $out"
echo "$out" | grep -q "^nope: unknown case; skipped" || fail "run_all.sh did not skip the unknown case: $out"
mc="$run/bolund_270_mc"; rans="$run/bolund_270_rans_van_leer"
[ -f "$mc/done" ] && [ -f "$rans/done" ] || fail "run_all.sh: done missing"
[ "$(cat "$mc/args.txt")" = "$(printf 'inputs\nterrain.file=%s/bolund.xyz' "$run")" ] \
  || fail "run_all.sh mc arguments: $(cat "$mc/args.txt")"
[ "$(cat "$rans/args.txt")" = "$(printf 'inputs\nterrain.file=%s/bolund.xyz\nconvection.limiter=van_leer' "$run")" ] \
  || fail "run_all.sh rans arguments: $(cat "$rans/args.txt")"
grep -qE '^real [0-9]+\.[0-9]+$' "$mc/time.txt" || fail "run_all.sh: time.txt is not 'real <s>'"

for s in run_all.sh sweep.sh; do
  bash -n "$D/$s" || fail "$s does not parse"
  if grep -q "/usr/bin/time" "$D/$s"; then fail "$s still calls /usr/bin/time"; fi
done
echo "wn2019_timed_run: passed"
