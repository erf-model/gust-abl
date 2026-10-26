# Sourced by run_all.sh and sweep.sh.
# timed_run DIR CMD [ARG...]: runs CMD in DIR, stdout to DIR/log.txt, stderr to DIR/err.txt, and the wall
# time as one line "real <seconds>" to DIR/time.txt. It uses bash's time keyword, which behaves the same on
# Linux and macOS (/usr/bin/time -l is BSD only; GNU time refuses -l). Each argument is passed as given, so a
# path with spaces stays one argument. On success DIR/done is touched; otherwise the exit status is printed
# and returned.
timed_run () {
  local dir=$1
  shift
  ( cd "$dir" || exit 1
    TIMEFORMAT='real %R'
    { time "$@" < /dev/null > log.txt 2> err.txt ; } 2> time.txt )
  local rc=$?
  if [ $rc -eq 0 ]; then
    touch "$dir/done"
  else
    echo "$dir: exit status $rc (see $dir/err.txt)" >&2
  fi
  return $rc
}
