#!/usr/bin/env python3
"""Checks of bc.inflow_profile = column (Tests/Canonical_RANS/Column_Inflow/README.md).

    check_inflow_column.py same-column TWO_RUN_TABLE TWO_RUN_LOG ONE_RUN_TABLE ONE_RUN_LOG
        The single run's column is the two-run column: the two tables (io.column_guess_file, every value in %.17g)
        are equal value for value, and the forcing both logs print (forcing.ug and forcing.vg when the drive is
        geostrophic, forcing.pressure_gradient; %.9g) is the same. The column run prints the pressure gradient only
        when it is not zero; the lines both logs have are compared, and there must be at least one.
    check_inflow_column.py record CHECKPOINT TABLE LOG Z0
        The checkpoint's InflowColumn record holds the table to the last bit, the faces' roughness Z0 (to 1e-6), and the
        forcing the log printed (to its 9 digits).
    check_inflow_column.py drop-record CHECKPOINT COPY
        Copy the checkpoint without its InflowColumn record (for the restart that must refuse it).

Exit status 0 when the check holds; 1 with a message otherwise.
"""
import os
import re
import shutil
import sys

FORCING = re.compile(r"^\s*forcing\.(ug|vg|pressure_gradient) = (.+?)\s*$")


def table(path):
    """The rows (tuples of floats) of a column table; '#' lines skipped."""
    rows = []
    with open(path) as f:
        for line in f:
            if line.strip() and not line.lstrip().startswith("#"):
                rows.append(tuple(float(t) for t in line.split()))
    if len(rows) < 2:
        sys.exit("%s: fewer than two rows" % path)
    return rows


def forcing(path):
    """The forcing a log prints (strings, as printed): ug and vg when the drive is geostrophic, pressure_gradient."""
    out = {}
    with open(path) as f:
        for line in f:
            m = FORCING.match(line)
            if m:
                out[m.group(1)] = m.group(2)
    if not out or ("ug" in out) != ("vg" in out):
        sys.exit("%s: no 'forcing.<ug|vg|pressure_gradient> = ...' lines, or only one of forcing.ug and forcing.vg" % path)
    return out


def same_rows(a, b, what):
    if len(a) != len(b):
        sys.exit("%s: %d rows against %d" % (what, len(a), len(b)))
    for i, (ra, rb) in enumerate(zip(a, b)):
        if ra != rb:
            sys.exit("%s: row %d differs: %r against %r" % (what, i, ra, rb))


def same_column(two_table, two_log, one_table, one_log):
    same_rows(table(two_table), table(one_table), "the tables %s and %s" % (two_table, one_table))
    fa, fb = forcing(two_log), forcing(one_log)
    common = sorted(set(fa) & set(fb))
    if not common:
        sys.exit("no forcing line in both %s and %s" % (two_log, one_log))
    fa = {k: fa[k] for k in common}
    fb = {k: fb[k] for k in common}
    if fa != fb:
        sys.exit("the geostrophic wind differs: %s in %s, %s in %s" % (fa, two_log, fb, one_log))
    print("same column: %d rows equal to the last bit; the same forcing in both runs: %s"
          % (len(table(one_table)), ", ".join("forcing.%s = %s" % kv for kv in sorted(fb.items()))))


def record(chk, table_path, log, z0):
    path = os.path.join(chk, "InflowColumn")
    with open(path) as f:
        tok = f.read().split()
    # GustInflowColumn 1 / z0 <z0> <map> / forcing <target> <ug> <vg> <px> <py> / coriolis <on> <lat> <period> / rows <n>
    if (tok[:2] != ["GustInflowColumn", "1"] or tok[2] != "z0" or tok[5] != "forcing" or tok[11] != "coriolis"
            or tok[15] != "rows"):
        sys.exit("%s: not a GustInflowColumn 1 record" % path)
    if abs(float(tok[3]) - float(z0)) > 1e-6 * float(z0):   # 1e-6: Z0 as a single-precision build stores it
        sys.exit("%s: z0 %s, expected %s" % (path, tok[3], z0))
    if tok[4] != "0":
        sys.exit("%s: the map flag is %s, expected 0 (a uniform bc.surface_z0)" % (path, tok[4]))
    if tok[6] != "1":
        sys.exit("%s: the forcing flag is %s, expected 1 (a target wind)" % (path, tok[6]))
    want = forcing(log)
    got = {}
    if "pressure_gradient" in want:   # z from the log: the record holds x and y
        got["pressure_gradient"] = "%.9g %.9g %.9g" % (float(tok[9]), float(tok[10]),
                                                      float(want["pressure_gradient"].split()[2]))
    if "ug" in want:
        got.update(ug="%.9g" % float(tok[7]), vg="%.9g" % float(tok[8]))
        if tok[12] != "1":   # a geostrophic wind: the column ran with Coriolis
            sys.exit("%s: the Coriolis flag is %s, expected 1" % (path, tok[12]))
    if got != want:
        sys.exit("%s: the forcing %s is not the log's %s" % (path, got, want))
    n = int(tok[16])
    vals = [float(t) for t in tok[17:]]
    if len(vals) != 5 * n:
        sys.exit("%s: %d values for %d rows" % (path, len(vals), n))
    rows = [tuple(vals[5 * i:5 * i + 5]) for i in range(n)]
    same_rows(rows, table(table_path), "the record %s and the table %s" % (path, table_path))
    print("record: %d rows equal to the table, z0 %s, the log's forcing: %s"
          % (n, tok[3], ", ".join("forcing.%s = %s" % kv for kv in sorted(got.items()))))


def drop_record(chk, copy):
    shutil.rmtree(copy, ignore_errors=True)
    shutil.copytree(chk, copy)
    os.remove(os.path.join(copy, "InflowColumn"))
    print("copied %s to %s without its InflowColumn record" % (chk, copy))


def main(argv):
    cmds = {"same-column": (same_column, 4), "record": (record, 4), "drop-record": (drop_record, 2)}
    if len(argv) < 2 or argv[1] not in cmds or len(argv) - 2 != cmds[argv[1]][1]:
        sys.exit(__doc__)
    cmds[argv[1]][0](*argv[2:])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
