#!/usr/bin/env python3
"""The bc.inflow_file table of a Nests_Flat_ABL deck, from that deck's steady column.

    make_column_table.py DECK --column-deck COLUMN_DECK
    gust COLUMN_DECK
    make_column_table.py DECK --table column_guess.dat --out TABLE

DECK is one of the case's decks (inputs_neutral_s, inputs_buoyant_b, inputs_stable_gabls1). --column-deck writes
a copy of DECK without its nest, its inflow table and its layer on u and v (bc.rayleigh_fields = w theta), which
stops after one iteration: its run writes the steady column's table (io.column_guess_file, column_guess.dat).

TABLE keeps the column's u, v and k and replaces its theta by the deck's sounding (bc.theta_profile_z,
bc.theta_profile, and bc.zhi_theta_grad above the last height). The column takes its start theta and its
theta target from the inflow table, so a table with the sounding gives the column the same theta it had without
the table: the column, with u and v damped towards the table, then returns the same solution, and the layer on
u and v only removes departures from it.
"""
import argparse
import os
import re
import sys


def deck_values(path):
    vals = {}
    with open(path) as fh:
        for line in fh:
            line = line.split("#", 1)[0].strip()
            m = re.match(r"([\w.]+)\s*=\s*(.*)", line)
            if m:
                vals[m.group(1)] = m.group(2).split()
    return vals


def sounding(z, zs, ts, grad):
    if z >= zs[-1]:
        return ts[-1] + grad * (z - zs[-1])
    for i in range(len(zs) - 1):
        if zs[i] <= z <= zs[i + 1]:
            return ts[i] + (ts[i + 1] - ts[i]) * (z - zs[i]) / (zs[i + 1] - zs[i])
    return ts[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("deck")
    ap.add_argument("--column-deck", help="write the column deck to this file and stop")
    ap.add_argument("--table", help="the column's table (column_guess.dat)")
    ap.add_argument("--out", help="the table to write")
    a = ap.parse_args()
    if a.column_deck:
        drop = ("bc.rayleigh_fields", "bc.inflow_profile", "bc.inflow_file", "nest.", "steady.max_iterations",
                "io.write_checkpoint")
        with open(a.deck) as fh, open(a.column_deck, "w") as out:
            out.write(f"# The column deck of {os.path.basename(a.deck)} (make_column_table.py --column-deck)\n")
            for line in fh:
                if not line.lstrip().startswith(drop) and not line.startswith("#"):
                    out.write(line)
            out.write("bc.rayleigh_fields = w theta\nsteady.max_iterations = 1\n")
        print(f"wrote {a.column_deck}")
        return
    if not (a.table and a.out):
        sys.exit("give --column-deck, or --table and --out")
    d = deck_values(a.deck)
    zs = [float(x) for x in d["bc.theta_profile_z"]]
    ts = [float(x) for x in d["bc.theta_profile"]]
    grad = float(d.get("bc.zhi_theta_grad", ["0"])[0])
    header, rows = [], []
    with open(a.table) as fh:
        for line in fh:
            if line.startswith("#"):
                header.append(line.rstrip("\n"))
            elif line.strip():
                rows.append([float(x) for x in line.split()])
    if not rows or any(len(r) != 5 for r in rows):
        sys.exit(f"{a.table}: expected rows of z u v theta k")
    with open(a.out, "w") as fh:
        for h in header:
            fh.write(h + "\n")
        fh.write(f"# theta replaced by the sounding of {os.path.basename(a.deck)} (make_column_table.py)\n")
        for z, u, v, _, k in rows:
            fh.write(f"{z:.17g} {u:.17g} {v:.17g} {sounding(z, zs, ts, grad):.17g} {k:.17g}\n")
    print(f"wrote {a.out}: {len(rows)} rows")


if __name__ == "__main__":
    main()
