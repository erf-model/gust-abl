#!/usr/bin/env python3
"""Download what the case needs into a folder (default: beside this script): the terrain awaken.asc (AWS Terrain
Tiles, 36.05-36.68 N, 97.80-97.01 W: the raw HRRR domain of 56 km with its margin) and the HRRR crops of the five
hours (byte ranges of HRRR's wrfnat and wrfsfc files on AWS, through Herbie; 32 km half-width).

    python3 prepare.py [--out DIR] [--save GRIB_CACHE] [--hours day20 night09 ...]

Exit status 0 when every file is written; files already there are kept.
"""
import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.normpath(os.path.join(HERE, "..", "..", "..", "tools"))
LAT, LON = 36.362, -97.405
# tag: (cycle, forecast hour); the valid time is the cycle plus the forecast hour
HOURS = {"day20": ("2023-05-10T19:00Z", 1), "eve01": ("2023-05-11T00:00Z", 1), "jet05": ("2023-05-10T04:00Z", 1),
         "night09": ("2023-05-11T08:00Z", 1), "dawn12": ("2023-05-10T11:00Z", 1)}


def crop_name(tag):
    cycle, fxx = HOURS[tag]
    return f"awaken_{cycle[0:4]}{cycle[5:7]}{cycle[8:10]}{cycle[11:13]}_f{fxx:02d}.npz"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=HERE)
    ap.add_argument("--save", default=os.path.join(os.path.expanduser("~"), ".cache", "gust", "hrrr"))
    ap.add_argument("--hours", nargs="+", default=list(HOURS), choices=list(HOURS))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    dem = os.path.join(a.out, "awaken.asc")
    if not os.path.exists(dem):
        subprocess.run([sys.executable, os.path.join(TOOLS, "fetch_terrain.py"), "--box", "36.05", "36.68", "-97.80",
                        "-97.01", "--out", dem], check=True)
    for tag in a.hours:
        out = os.path.join(a.out, crop_name(tag))
        if os.path.exists(out):
            continue
        cycle, fxx = HOURS[tag]
        subprocess.run([sys.executable, os.path.join(TOOLS, "hrrr_outer.py"), "fetch", "--cycle", cycle, "--fxx",
                        str(fxx), "--lat", str(LAT), "--lon", str(LON), "--half-km", "32", "--save", a.save,
                        "--out", out], check=True)
        print(f"{tag}: {out}")


if __name__ == "__main__":
    main()
