#!/usr/bin/env python3
"""Download what the case needs into a folder (default: beside this script): the terrain awaken.asc (AWS Terrain
Tiles, 36.05-36.68 N, 97.80-97.01 W, as ../HRRR_Outer) and the HRRR crops of the C1 and C2 hours (byte ranges of HRRR's
wrfnat and wrfsfc files on AWS, through Herbie; 45 km half-width, which also covers the 40 km box of the theta_v
gradient in tools/hrrr_regime.py column).

    python3 prepare.py [--out DIR] [--save GRIB_CACHE] [--hours may10_20 ...] [--neighbours] [--columns DIR]

--neighbours also downloads the hour before and the hour after each hour, which the classifier needs to test the
transitions (C7, C8) and the tendency (T); --columns writes the column of each crop for hrrr_regime.py
catalogue. Exit status 0 when every file is written; files already there are kept.
"""
import argparse
import datetime as dt
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.normpath(os.path.join(HERE, "..", "..", "..", "tools"))
LAT, LON = 36.362, -97.405
HALF_KM = 45
# tag: valid time (UTC) of the C1 (convective) and C2 (weakly unstable) hours; each is forecast hour 1 of the
# cycle an hour before
HOURS = {"may10_19": "2023-05-10T19:00Z", "may10_20": "2023-05-10T20:00Z", "may10_21": "2023-05-10T21:00Z",
         "apr02_18": "2023-04-02T18:00Z", "apr04_18": "2023-04-04T18:00Z", "apr19_21": "2023-04-19T21:00Z"}
FXX = 1


def valid_time(tag):
    return dt.datetime.strptime(HOURS[tag], "%Y-%m-%dT%H:%MZ")


def crop_name_at(valid):
    return f"awaken_{valid - dt.timedelta(hours=FXX):%Y%m%d%H}_f{FXX:02d}.npz"


def crop_name(tag):
    return crop_name_at(valid_time(tag))


def fetch(valid, out_dir, save):
    out = os.path.join(out_dir, crop_name_at(valid))
    if not os.path.exists(out):
        cycle = valid - dt.timedelta(hours=FXX)
        subprocess.run([sys.executable, os.path.join(TOOLS, "hrrr_outer.py"), "fetch", "--cycle",
                        f"{cycle:%Y-%m-%dT%H:%MZ}", "--fxx", str(FXX), "--lat", str(LAT), "--lon", str(LON),
                        "--half-km", str(HALF_KM), "--save", save, "--out", out], check=True)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=HERE)
    ap.add_argument("--save", default=os.path.join(os.path.expanduser("~"), ".cache", "gust", "hrrr"))
    ap.add_argument("--hours", nargs="+", default=list(HOURS), choices=list(HOURS))
    ap.add_argument("--neighbours", action="store_true", help="also the hours before and after each hour")
    ap.add_argument("--columns", default=None, metavar="DIR",
                    help="also write the column of each crop (tools/hrrr_regime.py column) into DIR")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    dem = os.path.join(a.out, "awaken.asc")
    if not os.path.exists(dem):
        subprocess.run([sys.executable, os.path.join(TOOLS, "fetch_terrain.py"), "--box", "36.05", "36.68", "-97.80",
                        "-97.01", "--out", dem], check=True)
    for tag in a.hours:
        v = valid_time(tag)
        for d in ((-1, 0, 1) if a.neighbours else (0,)):
            valid = v + dt.timedelta(hours=d)
            crop = fetch(valid, a.out, a.save)
            print(f"{tag} {d:+d} h: {crop}")
            if a.columns:
                os.makedirs(a.columns, exist_ok=True)
                col = os.path.join(a.columns, f"hrrr_{valid:%Y%m%d%H}_f{FXX:02d}.json")
                subprocess.run([sys.executable, os.path.join(TOOLS, "hrrr_regime.py"), "column", crop, "--valid",
                                f"{valid:%Y-%m-%dT%H:%MZ}", "--lat", str(LAT), "--lon", str(LON), "--fxx", str(FXX),
                                "--out", col], check=True)


if __name__ == "__main__":
    main()
