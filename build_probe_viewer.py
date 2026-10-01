"""Build a self-contained probe I-V viewer from SMU buffer exports (current sourced, voltage measured).

Each CSV is one sweep. The viewer plots source current against measured voltage, plus dI/dV and
d2I/dV2 (with adjustable smoothing) for locating the plasma potential.

Usage:
    python build_probe_viewer.py <sweep.csv> [<sweep.csv> ...] [-o probe_viewer.html] [--fragment page_body.html]

`--fragment` writes the page without the <!doctype>/<head> wrapper, for hosts that add their own.
"""
import argparse
import json
import re
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent

ap = argparse.ArgumentParser()
ap.add_argument("csv", nargs="+")
ap.add_argument("-o", "--out", default="probe_viewer.html")
ap.add_argument("--fragment", help="also write the wrapper-free page to this path")
args = ap.parse_args()


def load(path):
    lines = Path(path).read_text(encoding="utf-8-sig").splitlines()
    header = next(n for n, line in enumerate(lines) if line.startswith("Index,"))
    df = pd.read_csv(path, skiprows=header, encoding="utf-8-sig")
    # Two columns are both called "Unit": pandas renames the second to "Unit.1"
    if not (df["Unit"] == "Volt DC").all() or not (df["Unit.1"] == "Amp DC").all():
        raise SystemExit(f"{path}: expected Reading in Volt DC and Value in Amp DC")
    stamp = pd.to_datetime(df["Date"] + " " + df["Time"], format="%m/%d/%Y %H:%M:%S")
    t = (stamp - stamp.iloc[0]).dt.total_seconds() + df["Fractional Seconds"] - df["Fractional Seconds"].iloc[0]
    name = re.sub(r"^[0-9a-f]{8}-", "", Path(path).name)  # drop an upload-id prefix if present
    start = stamp.iloc[0] + pd.to_timedelta(df["Fractional Seconds"].iloc[0], unit="s")
    return {
        "name": name,
        "start": start.strftime("%d %b %Y %H:%M:%S"),
        "v": df["Reading"].round(6).tolist(),
        "i": (df["Value"] * 1e6).round(6).tolist(),  # A -> uA
        "t": t.round(4).tolist(),
        "lim": (df["Source Limit"].astype(str).str.upper() == "T").astype(int).tolist(),
    }


sweeps = [load(p) for p in args.csv]
payload = json.dumps({"sweeps": sweeps}, separators=(",", ":")).replace("<", "\\u003c")

head, body = (HERE / "probe_viewer_template.html").read_text(encoding="utf-8").split("<!--BODY-->")
body = body.replace("__PROBE_DATA__", payload)

if args.fragment:
    Path(args.fragment).write_text(head + body, encoding="utf-8")

page = (
    '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
    '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
    + head + "</head>\n<body>\n" + body + "</body>\n</html>\n"
)
Path(args.out).write_text(page, encoding="utf-8")
print(f"wrote {args.out}: {len(sweeps)} sweep(s), " + ", ".join(f"{s['name']} ({len(s['v'])} readings)" for s in sweeps))
