"""Build a self-contained, interactive sweep viewer from a sweep log CSV.

Each sweep is one ramp of the DAC between turning points (0 <-> 65535). The viewer
shows I2 (uA) against V2 (V) for one sweep at a time, with axes fitted to that sweep.

Usage:
    python build_sweep_viewer.py <sweep_log.csv> [-o sweep_viewer.html] [--fragment page_body.html]

`--fragment` writes the page without the <!doctype>/<head> wrapper, for hosts that add their own.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent

ap = argparse.ArgumentParser()
ap.add_argument("csv")
ap.add_argument("-o", "--out", default="sweep_viewer.html")
ap.add_argument("--fragment", help="also write the wrapper-free page to this path")
args = ap.parse_args()

df = pd.read_csv(args.csv, parse_dates=["timestamp"])
t_ms = ((df["timestamp"] - df["timestamp"].iloc[0]).dt.total_seconds() * 1000).round().astype(int).to_numpy()
dac = df["DAC"].astype(int).to_numpy()
v2 = df["V2"].round(4).to_numpy()
i2_ua = (df["I2"] * 1e6).round(5).to_numpy()

# Turning points are the rows where the DAC changes direction. Each sweep runs from one
# turning point to the next (both ends included), so every full sweep spans 0 <-> 65535.
step = np.sign(np.diff(dac))
turns = np.flatnonzero(step[1:] * step[:-1] < 0) + 1
edges = np.concatenate([[0], turns, [len(df) - 1]])

sweeps = []
for k, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
    sl = slice(a, b + 1)
    sweeps.append({
        "partial": bool(k == 0 or k == len(edges) - 2),
        "t": t_ms[sl].tolist(),
        "dac": dac[sl].tolist(),
        "v": v2[sl].tolist(),
        "i": i2_ua[sl].tolist(),
    })

data = {"t0": df["timestamp"].iloc[0].strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3], "rows": len(df), "sweeps": sweeps}
payload = json.dumps(data, separators=(",", ":"))

head, body = (HERE / "sweep_viewer_template.html").read_text(encoding="utf-8").split("<!--BODY-->")
body = body.replace("__SWEEP_DATA__", payload)

if args.fragment:
    Path(args.fragment).write_text(head + body, encoding="utf-8")

page = (
    '<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
    '<meta name="viewport" content="width=device-width, initial-scale=1">\n'
    + head + "</head>\n<body>\n" + body + "</body>\n</html>\n"
)
Path(args.out).write_text(page, encoding="utf-8")
print(f"wrote {args.out}: {len(sweeps)} sweeps, {len(df):,} rows, {len(page) / 1e6:.1f} MB")
