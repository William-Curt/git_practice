#!/usr/bin/env python3
"""Build the Langmuir Sweep Bench applet from a sweep CSV.

    python build.py data/sweep_log.csv -o langmuir_applet.html

Reads the CSV (timestamp, DAC, GSE_I, V1, I1, I2, V2), packs the columns into
base64 typed arrays, and inlines them into the page together with the analysis
and UI scripts. The output is a single self-contained HTML file.
"""
import argparse, base64, json, os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D3_CDN = "https://cdnjs.cloudflare.com/ajax/libs/d3/7.9.0/d3.min.js"


def segment_bounds(dac):
    d = np.sign(np.diff(dac.astype(np.float64)))
    turns = np.where(np.diff(d) != 0)[0] + 1
    return [0] + turns.tolist() + [len(dac)]


def pack(csv_path):
    df = pd.read_csv(csv_path)
    for c in ("timestamp", "DAC", "GSE_I", "I1", "I2", "V2"):
        if c not in df.columns:
            sys.exit(f"column {c} missing from {csv_path}")
    t = pd.to_datetime(df["timestamp"])
    ts = (t - t.iloc[0]).dt.total_seconds().to_numpy(np.float32)
    dac = np.round(df["DAC"].to_numpy(np.float64)).astype(np.uint16)
    gse = np.round((df["GSE_I"].to_numpy(np.float64) - 0.05) * 1e6).clip(0, 65535).astype(np.uint16)
    b64 = lambda a: base64.b64encode(np.ascontiguousarray(a).tobytes()).decode("ascii")
    payload = {
        "meta": {
            "source": os.path.basename(csv_path),
            "n": int(len(df)),
            "t0_iso": t.iloc[0].isoformat(timespec="milliseconds"),
            "columns": ["ts", "V2", "I1", "I2", "DAC", "GSE_I"],
            "encodings": {
                "ts": {"dtype": "float32", "unit": "s since t0"},
                "V2": {"dtype": "float32", "unit": "V"},
                "I1": {"dtype": "float32", "unit": "A", "nan": "preserved"},
                "I2": {"dtype": "float32", "unit": "A"},
                "DAC": {"dtype": "uint16"},
                "GSE_I": {"dtype": "uint16", "scale": 1e-6, "offset": 0.05, "unit": "A"},
            },
        },
        "bounds": segment_bounds(dac),
        "cols": {
            "ts": b64(ts),
            "V2": b64(df["V2"].to_numpy(np.float32)),
            "I1": b64(df["I1"].to_numpy(np.float32)),
            "I2": b64(df["I2"].to_numpy(np.float32)),
            "DAC": b64(dac),
            "GSE_I": b64(gse),
        },
    }
    return payload


def build(csv_path, out_path, d3_src=D3_CDN, fragment=False, template=None):
    payload = pack(csv_path)
    read = lambda name: open(os.path.join(HERE, name), encoding="utf-8").read()
    page = template if template is not None else read("template.html")
    analysis = read("analysis.js")
    app = read("app.js")
    data_js = "window.SWEEP_DATA = " + json.dumps(payload, separators=(",", ":")) + ";"
    scripts = (
        f'<script src="{d3_src}"></script>\n'
        f"<script>{data_js}</script>\n"
        f"<script>{analysis}</script>\n"
        f"<script>{app}</script>\n"
    )
    body = page + "\n" + scripts
    if fragment:
        html = body
    else:
        html = (
            "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
            "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1, viewport-fit=cover\">\n"
            + body.replace("<title>", "<title>", 1)
        )
        # move <title>, <meta>, <link>, <style> into <head>; everything from the first <div into <body>
        cut = html.index('<div class="app"')
        head, rest = html[:cut], html[cut:]
        html = head + "</head>\n<body>\n" + rest + "</body>\n</html>\n"
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)
    n = payload["meta"]["n"]
    return {"rows": n, "segments": len(payload["bounds"]) - 1, "bytes": len(html.encode("utf-8"))}


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv")
    ap.add_argument("-o", "--out", default="langmuir_applet.html")
    ap.add_argument("--d3", default=D3_CDN, help="script src for d3 (default: cdnjs 7.9.0)")
    ap.add_argument("--fragment", action="store_true", help="emit a body fragment (no <html>/<head>/<body>) for artifact publishing")
    a = ap.parse_args()
    info = build(a.csv, a.out, a.d3, a.fragment)
    print(f"wrote {a.out}: {info['rows']} rows, {info['segments']} half-sweeps, {info['bytes']/1e6:.2f} MB")
