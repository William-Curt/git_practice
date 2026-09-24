# Langmuir Sweep Bench

An interactive, self-contained viewer and analysis bench for the Langmuir-probe
sweep log recorded on 2026-09-17 (`data/sweep_log.csv`, 116,729 samples,
81 half-sweeps of a ±10 V triangle bias).

**Open `langmuir_applet.html` in any modern browser.** Everything (the data, the
analysis and the charts) is embedded in that one file; it only fetches the D3
library and two fonts from a CDN.

## What the applet does

- **Animated sweep scrubber** – drag the slider, press play (½× to 4×, looping),
  use ← / → / space, or click anywhere on the whole-log strip or on a trend
  point. The I–V curve morphs smoothly from sweep to sweep and the axes rescale.
- **Plot options** – current channel (I2 / I1), overlay (single sweep, the
  up + down pair, or all 81 sweeps ghosted), linear or log |I|, measured current
  or electron current, layers (samples, smoothed line, fits and markers, lower
  panel), and a lower panel showing dI/dV or the local slope temperature Tₑ(V).
- **Baseline** – by default a straight line fitted to the plasma-off sweeps
  (75–80) is subtracted; this removes the ≈148 nA instrument offset and the
  ≈1.7 nA/V leakage ramp so the ion branch becomes visible. A constant or no
  baseline can be chosen instead.
- **Analysis on the plot** – ion-saturation fit, exponential (Tₑ) fit and its
  window, floating potential V_f with uncertainty, plasma potential V_p (or the
  sheath-theory estimate when no knee lies inside the sweep), up/down hysteresis,
  and a readout block with the fitted numbers.
- **Stat tiles and quality flags** – Tₑ (with R² and the lower/upper-half check),
  V_f, V_p, electron and ion currents, nₑ, nᵢ, λ_D, hysteresis, and explicit
  flags for anything that could not be determined (no knee, over-range event,
  no plasma signal, …).
- **Fit parameters** – smoothing window, ion-region fraction, automatic or
  manual Tₑ window, probe area and working gas; everything recomputes live.
- **Trends** – six small multiples over the whole 24-minute log (click to jump).
- **Results table** – one row per sweep, with a copy-as-CSV button.

## Rebuilding with new data

```
pip install numpy pandas
python src/build.py data/sweep_log.csv -o langmuir_applet.html
```

The CSV must have the columns `timestamp, DAC, GSE_I, V1, I1, I2, V2`
(V2 = probe bias, I1/I2 = probe current in A, electron collection positive).

## Reference analysis in Python

`analysis/langmuir_analysis.py` is a numpy-only implementation of exactly the
same algorithm (spec in `analysis/ANALYSIS_SPEC_v2.md`, amendments at the end);
it gives a scriptable way to export per-sweep results and is used to
cross-check the JavaScript port:

```
python analysis/langmuir_analysis.py data/sweep_log.csv --csv-out results.csv
python analysis/langmuir_analysis.py data/sweep_log.csv --baseline none --te-manual -4,2 --out results.json
```

`analysis/results_default.csv` is the output for the applet's default settings
(channel I2, plasma-off baseline from sweeps 75–80, window 31, ion region 20 %,
Tₑ window 2–30 %, area 10 mm², argon).

To verify that the JavaScript in the applet reproduces the Python numbers
(they agree to machine precision on this log):

```
python analysis/langmuir_analysis.py data/sweep_log.csv --quantize-float32 --out ref.json --dump-segments 1,10,30
node analysis/compare_js_py.mjs ref.json          # or: ... ref.json none / const
```

`--quantize-float32` rounds the inputs the way the applet stores them.

## Layout

| Path | Contents |
|---|---|
| `langmuir_applet.html` | the built applet (open this) |
| `data/sweep_log.csv` | the raw sweep log |
| `src/template.html`, `src/app.js`, `src/analysis.js` | page, UI/charts, analysis (JS) |
| `src/build.py` | packs the CSV and assembles the applet |
| `analysis/langmuir_analysis.py` | Python reference implementation + CLI |
| `analysis/ANALYSIS_SPEC_v2.md` | the algorithm, step by step |
| `analysis/compare_js_py.mjs` | checks the JS port against the Python output |
