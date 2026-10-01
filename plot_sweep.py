"""Plot a sweep log: current (I2) vs. voltage (V2), one line per ramp, colored by time.

Usage: python plot_sweep.py <sweep_log.csv> [output.png]
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable

SURFACE, INK, INK_2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, BASELINE = "#e1e0d9", "#c3c2b7"
# One-hue sequential ramp (blue 250 -> 700), light = early, dark = late
RAMP = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"]

src = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else "sweep_plot.png"

df = pd.read_csv(src, parse_dates=["timestamp"])
t_min = (df["timestamp"] - df["timestamp"].iloc[0]).dt.total_seconds().to_numpy() / 60
v = df["V2"].to_numpy()
i_ua = df["I2"].to_numpy() * 1e6  # A -> uA

# Split into ramps at the DAC turning points (each ramp is one up or down pass)
step = np.sign(np.diff(df["DAC"].to_numpy()))
turns = np.flatnonzero(step[1:] * step[:-1] < 0) + 1
bounds = np.concatenate([[0], turns + 1, [len(df)]])
segments = [np.column_stack([v[a:b], i_ua[a:b]]) for a, b in zip(bounds[:-1], bounds[1:]) if b - a > 1]
seg_t = np.array([t_min[a:b].mean() for a, b in zip(bounds[:-1], bounds[1:]) if b - a > 1])
n_cycles = len(turns) // 2

cmap = LinearSegmentedColormap.from_list("blue_seq", RAMP)
norm = Normalize(t_min.min(), t_min.max())

fig, ax = plt.subplots(figsize=(10, 6), dpi=160, facecolor=SURFACE)
ax.set_facecolor(SURFACE)
ax.add_collection(LineCollection(segments, colors=cmap(norm(seg_t)), linewidths=0.9, alpha=0.8,
                                 capstyle="round", joinstyle="round"))
ax.set_xlim(-10.5, 10.5)
ax.set_ylim(i_ua.min() - 0.15, i_ua.max() + 0.25)

ax.axhline(0, color=BASELINE, lw=1, zorder=0)
ax.grid(True, color=GRID, lw=0.8)
ax.set_axisbelow(True)
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
for side in ("left", "bottom"):
    ax.spines[side].set_color(BASELINE)
ax.tick_params(colors=MUTED, labelcolor=INK_2, length=0, labelsize=10)

ax.set_xlabel("Swept voltage, V2 (V)", color=INK_2, fontsize=11, labelpad=8)
ax.set_ylabel("Measured current, I2 (µA)", color=INK_2, fontsize=11, labelpad=8)

# Call out the ramps that dip below zero, the one thing that stands out from the rest
neg = i_ua < 0
if neg.any():
    ax.text(6.6, i_ua.min() + 0.35, "Dips below zero at V2 ≈ 7–10 V\n(cycles at ~2.5 and ~17 min)",
            color=INK_2, fontsize=9.5, ha="right", va="center")

cb = fig.colorbar(ScalarMappable(norm, cmap), ax=ax, pad=0.02, fraction=0.04)
cb.set_label("Elapsed time (min)", color=INK_2, fontsize=11, labelpad=8)
cb.outline.set_visible(False)
cb.ax.tick_params(colors=MUTED, labelcolor=INK_2, length=0, labelsize=10)

start, end = df["timestamp"].iloc[0], df["timestamp"].iloc[-1]
fig.text(0.07, 0.955, "Current vs. voltage across repeated sweeps", color=INK, fontsize=15,
         fontweight="bold", ha="left", va="top")
fig.text(0.07, 0.905,
         f"{n_cycles} up/down cycles ({len(df):,} readings), {start:%d %b %Y %H:%M}–{end:%H:%M}. "
         "Each line is one ramp.",
         color=INK_2, fontsize=10.5, ha="left", va="top")

fig.subplots_adjust(left=0.09, right=0.88, top=0.86, bottom=0.11)
fig.savefig(out, facecolor=SURFACE)
print(f"saved {out}: {len(segments)} ramps, {n_cycles} cycles")
