"""Probe current-sweep plot: voltage (y, linear) vs. current (x, log)."""
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import FuncFormatter

HERE = Path(__file__).parent

SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e6e5e1"
BLUE, ORANGE = "#2a78d6", "#eb6834"

FILES = [
    # (csv, legend label, colour); the coarse stepped sweep is drawn first (behind)
    ("1.5u_66u_0.5uStep.csv", "1.5u_66u_0.5uStep", ORANGE),
    ("defbuffer1.csv", "defbuffer1", BLUE),
]


def load(name):
    # Keithley buffer export: 8 header lines, then a column row.
    # Col 1 = measured voltage (V); col 14 = sourced current (A).
    df = pd.read_csv(HERE / "data" / name, skiprows=8)
    return pd.DataFrame({"V": df.iloc[:, 1], "I_uA": df.iloc[:, 14] * 1e6})


fig, ax = plt.subplots(figsize=(8, 5), dpi=200, facecolor=SURFACE)
ax.set_facecolor(SURFACE)

dropped = {}
for z, (name, label, colour) in enumerate(FILES):
    df = load(name)
    pos = df[df.I_uA > 0]  # a log axis can only show I > 0
    dropped[label] = len(df) - len(pos)
    dense = len(pos) > 1000
    # Dense sweep -> line. Sparse stepped sweep -> markers only, so no line is
    # drawn across the gap between the ~0 uA setpoint and the first real step.
    ax.plot(
        pos.I_uA, pos.V,
        color=colour, zorder=2 + z,
        lw=1.4 if dense else 0, marker=None if dense else "o", markersize=4,
        markeredgecolor=SURFACE, markeredgewidth=0.5,
        label=f"{label}  ({len(pos):,} pts)",
    )

ax.set_xscale("log")
ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
ax.set_xlabel("Current (µA)", color=INK_2)
ax.set_ylabel("Voltage (V)", color=INK_2)
ax.set_title("Probe current sweep: voltage vs. current", loc="left",
             color=INK, fontsize=13, fontweight="bold", pad=12)

ax.grid(True, which="major", color=GRID, lw=0.8)
ax.grid(True, which="minor", axis="x", color=GRID, lw=0.4, alpha=0.7)
ax.set_axisbelow(True)
ax.tick_params(colors=INK_2, length=0)
for side in ("top", "right", "left"):
    ax.spines[side].set_visible(False)
ax.spines["bottom"].set_color(GRID)

leg = ax.legend(loc="lower right", frameon=False, labelcolor=INK)

note = ("Log axis shows I > 0 only; omitted non-positive points: "
        + ", ".join(f"{k} {v:,}" for k, v in dropped.items()))
fig.text(0.015, 0.012, note, color=INK_2, fontsize=7.5, ha="left", va="bottom")

fig.tight_layout(rect=(0, 0.03, 1, 1))
out = HERE / "iv_plot.png"
fig.savefig(out, facecolor=SURFACE)
print("wrote", out, dropped)
