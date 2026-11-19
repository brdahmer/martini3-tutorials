#!/usr/bin/env python3
"""Plot KALP helix tilt over time for the Insane, COBY and self-assembly routes in POPC,
and for the Insane route in DLPC (section 5.1).

gmx gangle reports the angle between the residue 3 -> 21 backbone vector and +z.
The vector points down, so the tilt from the membrane normal is 180 - angle.

Run from this folder:  python3 plot_tilt.py
Writes ../figures/tilt.png
"""
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Colors from the reference dataviz palette (light mode), categorical slots 1-4.
SURFACE, INK, INK_2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
ROUTES = [
    ("POPC · Insane (water bias)", "insane-angles.xvg", "#2a78d6"),
    ("POPC · COBY", "coby-angles.xvg", "#eb6834"),
    ("POPC · self-assembly (water bias)", "self-assembly-angles.xvg", "#1baf7a"),
    ("DLPC · Insane (water bias)", "dlpc-angles.xvg", "#eda100"),
]


def load(path):
    d = np.loadtxt(path, comments=("#", "@"))
    return d[:, 0] / 1000.0, 180.0 - d[:, 1]  # ps -> ns, angle -> tilt


plt.rcParams.update({"font.size": 10, "font.family": "sans-serif"})
fig, axes = plt.subplots(len(ROUTES), 1, figsize=(7.0, 8.4), dpi=200, sharex=True, sharey=True)
fig.patch.set_facecolor(SURFACE)

for ax, (name, path, color) in zip(axes, ROUTES):
    t, tilt = load(path)
    ax.set_facecolor(SURFACE)
    # The self-assembly run saved one frame per ns, so mark its points.
    ax.plot(t, tilt, color=color, linewidth=1.0, marker="o" if len(t) < 100 else None, markersize=3)
    mean, sd = tilt.mean(), tilt.std()
    ax.axhline(mean, color=INK_2, linewidth=1.0, linestyle="--")
    # Name the panel inside the plot, and label the mean line directly.
    ax.text(0.01, 0.92, name, transform=ax.transAxes, color=INK, fontweight="bold", va="top")
    ax.text(0.99, 0.92, f"mean {mean:.1f}° ± {sd:.1f}°", transform=ax.transAxes, color=INK_2, ha="right", va="top")
    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelcolor=INK_2)
    ax.margins(x=0)
    ax.set_ylim(0, 72)  # headroom so the panel labels sit above the data
    ax.set_yticks([0, 20, 40, 60])

fig.supylabel("Tilt from membrane normal (°)", color=INK_2, fontsize=10, x=0.015)
axes[-1].set_xlabel("Time (ns)", color=INK_2)
fig.suptitle("KALP helix tilt by build route and lipid", x=0.02, y=0.985, ha="left", color=INK, fontsize=12, fontweight="bold")
fig.text(0.02, 0.952, "Backbone residue 3 → 21 vector vs. z, 20 ns production each (self-assembly: 1 frame per ns)", color=INK_2, fontsize=9)
fig.subplots_adjust(top=0.925, bottom=0.06, left=0.1, right=0.97, hspace=0.15)

out = os.path.join("..", "figures", "tilt.png")
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, facecolor=SURFACE)
print(f"wrote {out}")
