#!/usr/bin/env python3
"""Plot the KALP dimer run (section 5.2): closest distance between the two
peptides, and the angle between their helix axes, over time.

Run from this folder:  python3 plot_dimer.py
Writes ../figures/dimer.png (same style as tilt.png) and
../figures-publication/dimer.png/.pdf (same style as plot_tilt_pub.py).
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "tools"))
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

CONTACT = 0.6  # nm; closer than this counts as touching
WINDOW = 250  # frames = 5 ns at 20 ps per frame


def load(path):
    d = np.loadtxt(path, comments=("#", "@"))
    return d[:, 0] / 1000.0, d[:, 1]


def running_mean(t, y, n):
    k = np.ones(n) / n
    return t[n // 2:n // 2 + len(y) - n + 1], np.convolve(y, k, mode="valid")


t, dist = load("dimer-mindist.xvg")
_, angle = load("dimer-angle.xvg")


def draw(axes, color, raw_alpha, lw, ink, note_kw):
    ax1, ax2 = axes
    ax1.axhspan(0, CONTACT, color=color, alpha=0.12, lw=0)
    ax1.plot(t, dist, color=color, alpha=raw_alpha, linewidth=0.6)
    ax1.plot(*running_mean(t, dist, WINDOW), color=color, linewidth=lw)
    ax1.text(0.99, CONTACT + 0.1, f"contact (< {CONTACT} nm)", transform=ax1.get_yaxis_transform(),
             ha="right", va="bottom", color=ink, **note_kw)
    ax1.text(0.99, 0.97, f"closest {dist.min():.1f} nm", transform=ax1.transAxes,
             ha="right", va="top", color=ink, **note_kw)
    ax1.set_ylim(0, 7)
    ax1.set_ylabel("Min. distance (nm)")

    ax2.plot(t, angle, color=color, alpha=raw_alpha, linewidth=0.6)
    ax2.plot(*running_mean(t, angle, WINDOW), color=color, linewidth=lw)
    ax2.text(0.99, 0.97, f"mean {angle.mean():.0f}° ± {angle.std():.0f}°", transform=ax2.transAxes,
             ha="right", va="top", color=ink, **note_kw)
    ax2.set_ylim(0, 90)
    ax2.set_yticks([0, 30, 60, 90])
    ax2.set_ylabel("Helix–helix angle (°)")
    ax2.set_xlabel("Time (ns)")
    for ax in axes:
        ax.set_xlim(0, t[-1])


# Original style, matching plot_tilt.py ----------------------------------------
SURFACE, INK, INK_2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
with plt.rc_context({"font.size": 10, "font.family": "sans-serif"}):
    fig, axes = plt.subplots(2, 1, figsize=(7.0, 5.0), dpi=200, sharex=True)
    fig.patch.set_facecolor(SURFACE)
    draw(axes, "#2a78d6", 0.35, 1.5, INK_2, {})
    for ax in axes:
        ax.set_facecolor(SURFACE)
        ax.grid(True, color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(AXIS)
        ax.tick_params(colors=MUTED, labelcolor=INK_2)
        ax.yaxis.label.set_color(INK_2)
        ax.xaxis.label.set_color(INK_2)
    fig.suptitle("Two KALP peptides in POPC do not dimerize in 250 ns", x=0.02, y=0.98, ha="left",
                 color=INK, fontsize=12, fontweight="bold")
    fig.text(0.02, 0.925, "Raw values (light) and 5 ns running mean (bold), one replica", color=INK_2, fontsize=9)
    fig.subplots_adjust(top=0.88, bottom=0.1, left=0.1, right=0.97, hspace=0.15)
    os.makedirs(os.path.join("..", "figures"), exist_ok=True)
    out = os.path.join("..", "figures", "dimer.png")
    fig.savefig(out, facecolor=SURFACE)
    plt.close(fig)
    print(f"wrote {out}")

# Publication style -------------------------------------------------------------
from pubstyle import PALETTE, apply_style, finalize  # noqa: E402

apply_style(font_size=15, axes_linewidth=2.0)
fig, axes = plt.subplots(2, 1, figsize=(7.5, 6.0), sharex=True)
draw(axes, PALETTE["blue_main"], 0.3, 2.5, PALETTE["dark"], {"fontsize": 13})
for ax, letter in zip(axes, "ab"):
    ax.text(0.01, 0.97, f"({letter})", transform=ax.transAxes, va="top", fontweight="bold")
finalize(fig, os.path.join("..", "figures-publication", "dimer"), pad=1.0)
