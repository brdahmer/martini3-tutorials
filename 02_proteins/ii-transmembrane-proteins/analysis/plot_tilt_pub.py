#!/usr/bin/env python3
"""Publication-style version of plot_tilt.py.

Each panel shows the raw tilt (light), a 1 ns running mean (bold), and the
3-20 ns mean ± SD as a dashed line and band. The first 3 ns are shaded as
equilibration, matching the averages in the README.

Run from this folder:  python3 plot_tilt_pub.py
Writes ../figures-publication/tilt.png and tilt.pdf
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "tools"))
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from pubstyle import PALETTE, apply_style, finalize  # noqa: E402

ROUTES = [
    ("POPC, Insane", "insane-angles.xvg", PALETTE["blue_main"]),
    ("POPC, COBY", "coby-angles.xvg", PALETTE["teal"]),
    ("POPC, TS2CG", "ts2cg-angles.xvg", PALETTE["red_strong"]),
    ("DLPC, Insane", "dlpc-angles.xvg", PALETTE["violet"]),
]
T_EQ = 3.0  # ns excluded from the average
WINDOW = 50  # frames = 1 ns at 20 ps per frame


def load(path):
    d = np.loadtxt(path, comments=("#", "@"))
    return d[:, 0] / 1000.0, 180.0 - d[:, 1]  # ps -> ns, angle -> tilt


def running_mean(y, n):
    kernel = np.ones(n) / n
    return np.convolve(y, kernel, mode="valid")


apply_style(font_size=15, axes_linewidth=2.0)
fig, axes = plt.subplots(2, 2, figsize=(10, 6.4), sharex=True, sharey=True)

for ax, letter, (name, path, color) in zip(axes.flat, "abcd", ROUTES):
    t, tilt = load(path)
    keep = t >= T_EQ
    mean, sd = tilt[keep].mean(), tilt[keep].std()

    ax.axvspan(0, T_EQ, color=PALETTE["neutral"], alpha=0.35, lw=0)
    ax.fill_between([T_EQ, t[-1]], mean - sd, mean + sd, color=color, alpha=0.12, lw=0)
    ax.plot(t, tilt, color=color, alpha=0.3, linewidth=0.8)
    ax.plot(t[WINDOW // 2:WINDOW // 2 + len(tilt) - WINDOW + 1], running_mean(tilt, WINDOW), color=color, linewidth=2.5)
    ax.hlines(mean, T_EQ, t[-1], color=PALETTE["dark"], linewidth=1.5, linestyle="--")

    ax.text(0.03, 0.96, f"({letter}) {name}", transform=ax.transAxes, va="top", fontweight="bold")
    ax.text(0.97, 0.96, f"{mean:.1f}° ± {sd:.1f}°", transform=ax.transAxes, va="top", ha="right",
            color=PALETTE["dark"])
    ax.set_xlim(0, t[-1])
    ax.set_ylim(0, 60)
    ax.set_yticks([0, 20, 40, 60])

for ax in axes[-1]:
    ax.set_xlabel("Time (ns)")
for ax in axes[:, 0]:
    ax.set_ylabel("Tilt (°)")

finalize(fig, os.path.join("..", "figures-publication", "tilt"), pad=1.0)
