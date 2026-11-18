#!/usr/bin/env python3
"""Plot a GROMACS .xvg file as a PNG.

Titles, axis labels and legends are read from the xvg header when it has one.
Any of them can be overridden on the command line, which is needed for files
written with `-xvg no` or by numpy (they have no header).

Example:
    python3 plot_xvg.py area.xvg -o ../figures/area-per-lipid.png \
        --title "Area per lipid" --ylabel "Area per lipid (nm$^2$)" --ns
"""
import argparse
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Colors from the reference dataviz palette (light mode).
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SERIES = ["#2a78d6", "#eb6834"]  # categorical slots 1 and 2, in fixed order
DASHES = ["-", "--"]  # second encoding so series are not told apart by color alone


def grace_to_mathtext(text):
    """Convert xmgrace markup like nm\\S2\\N or m\\S-3\\N to matplotlib mathtext."""
    text = text.replace("\\\\", "\\")
    text = re.sub(r"\\S(.*?)\\N", r"$^{\1}$", text)
    text = re.sub(r"\\s(.*?)\\N", r"$_{\1}$", text)
    return text


def read_xvg(path):
    meta = {"legends": {}}
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if line.startswith("@"):
                m = re.match(r'@\s+(title|subtitle)\s+"(.*)"', line)
                if m:
                    meta[m.group(1)] = grace_to_mathtext(m.group(2))
                m = re.match(r'@\s+([xy])axis\s+label\s+"(.*)"', line)
                if m:
                    meta[m.group(1) + "label"] = grace_to_mathtext(m.group(2))
                m = re.match(r'@\s+s(\d+)\s+legend\s+"(.*)"', line)
                if m:
                    meta["legends"][int(m.group(1))] = grace_to_mathtext(m.group(2))
                continue
            if line.strip():
                rows.append([float(v) for v in line.split()])
    return np.array(rows), meta


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("xvg")
    p.add_argument("-o", "--output", required=True, help="output .png path")
    p.add_argument("--title")
    p.add_argument("--subtitle")
    p.add_argument("--xlabel")
    p.add_argument("--ylabel")
    p.add_argument("--legend", nargs="+", help="one name per data column after x")
    p.add_argument("--ns", action="store_true", help="convert a time axis in ps to ns")
    p.add_argument("--columns", type=int, nargs="+", help="data columns to plot (1 = first after x); default all")
    a = p.parse_args()

    data, meta = read_xvg(a.xvg)
    x, ys = data[:, 0], data[:, 1:]
    if a.columns:
        ys = ys[:, [c - 1 for c in a.columns]]
        meta["legends"] = {i: meta["legends"].get(c - 1, "") for i, c in enumerate(a.columns)}
    xlabel = a.xlabel or meta.get("xlabel", "")
    if a.ns:
        x = x / 1000.0
        xlabel = a.xlabel or re.sub(r"\(ps\)", "(ns)", xlabel)
    legends = a.legend or [meta["legends"].get(i, f"Series {i + 1}") for i in range(ys.shape[1])]
    if ys.shape[1] > len(SERIES):
        raise SystemExit(f"{a.xvg}: {ys.shape[1]} series; this script styles at most {len(SERIES)}")

    plt.rcParams.update({"font.size": 10, "font.family": "sans-serif"})
    fig, ax = plt.subplots(figsize=(6.4, 4.0), dpi=200)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    for i in range(ys.shape[1]):
        ax.plot(x, ys[:, i], color=SERIES[i], linestyle=DASHES[i], linewidth=1.5, label=legends[i])

    title = a.title or meta.get("title", "")
    subtitle = a.subtitle if a.subtitle is not None else meta.get("subtitle", "")
    ax.set_title(title, loc="left", color=INK, fontsize=12, fontweight="bold", pad=18 if subtitle else 8)
    if subtitle:
        ax.text(0, 1.02, subtitle, transform=ax.transAxes, color=INK_2, fontsize=9, va="bottom")
    ax.set_xlabel(xlabel, color=INK_2)
    ax.set_ylabel(a.ylabel or meta.get("ylabel", ""), color=INK_2)

    ax.grid(True, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelcolor=INK_2)
    ax.margins(x=0)

    if ys.shape[1] > 1:  # a single series is named by the title, so it gets no legend
        ax.legend(frameon=False, labelcolor=INK_2)

    fig.tight_layout()
    fig.savefig(a.output, facecolor=SURFACE)
    print(f"wrote {a.output}")


if __name__ == "__main__":
    main()
