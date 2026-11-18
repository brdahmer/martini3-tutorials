#!/usr/bin/env python3
"""Plot a GROMACS .xvg file in publication style (PNG + PDF).

Same idea as plot_xvg.py, but with no title on the figure (the caption carries
it), thicker axes, no grid, and vector export. Axis labels and legends come
from the xvg header unless given on the command line.

Example:
    python3 plot_xvg_pub.py area.xvg -o ../figures-publication/area-per-lipid \
        --ylabel "Area per lipid (nm$^2$)" --ns
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib.pyplot as plt  # noqa: E402

from plot_xvg import read_xvg  # noqa: E402
from pubstyle import DASHES, PALETTE, SERIES, apply_style, finalize  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("xvg")
    p.add_argument("-o", "--output", required=True, help="output path; .png and .pdf are both written")
    p.add_argument("--xlabel")
    p.add_argument("--ylabel")
    p.add_argument("--legend", nargs="+", help="one name per data column after x")
    p.add_argument("--ns", action="store_true", help="convert a time axis in ps to ns")
    p.add_argument("--columns", type=int, nargs="+", help="data columns to plot (1 = first after x)")
    p.add_argument("--note", help="short in-plot annotation, e.g. a fitted value")
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

    apply_style()
    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    for i in range(ys.shape[1]):
        ax.plot(x, ys[:, i], color=SERIES[i], linestyle=DASHES[i], linewidth=2.5, label=legends[i])

    ax.set_xlabel(xlabel)
    ax.set_ylabel(a.ylabel or meta.get("ylabel", ""))
    ax.margins(x=0)
    # Legend and note sit above the axes so they never cover data.
    y_top = 1.02
    if ys.shape[1] > 1:
        ax.legend(loc="lower left", bbox_to_anchor=(0, y_top), ncol=ys.shape[1], handlelength=2.2,
                  borderaxespad=0, columnspacing=1.5, fontsize=14)
        y_top = 1.17
    if a.note:
        ax.text(0, y_top, a.note, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=13, color=PALETTE["dark"])
    finalize(fig, a.output)


if __name__ == "__main__":
    main()
