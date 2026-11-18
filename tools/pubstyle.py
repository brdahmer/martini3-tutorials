"""Shared publication style for the figures in figures-publication/ folders.

Follows the scientific-figure-making house style: Arial/Helvetica, top and
right spines off, no grid, frameless legends, a blue/red/teal/violet palette,
tight_layout(pad=2), and export as 300 dpi PNG plus vector PDF.
"""
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

PALETTE = {
    "blue_main": "#0F4D92",
    "blue_secondary": "#3775BA",
    "red_strong": "#B64342",
    "red_1": "#F6CFCB",
    "teal": "#42949E",
    "violet": "#9A4D8E",
    "neutral": "#CFCECE",
    "dark": "#272727",
    "gray": "#767676",
}
# Fixed series order: key series first, then contrast. Dashes keep series
# distinct in grayscale print.
SERIES = [PALETTE["blue_main"], PALETTE["red_strong"]]
DASHES = ["-", "--"]


def apply_style(font_size=16, axes_linewidth=2.0):
    plt.rcParams.update({
        "font.family": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": font_size,
        "axes.linewidth": axes_linewidth,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "xtick.major.width": axes_linewidth,
        "ytick.major.width": axes_linewidth,
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "legend.frameon": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,  # embed TrueType so text stays editable
        "mathtext.default": "regular",
    })


def finalize(fig, out, pad=2.0):
    """Save `out` (path with or without extension) as both .png and .pdf."""
    base = Path(out).with_suffix("")
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(pad=pad)
    for ext in ("png", "pdf"):
        fig.savefig(base.with_suffix("." + ext), dpi=300, facecolor="white")
    plt.close(fig)
    print(f"wrote {base}.png/.pdf")
