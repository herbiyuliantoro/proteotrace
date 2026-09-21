#!/usr/bin/env python3
"""Shared journal figure style for ProteoTrace.

One place defines type, colour, sizing and export so the four main figures
read as a single system. Import and call `apply()` before building a figure.

Type          Arial (registered from the Windows font directory under WSL),
              8 pt base. Nothing in a main figure is set below 7.5 pt: any
              explanation that would need smaller type belongs in the caption.
Widths        `WIDTH_1COL` 85 mm, `WIDTH_15COL` 140 mm, `WIDTH_2COL` 180 mm.
Colour        Categorical slots 1-3 of the validated reference palette
              (blue/orange/aqua). That three-hue subset clears the all-pairs
              colour-vision gates (worst CVD dE 9.2, normal-vision dE 24.0),
              so it is safe for scatter and small-multiple use. Greys are
              non-series: they carry "baseline", "unrelated" or "raw", never
              a fourth category competing with the hues.
Export        600 dpi PNG plus vector PDF with TrueType (type 42) fonts, which
              is what publishers require for editable text.
"""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/proteotrace-mpl")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.transforms import ScaledTranslation

# ---------------------------------------------------------------- dimensions
MM = 1.0 / 25.4
WIDTH_1COL = 85 * MM
WIDTH_15COL = 140 * MM
WIDTH_2COL = 180 * MM

# -------------------------------------------------------------------- colour
BLUE = "#2a78d6"     # categorical slot 1
ORANGE = "#eb6834"   # categorical slot 2
AQUA = "#1baf7a"     # categorical slot 3
INK = "#1a1a19"      # primary text and mark outlines
MUTED = "#6b6a66"    # secondary text
GREY = "#9a9894"     # non-series grey: baseline / raw / unrelated
PALE = "#d8d7d2"     # de-emphasised marks
GRID = "#e6e5e1"     # recessive grid
SURFACE = "#ffffff"

SEQ_BLUE = "#eef4fc", "#2a78d6"  # light -> dark, single hue


def _register_arial() -> str:
    """Register Arial from the Windows font directory; fall back if absent."""
    candidates = [
        "/mnt/c/Windows/Fonts/arial.ttf",
        "/mnt/c/Windows/Fonts/arialbd.ttf",
        "/mnt/c/Windows/Fonts/ariali.ttf",
        "/mnt/c/Windows/Fonts/arialbi.ttf",
    ]
    found = False
    for path in candidates:
        if Path(path).exists():
            fm.fontManager.addfont(path)
            found = True
    return "Arial" if found else "DejaVu Sans"


def apply(base: float = 8.0) -> str:
    """Install the shared rcParams. Returns the resolved font family."""
    family = _register_arial()
    plt.rcParams.update({
        "font.family": family,
        "font.size": base,
        "axes.titlesize": base + 1,
        "axes.titleweight": "bold",
        "axes.titlepad": 4.0,
        "figure.titlesize": base + 1,
        "axes.titlelocation": "left",
        "axes.labelsize": base,
        "axes.labelcolor": INK,
        "axes.edgecolor": "#b9b8b3",
        "axes.linewidth": 0.6,
        "axes.labelpad": 2.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.labelsize": base - 0.5,
        "ytick.labelsize": base - 0.5,
        "xtick.color": INK,
        "ytick.color": INK,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.4,
        "ytick.major.size": 2.4,
        "xtick.major.pad": 2.0,
        "ytick.major.pad": 2.0,
        "legend.fontsize": base - 0.5,
        "legend.frameon": False,
        "legend.handlelength": 1.4,
        "legend.handletextpad": 0.5,
        "legend.borderaxespad": 0.2,
        "legend.labelspacing": 0.35,
        "legend.columnspacing": 1.1,
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "lines.linewidth": 1.2,
        "lines.markersize": 4.0,
        "patch.linewidth": 0.6,
        "figure.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "savefig.dpi": 600,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })
    return family


def _offset(ax, dx_pt: float, dy_pt: float):
    """Transform anchored at the axes top-left, displaced in typographic points.

    Offsets in points keep the header block identically spaced whatever the
    panel height is; fractions of the axes do not.
    """
    return ax.transAxes + ScaledTranslation(
        dx_pt / 72.0, dy_pt / 72.0, ax.figure.dpi_scale_trans)


def head(ax, letter: str | None, title: str | None = None,
         sub: str | None = None, letter_dx: float = -30.0) -> None:
    """Panel letter, title and one recessive subtitle, stacked without collision.

    Use this instead of `ax.set_title` so the three elements share one ruler.
    """
    if letter:
        ax.text(0, 1, letter, transform=_offset(ax, letter_dx, 9.0), fontsize=11,
                fontweight="bold", color=INK, va="baseline", ha="left")
    if title:
        ax.text(0, 1, title, transform=_offset(ax, 0.0, 9.0), fontsize=9,
                fontweight="bold", color=INK, va="baseline", ha="left")
    if sub:
        # Retained for compatibility. The main figures carry no subtitle line:
        # design facts and interval definitions belong in the caption.
        ax.text(0, 1, sub, transform=_offset(ax, 0.0, 6.0), fontsize=7,
                color=MUTED, va="baseline", ha="left")


def panel(ax, letter: str, dx: float = -30.0) -> None:
    head(ax, letter, letter_dx=dx)


def subtitle(ax, text: str, dy: float = 6.0) -> None:
    ax.text(0, 1, text, transform=_offset(ax, 0.0, dy), fontsize=6,
            color=MUTED, va="baseline", ha="left")


def ygrid(ax, axis: str = "y") -> None:
    ax.grid(axis=axis, color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)


def save(fig, out_dir: Path, stem: str) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / f"{stem}.png"
    pdf = out_dir / f"{stem}.pdf"
    fig.savefig(png, dpi=600)
    fig.savefig(pdf)
    plt.close(fig)
    return png, pdf
