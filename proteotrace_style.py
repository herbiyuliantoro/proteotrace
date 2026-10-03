#!/usr/bin/env python3
"""Shared journal figure style for ProteoTrace.

One place defines type, colour, sizing and export so the main figures read
as a single system. Import and call `apply()` before building a figure.

Sizes follow the Journal of Proteome Research (ACS) author guidelines:
single-column graphics up to 240 pt (3.33 in) wide, double-column between
300 and 504 pt (7 in), at most 660 pt deep including the caption, lettering
no smaller than 4.5 pt in Helvetica or Arial, lines no thinner than 0.5 pt,
colour art at least 300 dpi. `save()` checks these before writing anything.

Type          Arial (registered from the Windows font directory under WSL),
              8 pt base. Nothing in a main figure is set below 7.5 pt: any
              explanation that would need smaller type belongs in the caption.
Widths        `WIDTH_1COL` 84 mm (238 pt), `WIDTH_2COL` 176 mm (499 pt).
Axes          Black left and bottom spines, outward ticks, no grid.
Colour        Categorical slots 1-3 of the validated reference palette
              (blue/orange/aqua). That three-hue subset clears the all-pairs
              colour-vision gates (worst CVD dE 9.2, normal-vision dE 24.0),
              so it is safe for scatter and small-multiple use. Greys are
              non-series: they carry "baseline", "unrelated" or "raw", never
              a fourth category competing with the hues.
Export        Vector PDF with TrueType (type 42) fonts for editable text,
              600 dpi TIFF (LZW) for submission, and 600 dpi PNG for viewing.
"""
from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/proteotrace-mpl")

import matplotlib

matplotlib.use("Agg")
import matplotlib.collections
import matplotlib.lines
import matplotlib.text
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager as fm
from matplotlib.transforms import ScaledTranslation

# ---------------------------------------------------------------- dimensions
MM = 1.0 / 25.4
WIDTH_1COL = 84 * MM          # JPR single column: up to 240 pt
WIDTH_2COL = 176 * MM         # JPR double column: 300-504 pt
MAX_DEPTH_PT = 660            # JPR maximum depth, caption included
MIN_FONT_PT = 4.5             # JPR minimum lettering
MIN_LINE_PT = 0.5             # JPR minimum line weight
MIN_MARGIN_MM = 1.5           # house rule: nothing closer to the edge than this

# -------------------------------------------------------------------- colour
BLUE = "#2a78d6"     # categorical slot 1
ORANGE = "#eb6834"   # categorical slot 2
AQUA = "#1baf7a"     # categorical slot 3
INK = "#1a1a19"      # primary text and mark outlines
MUTED = "#6b6a66"    # secondary text
GREY = "#9a9894"     # non-series grey: baseline / raw / unrelated
PALE = "#d8d7d2"     # de-emphasised marks
GRID = "#e6e5e1"     # dividers only; main figures carry no grid
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
        "axes.edgecolor": INK,
        "axes.linewidth": 0.6,
        "axes.grid": False,
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
        "xtick.direction": "out",
        "ytick.direction": "out",
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
    """Retained for the archived figure scripts; the main figures use no grid."""
    ax.grid(axis=axis, color=GRID, linewidth=0.5)
    ax.set_axisbelow(True)


def check_journal(fig) -> None:
    """Refuse a figure that breaks the JPR artwork rules.

    Width must fit one column (<= 240 pt) or two (300-504 pt), depth must
    leave room for a caption under 660 pt, every visible text must be at least
    4.5 pt and every visible stroked line at least 0.5 pt, and nothing may sit
    closer than 1.5 mm to the edge.
    """
    width_pt, height_pt = (v * 72 for v in fig.get_size_inches())
    if not (width_pt <= 240.5 or 300 <= width_pt <= 504):
        raise ValueError(f"figure width {width_pt:.1f} pt fits neither JPR column format")
    if height_pt > MAX_DEPTH_PT - 60:
        raise ValueError(f"figure depth {height_pt:.1f} pt leaves no room for a caption")
    small = [t.get_text() for t in fig.findobj(matplotlib.text.Text)
             if t.get_visible() and t.get_text().strip() and t.get_fontsize() < MIN_FONT_PT]
    if small:
        raise ValueError(f"text below {MIN_FONT_PT} pt: {small[:3]}")
    thin = [line for line in fig.findobj(matplotlib.lines.Line2D)
            if line.get_visible() and line.get_linestyle() not in ("None", "")
            and 0 < line.get_linewidth() < MIN_LINE_PT and len(line.get_xdata()) > 1]
    for coll in fig.findobj(matplotlib.collections.LineCollection):
        drawn = coll.get_visible() and len(coll.get_segments()) > 0   # colorbar dividers are empty
        if drawn and any(0 < w < MIN_LINE_PT for w in coll.get_linewidths()):
            thin.append(coll)
    if thin:
        raise ValueError(f"{len(thin)} line(s) thinner than {MIN_LINE_PT} pt")

    # Rendered margins: nothing may come closer to the edge than MIN_MARGIN_MM,
    # or the figure reads as cut off and may be clipped in production.
    dpi = 200
    fig.set_dpi(dpi)
    fig.canvas.draw()
    rgba = np.asarray(fig.canvas.buffer_rgba())
    ink = (rgba[..., :3].min(axis=2) < 200) & (rgba[..., 3] > 0)
    cols, rows = np.where(ink.any(axis=0))[0], np.where(ink.any(axis=1))[0]
    mm = 25.4 / dpi
    margins = {"left": cols.min() * mm, "right": (ink.shape[1] - 1 - cols.max()) * mm,
               "top": rows.min() * mm, "bottom": (ink.shape[0] - 1 - rows.max()) * mm}
    tight = {k: round(v, 1) for k, v in margins.items() if v < MIN_MARGIN_MM}
    if tight:
        raise ValueError(f"content within {MIN_MARGIN_MM} mm of the edge (mm): {tight}")


def save(fig, out_dir: Path, stem: str) -> tuple[Path, Path]:
    check_journal(fig)
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / f"{stem}.png"
    pdf = out_dir / f"{stem}.pdf"
    tif = out_dir / f"{stem}.tif"
    fig.savefig(pdf)
    fig.savefig(png, dpi=600)
    fig.savefig(tif, dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    return png, pdf
