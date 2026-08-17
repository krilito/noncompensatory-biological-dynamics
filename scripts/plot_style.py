"""Central figure style configuration — PUBLIC_FIGURE_STYLE_V1.

One module shared by every figure so cohort colour, response encoding and
analysis-unit grammar cannot drift between panels.

Design constraints honoured here:
  * colourblind-safe; colour is NEVER the only encoding for response
  * patient-level results are visually distinct from sample-level results
  * grayscale-legible (lightness of the two response colours differs by ~25 L*)
  * no figure-level "Figure X" title (panel letters only)
  * PNG previews <= 200 dpi so the whole figure is visible on one screen;
    PDF/SVG stay vector for print
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

STYLE_ID = "PUBLIC_FIGURE_STYLE_V1"

# ------------------------------------------------------------------ canvas geometry
MM = 1.0 / 25.4
W_SINGLE = 89 * MM          # 3.504 in — Nature single column
W_ONEHALF = 136 * MM        # 5.354 in
W_DOUBLE = 183 * MM         # 7.205 in — Nature double column

# ------------------------------------------------------------------ typography
FONT_FAMILY = "sans-serif"
FONT_STACK = ["DejaVu Sans", "Arial", "Helvetica", "Liberation Sans"]
FS_PANEL_LABEL = 9          # a, b, c ...
FS_TITLE = 7.5              # panel sub-title
FS_AXIS_LABEL = 7.0
FS_TICK = 6.5
FS_ANNOT = 6.5
FS_LEGEND = 6.5
FS_MIN = 6.5                # nothing smaller than this is permitted

# ------------------------------------------------------------------ line / marker
LW_SPINE = 0.8
LW_MAIN = 1.0
LW_CI = 0.9
LW_REF = 0.8                # zero / 0.5 reference lines
LW_LINK = 0.6               # paired connector lines
MS_POINT = 4.5
MS_SMALL = 3.2

# ------------------------------------------------------------------ semantic colour
# Response color convention: responder cool, non-responder warm.
C_RESPONDER = "#1E6E8C"      # cool teal-blue     L* ~ 42
C_NONRESPONDER = "#D17A22"   # warm ochre         L* ~ 58
C_RESPONSE = {"R": C_RESPONDER, "NR": C_NONRESPONDER,
              "Responder": C_RESPONDER, "Nonresponder": C_NONRESPONDER,
              1: C_RESPONDER, 0: C_NONRESPONDER}

# Redundant, colour-independent encoding for response
MK_RESPONSE = {"R": "o", "NR": "D", "Responder": "o", "Nonresponder": "D", 1: "o", 0: "D"}
LS_RESPONSE = {"R": "-", "NR": (0, (3, 1.4)), "Responder": "-", "Nonresponder": (0, (3, 1.4))}

# Ink / structure
C_INK = "#1A1A1A"
C_MID = "#5A5A5A"
C_GREY = "#8C8C8C"
C_LIGHT = "#D8D8D8"
C_PANEL_BG = "#F4F4F2"
C_REF = "#4D4D4D"            # zero / chance reference line

# Direction accents for delta_L and other signed quantities.
# Deliberately NOT red/green.
C_HELPED = "#2F6F5E"         # adding change lowered loss  (delta_L > 0)
C_HURT = "#8C3B2E"           # adding change raised loss   (delta_L < 0)

# Cohort colours — FIXED across every figure
C_COHORT = {
    "GSE91061": "#4C4C4C",
    "PRJEB23709": "#1E6E8C",
    "MORRISON-1-public": "#2F6F5E",
    "MGH": "#7A5C99",
    "Sade-Feldman": "#B5651D",
    "GSE78220": "#8C8C8C",
    "GSE282471": "#6E6E6E",
    "IMvigor210": "#8C3B2E",
}

# MELD axis colours
C_AXIS = {"T": "#8C3B2E", "E": "#1E6E8C", "X": "#8C8C8C", "C": "#2F6F5E"}

# ------------------------------------------------------------------ analysis-unit grammar
# Sample-level  : open marker, thin stroke
# Patient-level : FILLED SQUARE, heavier stroke  -> never confusable with a sample row
# Record-level  : open marker with a vertical tick
# Not estimable : grey open circle with an x, no interval
UNIT_STYLE = {
    "sample": dict(marker="o", markersize=MS_POINT, markerfacecolor="white",
                   markeredgewidth=1.0),
    "record": dict(marker="v", markersize=MS_POINT, markerfacecolor="white",
                   markeredgewidth=1.0),
    "patient": dict(marker="s", markersize=MS_POINT + 0.6, markeredgewidth=1.3),
    "pooled_cell": dict(marker="h", markersize=MS_POINT, markerfacecolor="white",
                        markeredgewidth=1.0),
    "not_estimable": dict(marker="x", markersize=MS_POINT, markeredgewidth=1.2,
                          color=C_GREY),
}

# ------------------------------------------------------------------ status grammar
# Text + shape/line, never colour alone.
STATUS_GLYPH = {
    "SUPPORTED": "●",
    "NOT_SUPPORTED": "○",
    "CONFLICT": "▲",
    "NOT_TESTED": "–",
    "DESCRIPTIVE": "◆",
    "DECLARED": "□",
}

LEGEND_ORDER = ["Responder", "Non-responder", "sample-level", "record-level",
                "patient-level", "not estimable"]


def apply() -> None:
    """Install the global rcParams. Call once at the top of every figure script."""
    mpl.rcParams.update({
        "figure.dpi": 100,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "svg.hashsalt": "meld-icb-submission-sot-v1",
        "figure.facecolor": "white",
        "axes.facecolor": "white",

        "font.family": FONT_FAMILY,
        "font.sans-serif": FONT_STACK,
        "font.size": FS_TICK,
        "axes.labelsize": FS_AXIS_LABEL,
        "axes.titlesize": FS_TITLE,
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "legend.fontsize": FS_LEGEND,

        "axes.labelcolor": C_INK,
        "axes.edgecolor": C_MID,
        "axes.linewidth": LW_SPINE,
        "axes.titlelocation": "left",
        "axes.titlepad": 4.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,

        "xtick.color": C_MID,
        "ytick.color": C_MID,
        "xtick.major.width": LW_SPINE,
        "ytick.major.width": LW_SPINE,
        "xtick.major.size": 2.4,
        "ytick.major.size": 2.4,
        "xtick.direction": "out",
        "ytick.direction": "out",

        "lines.linewidth": LW_MAIN,
        "lines.markersize": MS_POINT,
        "legend.frameon": False,
        "legend.handlelength": 1.4,
        "legend.handletextpad": 0.5,
        "legend.columnspacing": 1.1,
        "legend.borderaxespad": 0.2,

        "pdf.fonttype": 42,      # embed TrueType so text stays editable/searchable
        "ps.fonttype": 42,
        "svg.fonttype": "none",  # keep SVG text as text, not paths
    })


def panel_label(ax, letter: str, dx: float = -0.085, dy: float = 1.06) -> None:
    """Bold lowercase panel letter in the top-left, in axes coordinates."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=FS_PANEL_LABEL,
            fontweight="bold", va="bottom", ha="left", color=C_INK)


def panel_title(ax, text: str) -> None:
    ax.set_title(text, fontsize=FS_TITLE, color=C_INK, loc="left", pad=4)


def subtle(ax) -> None:
    """Trim an axes to the restrained house style."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(C_MID)
        ax.spines[side].set_linewidth(LW_SPINE)


def refline(ax, value: float, axis: str = "x", label: str | None = None) -> None:
    """Neutral, strong reference line (zero for deltas, 0.5 for AUC)."""
    if axis == "x":
        ax.axvline(value, color=C_REF, lw=LW_REF, ls=(0, (2.5, 2.0)), zorder=0)
    else:
        ax.axhline(value, color=C_REF, lw=LW_REF, ls=(0, (2.5, 2.0)), zorder=0)
    if label:
        ax.annotate(label, xy=(value, 1.0), xycoords=("data", "axes fraction"),
                    xytext=(2, -8), textcoords="offset points",
                    fontsize=FS_ANNOT, color=C_REF, ha="left", va="top")


def export(fig, outdir: Path, stem: str, png_dpi: int = 200) -> dict:
    """Write PDF + screen-friendly PNG. Returns {ext: path}."""
    outdir.mkdir(parents=True, exist_ok=True)
    paths = {}
    p = outdir / f"{stem}.pdf"
    fig.savefig(p, format="pdf", metadata={"CreationDate": None, "ModDate": None})
    paths["pdf"] = p
    p = outdir / f"{stem}.png"
    fig.savefig(p, format="png", dpi=png_dpi)
    paths["png"] = p
    plt.close(fig)
    return paths


# ------------------------------------------------------------------ wireframe mode
def apply_wireframe() -> None:
    """Grayscale, low-style structural mode for the pre-styling QC pass."""
    apply()
    mpl.rcParams.update({
        "axes.prop_cycle": mpl.cycler(color=[C_MID]),
        "savefig.dpi": 150,
    })


WIREFRAME_PALETTE = dict(
    responder=C_MID, nonresponder=C_GREY, helped=C_MID, hurt="#3A3A3A",
)
