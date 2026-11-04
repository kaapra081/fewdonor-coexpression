"""Shared figure helpers for paper2/make_fig*.py (stand-alone: no external plotting skill needed).

apply_figure_style(sizes=(base, secondary, tick)) : rcParams for the paper figures (DejaVu Sans, Type-42 fonts, open spines).
                                                   Unlike the 'tight' default of the earlier export, figures are saved at their
                                                   *design* size (6.5 in wide = \\linewidth), so printed size == nominal size.
panel_letter(ax, letter, ...)                     : bold lower-case panel letter outside the top-left of an axes.
audit_pdf(path, width_in=6.5, min_pt=7.0)          : reads the exported PDF with PyMuPDF and reports page size, the effective
                                                   (printed) size of every text span, and any span below min_pt.
overlap_report(fig)                                : pairs of visible text boxes that overlap, and texts outside the canvas.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt

PAGE_W = 6.5   # in; \linewidth of the paper (article, 10 pt, 1-in margins, US Letter)


def apply_figure_style(sizes=(8, 7, 7), font=None):
    base, secondary, tick = sizes
    rc = {
        "font.family": "sans-serif", "font.size": base, "axes.labelsize": base, "axes.titlesize": base,
        "legend.fontsize": secondary, "xtick.labelsize": tick, "ytick.labelsize": tick,
        "axes.linewidth": 0.6, "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False, "legend.frameon": False,
        "figure.dpi": 200, "savefig.dpi": 300, "savefig.bbox": "standard",   # NOT 'tight': keep the design width
        "axes.titleweight": "normal", "axes.titlelocation": "left", "axes.labelweight": "normal",
        "lines.linewidth": 1.2, "patch.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42,
        "mathtext.fontset": "dejavusans",
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
    mpl.rcParams.update(rc)


def panel_letter(ax, letter, dx=-0.18, dy=1.02, fontsize=None):
    fontsize = fontsize or mpl.rcParams["font.size"] + 2
    ax.text(dx, dy, letter, transform=ax.transAxes, fontweight="bold", fontsize=fontsize, va="bottom", ha="left")


def overlap_report(fig):
    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
    ov = [(a.get_text()[:30], b.get_text()[:30]) for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:] if ba.overlaps(bb)]
    fb = fig.bbox
    out = [t.get_text()[:30] for t, bb in texts if bb.x0 < fb.x0 - 0.5 or bb.y0 < fb.y0 - 0.5 or bb.x1 > fb.x1 + 0.5 or bb.y1 > fb.y1 + 0.5]
    return ov, out


def audit_pdf(path, width_in=PAGE_W, min_pt=7.0, verbose=True):
    import pymupdf
    pg = pymupdf.open(path)[0]
    W, H = pg.rect.width / 72, pg.rect.height / 72
    sc = width_in / W
    spans = [(s["size"] * sc, s["text"].strip()) for b in pg.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"] if s["text"].strip()]
    sizes = sorted({round(a, 2) for a, _ in spans})
    low = [(round(a, 2), t) for a, t in spans if a < min_pt - 0.05]
    res = dict(page_in=(round(W, 3), round(H, 3)), scale=round(sc, 4), printed_height_in=round(H * sc, 3), min_eff_pt=round(min(a for a, _ in spans), 2),
               sizes_eff_pt=sizes, n_spans=len(spans), below_min=low)
    if verbose:
        print(path, res)
    return res


def enlarge_mathtext_subscripts(factor=0.85):
    """mathtext shrinks sub/superscripts to 0.7 x the parent size (hard-wired module constant). With factor 0.85, 8.5-pt math has 7.2-pt
    subscripts, so the figure meets a 7-pt floor without 10-pt parent text. Silently does nothing if the constant is absent (then audit_pdf flags it)."""
    try:
        import matplotlib._mathtext as mt
        mt.SHRINK_FACTOR = factor
    except Exception:
        pass
