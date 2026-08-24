"""Shared plotting helpers for the paper figures."""
import numpy as np
from matplotlib.ticker import FuncFormatter
from matplotlib.patches import Patch

# plain (non-scientific) tick labels for log axes, e.g. 0.01, 0.1, 1, 10, 100
plain_log_fmt = FuncFormatter(lambda y, _: f"{y:g}")

def gradient_vspan(ax, core, fade, color="orange", alpha_max=0.3, n=200, zorder=0):
    """Shade a solid band over `core`=(lo, hi), fading linearly to transparent over `fade`
    x-units on each side. Returns a proxy Patch for the legend."""
    lo, hi = core
    flo, fhi = lo - fade, hi + fade
    edges = np.linspace(flo, fhi, n + 1)
    for x0, x1 in zip(edges[:-1], edges[1:]):
        xm = 0.5 * (x0 + x1)
        if xm < lo:
            a = alpha_max * (xm - flo) / (lo - flo)
        elif xm > hi:
            a = alpha_max * (fhi - xm) / (fhi - hi)
        else:
            a = alpha_max
        ax.axvspan(x0, x1, color=color, alpha=max(a, 0.0), lw=0, zorder=zorder)
    return Patch(facecolor=color, alpha=alpha_max)
