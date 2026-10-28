"""Figures 1-5 of the paper (Section 2, illustration of the initial transformation).

Vector PDFs typeset with LaTeX (Computer Modern), 4.7 in wide, for
\\includegraphics[width=0.7\\textwidth]{...} or similar.

Requirements: Python 3, numpy, scipy, matplotlib (pgf backend) and a TeX
installation with pdflatex and amsmath.

    python3 make_figures.py            # writes Figure_1.pdf ... Figure_5.pdf
"""
import numpy as np
import matplotlib

matplotlib.use("pgf")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MultipleLocator
from scipy.optimize import brentq
from scipy.special import erfc

BLUE = "#2a78d6"   # curves
INK = "#0b0b0b"    # axis labels, zero markers
MUTED = "#52514e"  # tick labels
AXIS = "#8a8984"   # spines, ticks, zero line
GRID = "#e4e3df"   # grid

plt.rcParams.update({
    "pgf.texsystem": "pdflatex",
    "pgf.rcfonts": False,
    "pgf.preamble": r"\usepackage{amsmath}",
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 9,
    "axes.labelsize": 10,
    "axes.labelcolor": INK,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "xtick.labelcolor": MUTED,
    "ytick.labelcolor": MUTED,
    "xtick.color": AXIS,
    "ytick.color": AXIS,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.color": GRID,
    "grid.linewidth": 0.5,
    "lines.linewidth": 1.0,
})

WIDTH = 4.7  # inches


def f(x):
    """Test function (2.6) of the paper."""
    return np.sin(np.pi * x / 2) + np.sin(4 * np.pi * x) ** 2 + 0.1


def real_zeros(func, a, b, n=200_001):
    """All sign changes of func on [a, b], refined by Brent's method."""
    x = np.linspace(a, b, n)
    y = func(x)
    idx = np.nonzero(np.sign(y[:-1]) * np.sign(y[1:]) < 0)[0]
    return np.array([brentq(func, x[i], x[i + 1], xtol=1e-14) for i in idx])


label = FuncFormatter(lambda v, _: f"${v:g}$")


def new_axes(height, xlim, ylim, xstep, ystep, ylabel):
    fig, ax = plt.subplots(figsize=(WIDTH, height))
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.xaxis.set_major_locator(MultipleLocator(xstep))
    ax.yaxis.set_major_locator(MultipleLocator(ystep))
    ax.xaxis.set_major_formatter(label)
    ax.yaxis.set_major_formatter(label)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(ylabel)
    return fig, ax


def save(fig, name):
    fig.tight_layout(pad=0.16)
    fig.savefig(name)
    plt.close(fig)


# Figure 1: erfc on the real axis
x = np.linspace(-3, 3, 4001)
fig, ax = new_axes(1.9, (-3, 3), (0, 2.05), 1, 0.5, r"$\operatorname{erfc}(x)$")
ax.plot(x, erfc(x), color=BLUE)
save(fig, "Figure_1.pdf")

# Figure 2: the test function f and its real zeros
x = np.linspace(0, 4, 8001)
z = real_zeros(f, 0, 4)
fig, ax = new_axes(2.1, (0, 4), (-1, 2.2), 0.5, 1, r"$f(x)$")
ax.axhline(0, color=AXIS, linewidth=0.6)
ax.plot(x, f(x), color=BLUE)
ax.plot(z, np.zeros_like(z), "o", ms=2.6, color=INK, zorder=3)
save(fig, "Figure_2.pdf")

# Figures 3-5: q_alpha(x) = erfc(alpha f(x)^2) for alpha = 1, 10, 100
for number, alpha in ((3, 1), (4, 10), (5, 100)):
    fig, ax = new_axes(1.6, (0, 4), (0, 1.06), 0.5, 0.5, rf"$q_{{{alpha}}}(x)$")
    ax.plot(x, erfc(alpha * f(x) ** 2), color=BLUE)
    save(fig, f"Figure_{number}.pdf")

print(f"{len(z)} real zeros of f on [0, 4]:", np.round(z, 4))
