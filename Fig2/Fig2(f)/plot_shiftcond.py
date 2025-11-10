# %%
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import cycler
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import math
import os


def get_cycle(cmap, N, use_index="auto"):
    if isinstance(cmap, str):
        if use_index == "auto":
            if cmap in [
                "Pastel1",
                "Pastel2",
                "Paired",
                "Accent",
                "Dark2",
                "Set1",
                "Set2",
                "Set3",
                "tab10",
                "tab20",
                "tab20b",
                "tab20c",
            ]:
                use_index = True
                cmap = mpl.colormaps[cmap]
                if not N:
                    N = cmap.N
            else:
                use_index = False
                cmap = mpl.colormaps[cmap]
                if not N:
                    N = 10

    if use_index == "auto":
        if cmap.N > 100:
            use_index = False
        elif isinstance(cmap, LinearSegmentedColormap):
            use_index = False
        elif isinstance(cmap, ListedColormap):
            use_index = True
    if use_index:
        ind = np.arange(int(N)) % cmap.N
        return cycler("color", cmap(ind))
    else:
        colors = cmap(np.linspace(0, 1, N))
        return cycler("color", colors)


def mpl_std_Params(x, y, cmap=None, N=None):
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    mpl.rcParams["font.family"] = "Times new roman"  # "Arial"
    mpl.rcParams["mathtext.fontset"] = "cm"
    mpl.rcParams["font.size"] = 10
    mpl.rcParams["xtick.direction"] = "in"
    mpl.rcParams["ytick.direction"] = "in"
    mpl.rcParams["ytick.labelsize"] = 8
    mpl.rcParams["xtick.labelsize"] = 8
    mpl.rcParams["figure.figsize"] = [x / 2.54, y / 2.54]  # centimeters in inches
    mpl.rcParams["lines.linewidth"] = 1
    mpl.rcParams["lines.marker"] = "None"
    mpl.rcParams["lines.markersize"] = 0.5

    if cmap:
        if type(cmap) == str:
            mpl.rcParams["axes.prop_cycle"] = get_cycle(cmap, N)
        elif type(cmap) == list:
            mpl.rcParams["axes.prop_cycle"] = cycler("color", cmap)


# %%

mpl_std_Params(
    4 * 217 / 206,
    4 * 221 / 200,
    cmap=["maroon", "lightcoral", "cornflowerblue", "midnightblue", "lightgrey"],
)

fig = plt.figure()
ax = fig.add_subplot(111)


for dirs in ["112", "333", "311", "113"]:
    filename = "./" + dirs + "/shiftcond_" + dirs + ".dat"
    sc = np.loadtxt(filename, unpack=True)
    omegas = sc[0, :]
    ax.plot(omegas, np.sum(sc[1:], axis=0))

ax.set_xlabel(r"$\omega $ (eV)")
ax.set_ylabel(r"Shift conductivity $ (\mathrm{\mu A / V^2}$)")
ax.set_xlim(1.0, 1.5)
ax.set_xticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5])

# ax2 = ax.twinx()
# ax2.plot()
plt.savefig("result.png", dpi=1000, bbox_inches="tight")
plt.savefig("result.svg", bbox_inches="tight", transparent=True)


# %%
