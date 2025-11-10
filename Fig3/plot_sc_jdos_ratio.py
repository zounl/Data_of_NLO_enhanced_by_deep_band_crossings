# %%
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import re
from matplotlib.pyplot import cycler
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import pandas as pd
from adjustText import adjust_text

# %%
df = pd.read_csv("filtered_sc_peak.csv")

peak_form = df["Formular"].values
scpeaks = df["peak"].values
omegas = df["omegas"].values
ratios = df["mirror_ratios"].values
rmirror = df["rmirror"].values

# %%


def formula_to_latex(formula):
    """将化学式转换为latex格式"""

    # 匹配所有数字
    pattern = re.compile(r"(\d+)")

    # 将每个匹配到的数字替换为下标
    latex = pattern.sub(r"$_{\1}$", formula)

    return latex


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
    mpl.rcParams["figure.figsize"] = [x / 2.54, y / 2.54]  # centimeters in inches
    mpl.rcParams["lines.linewidth"] = 1
    mpl.rcParams["lines.marker"] = "None"

    if cmap:
        if type(cmap) == str:
            mpl.rcParams["axes.prop_cycle"] = get_cycle(cmap, N)
        elif type(cmap) == list:
            mpl.rcParams["axes.prop_cycle"] = cycler("color", cmap)


cm = 1 / 2.54  # centimeters in inches

mpl_std_Params(8.7, 6.6)
fig1, ax1 = plt.subplots()

plt.rcParams.update(
    {"font.size": 10, "font.family": "Times new roman", "mathtext.fontset": "cm"}
)
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"


im = ax1.scatter(
    omegas[rmirror != np.nan],
    scpeaks[rmirror != np.nan],
    s=1 + 15 * ratios[rmirror != np.nan],
    c=ratios[rmirror != np.nan],
)
# fig1.colorbar(im, ax=ax1)
# fig1.colorbar()
# norm = plt.Normalize(ratios.min(), ratios.max())
texts = []
text_peak_threa = 200
for i in range(omegas.size):
    if rmirror[i] != np.nan:
        if scpeaks[i] > text_peak_threa:
            texts.append(
                ax1.text(
                    omegas[i], scpeaks[i], formula_to_latex(peak_form[i]), fontsize=6
                )
            )
            # texts.append(plt.text(omegas[i], scpeaks[i], formula_to_latex(peak_form[i]),fontsize = 6, color=plt.cm.viridis(norm(ratios[i]))))
ax1.set_xlabel("Frequency (eV)")
ax1.set_ylabel(r"Peak value ($\rm{\mu A / V^2}$)")
# ax1.set_ylim(150,650)
ax1.set_xscale("log")
ax1.set_xticks([1, 2, 3, 4, 5, 6])
ax1.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())


# plt.xlim(0,4.5)
adjust_text(texts, arrowprops=dict(arrowstyle="-", color="grey", lw=0.1))
plt.savefig("peak_filter_omega_ratio.png", bbox_inches="tight")
plt.savefig("peak_filter_omega_ratio.svg", bbox_inches="tight")
plt.close()
# %%
