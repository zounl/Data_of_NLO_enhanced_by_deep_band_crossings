# %%
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import cycler
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import os
from cmath import pi
import matplotlib.ticker as ticker


def read_file(file_path):
    with open(file_path, "r") as file:
        lines = file.readlines()

    names = []
    xs = []
    ys = []
    for i in range(len(lines)):
        line = lines[i]
        if line.startswith("calculation finish"):
            break
        if i % 3 == 0:
            names.append(line.strip())
        elif i % 3 == 1:
            xs.append([float(x) for x in line.strip().split()])
        elif i % 3 == 2:
            ys.append([float(x) for x in line.strip().split()])
    return names, xs, ys


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


def set_axlabel(ax, name, **kwargs):
    ax.set_ylabel(r"$\sigma^{xxy} (\rm{\mu A / V^2}$)", labelpad=-2)
    return "sigma"


def read_folder(file_folder):
    file_list = os.listdir(file_folder)
    dat_files = [file for file in file_list if file.endswith(".dat")]
    print("reading: ", dat_files)
    names = []
    xs = []
    ys = []
    for file in dat_files:
        names_tmp, xs_tmp, ys_tmp = read_file(file_folder + file)
        names += names_tmp
        xs += xs_tmp
        ys += ys_tmp
    return names, xs, ys


def eV_to_meV(x, pos):
    return f"{x * 1000:.0f}"


# %% get kp
kp_para = np.loadtxt("kp_parameter")
eps_calc = np.abs(kp_para[:, 1])
eps_typeI = np.where(eps_calc < 1, eps_calc, np.nan)
eps_typeII = np.where(eps_calc > 1, eps_calc, np.nan)
sc_const = 4 * pi * np.abs(kp_para[:, 0])  # 2*kz(2pi)
ntheta = 4000

eps, theta = np.meshgrid(
    eps_typeI,
    np.linspace(start=0, stop=2 * pi, num=ntheta),
)


kx = 1 / (1 - eps**2) * (np.sin(theta) + eps)
normk = eps * kx + 1
norm_dvg_f = np.sqrt(1 - 2 * eps * kx / normk + eps**2)
dl = np.sqrt(np.sin(theta) ** 2 * (1 - eps**2) + np.cos(theta) ** 2) / (1 - eps**2)
sigmas1 = sc_const * np.sum(kx / normk**2 / norm_dvg_f * dl, axis=0) / ntheta * 2 * pi


eps, theta = np.meshgrid(
    eps_typeII,
    np.linspace(start=-pi / 2, stop=pi / 2, num=ntheta + 2)[1:-1],
)
kx = (1 / np.cos(theta) - eps) / (eps**2 - 1)
normk = eps * kx + 1
norm_dvg_f = np.sqrt(1 - 2 * eps * kx / normk + eps**2)
dl = np.sqrt(
    1 / np.cos(theta) ** 4 / (eps**2 - 1)
    + (np.tan(theta) / np.cos(theta) / (eps**2 - 1)) ** 2
)
sigmas2 = sc_const * np.sum(kx / normk**2 / norm_dvg_f * dl, axis=0) / ntheta * pi

eps, theta = np.meshgrid(
    eps_typeII,
    np.linspace(start=pi / 2, stop=3 * pi / 2, num=ntheta + 2)[1:-1],
)
kx = (1 / np.cos(theta) - eps) / (eps**2 - 1)
normk = -eps * kx - 1
norm_dvg_f = np.sqrt(1 + 2 * eps * kx / normk + eps**2)
dl = np.sqrt(
    1 / np.cos(theta) ** 4 / (eps**2 - 1)
    + (np.tan(theta) / np.cos(theta) / (eps**2 - 1)) ** 2
)
sigmas3 = sc_const * np.sum(kx / normk**2 / norm_dvg_f * dl, axis=0) / ntheta * pi

# %%

file_folder = "./"
names, xs, ys = read_folder(file_folder)

# initialization mpl Params
cmap = ["lightseagreen", "lightcoral", "royalblue"]
mpl_std_Params(3, 3.5, cmap=cmap)

fig = plt.figure()
ax = fig.add_subplot(111)

# label and ticks
ax1_label = ""
ax2_label = ""
ax.set_xlabel(r"$\hbar \omega - \epsilon_0$ (meV)", fontsize=10, color="k")
ax.tick_params(axis="both", labelsize=10, color="k", labelcolor="k")
formatter = ticker.FuncFormatter(eV_to_meV)
plt.gca().xaxis.set_major_formatter(formatter)

ax.set_ylim(-1.5, 1.5)
# ax.set_yscale("symlog")
# ax.set_ylim(-0.3, 0.3)
ax.set_xlim(-0.025, 0.025)
plot_thread = 0.001

for i in [1, 2, 3]:
    print(i, eps_calc[i])
    xsp = np.array(xs[i])
    ysp = np.array(ys[i])
    if ax1_label == "":
        ax1_label = set_axlabel(ax, names[i], fontsize=10)
        ax.plot(
            xsp[((xsp < -plot_thread) | (xsp > plot_thread))],
            ysp[((xsp < -plot_thread) | (xsp > plot_thread))],
            label=names[i],
        )
    else:
        ax.plot(
            xsp[((xsp < -plot_thread) | (xsp > plot_thread))],
            ysp[((xsp < -plot_thread) | (xsp > plot_thread))],
            label=names[i],
        )

ax.set_prop_cycle(None)
for i in [1, 2, 3]:
    if eps_calc[i] < 1:
        if i != 2:
            ax.plot(xs[i], sigmas1[i] * np.sign(xs[i]), linestyle="--")
        else:
            ax.plot(xs[i], sigmas1[i] * np.sign(xs[i]), linestyle="--")
    else:
        ax.plot(xs[i], (sigmas3[i] + sigmas2[i]) * np.sign(xs[i]), linestyle="--")

# ax.hlines(sigmas1, -0.1, 0.1, color="k", linestyles="dashed")
# ax.legend(loc="upper left", fontsize=6, handlelength=0.5)
plt.savefig("result.png", dpi=1000, bbox_inches="tight")
plt.savefig("result.pdf", bbox_inches="tight")
plt.close()

# %%
