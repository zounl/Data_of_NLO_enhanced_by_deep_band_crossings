# %%
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import cycler
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
import math
from matplotlib.collections import LineCollection


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
    mpl.rcParams["ytick.labelsize"] = 8
    mpl.rcParams["xtick.labelsize"] = 8
    mpl.rcParams["xtick.direction"] = "in"
    mpl.rcParams["ytick.direction"] = "in"
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

filepath = ""
rlv = np.array(
    [
        [0.78590898800056, 0.453744765780669, 0.305330301010679],
        [-0.78590898800056, 0.453744765780669, 0.305330301010679],
        [0.0, -0.907489531561339, 0.305330301010679],
    ]
).T


degeneracy = np.loadtxt(filepath + "Ecp.dat") - np.loadtxt(filepath + "Ec.dat")
gap = np.loadtxt(filepath + "Ec.dat") - np.loadtxt(filepath + "Ev.dat")
epsilon_c = np.abs(np.loadtxt(filepath + "epsilon_c.dat")) / 2
epsilon_v = np.abs(np.loadtxt(filepath + "epsilon_v.dat"))


kx = np.loadtxt(filepath + "kx.dat")
ky = np.loadtxt(filepath + "ky.dat")
kz = np.loadtxt(filepath + "kz.dat")
# %%
line_sign = np.zeros(degeneracy.shape)
slice_vec = np.dot(rlv, [1, 1, 1])

omega_dirac = []
epec_dirac = []
k_dirac = []
deg_thres = 5e-4

for i in range(degeneracy.shape[0]):
    for j in range(degeneracy.shape[1]):
        if degeneracy[i, j] < deg_thres:
            k = np.array([kx[i, j], ky[i, j], kz[i, j]])
            k_car = np.dot(rlv, [kx[i, j], ky[i, j], kz[i, j]])
            slice_k = np.dot(k_car, slice_vec) / np.dot(slice_vec, slice_vec)

            if abs(k[2] - k[1]) > 1e-3:
                k_dirac.append(slice_k)
                omega_dirac.append(gap[i, j])
                epec_dirac.append(epsilon_v[i, j] / epsilon_c[i, j])

filename = "shiftcond_112.dat"
sc = np.loadtxt(filename, unpack=True)
omegas = sc[0, :]
batch = 1
sc_batch = np.zeros((sc.shape[0] // batch, sc.shape[1]))
for i in range(1, sc.shape[0], batch):
    sc_batch[i // batch - 1, :] = np.sum(sc[i : i + batch, :], axis=0)

nsc = np.linalg.norm(sc_batch, axis=1)
snsc = sorted(nsc)
mpl_std_Params(6, 17 / 4, cmap="Set2")
# plt.style.use("classic")
plt.plot(nsc)
# plt.autoscale(enable=True, axis="both", tight=True)
plt.savefig(filename[:-4] + "_nsc.png")
plt.close()

# %%

mpl_std_Params(6.5, 4.5, cmap="Set2")
fig = plt.figure()
ax = fig.add_subplot(111)


ydis = 0.02

j = 0

for i in range(0, sc_batch.shape[0]):
    ax.plot(omegas, sc_batch[i, :] / 2000 + ydis * j, lw=1, c="royalblue")
    j = j + 1

ax.plot([1.03, 1.03], [0.35, 0.35 + 100 / 2000], lw=0.5, c="black")

scatter = ax.scatter(
    omega_dirac, k_dirac, c=epec_dirac, zorder=-10, vmin=1, vmax=9, cmap="Reds_r"
)
cbar = plt.colorbar(scatter, ax=ax, pad=0.15)
cbar.set_ticks([1, 3, 5, 7, 9])
cbar.ax.set_ylabel(r"$v_{\mathrm{v}} /v_{\mathrm{c}}$", rotation=270, labelpad=13)


ax.set_xlabel(r"$\omega $ (eV)")
ax.set_ylabel(r"$\sigma^{xxy}$/slice $(\mathrm{\mu A / V^2}$)")
y_positions = [0.3, 0.4, 0.5, 0.6, 0.7]
tick_labels = [0, 200, 400, 600, 800]
# ax.set_yticks(y_positions, tick_labels)
ax.set_yticks([])


ax.set_xticks([1.0, 1.1, 1.2, 1.3, 1.4, 1.5])
ax.set_xlim(1, 1.5)
ax.set_ylim(0.3, 0.7)


plt.savefig(filename + ".png", dpi=1000, bbox_inches="tight")
plt.savefig(filename + ".svg", bbox_inches="tight", transparent=True)
plt.close()

# %%
