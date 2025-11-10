# %%
from cmath import pi
from re import S
import numpy as np
import matplotlib.pyplot as plt

neps1 = 100
ntheta = 2000
eps_stop = 5


eps, theta = np.meshgrid(
    np.linspace(start=0, stop=1, num=neps1)[:-1],
    np.linspace(start=0, stop=2 * pi, num=ntheta),
)

kx = 1 / (1 - eps**2) * (np.sin(theta) + eps)
normk = eps * kx + 1
norm_dvg_f = np.sqrt(1 - 2 * eps * kx / normk + eps**2)
dl = np.sqrt(np.sin(theta) ** 2 * (1 - eps**2) + np.cos(theta) ** 2) / (1 - eps**2)
sigmas1 = np.sum(kx / normk**2 / norm_dvg_f * dl, axis=0) / ntheta * 2 * pi

# %%
neps = (eps_stop - 1) * neps1
ntheta = 1000
# theta_array = np.concatenate((np.linspace(start=-pi/2, stop=pi/2, num=ntheta), np.linspace(start=-pi/2, stop=pi/2, num=ntheta)))
eps, theta = np.meshgrid(
    np.linspace(start=1, stop=eps_stop, num=neps)[1:],
    np.linspace(start=-pi / 2, stop=pi / 2, num=ntheta + 2)[1:-1],
)
kx = (1 / np.cos(theta) - eps) / (eps**2 - 1)
normk = eps * kx + 1
norm_dvg_f = np.sqrt(1 - 2 * eps * kx / normk + eps**2)
dl = np.sqrt(
    1 / np.cos(theta) ** 4 / (eps**2 - 1)
    + (np.tan(theta) / np.cos(theta) / (eps**2 - 1)) ** 2
)
sigmas2 = np.sum(kx / normk**2 / norm_dvg_f * dl, axis=0) / ntheta * pi

eps, theta = np.meshgrid(
    np.linspace(start=1, stop=eps_stop, num=neps)[1:],
    np.linspace(start=pi / 2, stop=3 * pi / 2, num=ntheta + 2)[1:-1],
)
kx = (1 / np.cos(theta) - eps) / (eps**2 - 1)
normk = -eps * kx - 1
norm_dvg_f = np.sqrt(1 + 2 * eps * kx / normk + eps**2)
dl = np.sqrt(
    1 / np.cos(theta) ** 4 / (eps**2 - 1)
    + (np.tan(theta) / np.cos(theta) / (eps**2 - 1)) ** 2
)
sigmas3 = np.sum(kx / normk**2 / norm_dvg_f * dl, axis=0) / ntheta * pi

# %%
cm = 1 / 2.54  # centimeters in inches
fig = plt.figure(figsize=(4 * cm, 5.2 * cm))
plt.rcParams.update(
    {"font.size": 10, "font.family": "Times new roman", "mathtext.fontset": "cm"}
)
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
ax = fig.add_subplot(111)

ax.plot(
    np.linspace(start=0, stop=1, num=neps1)[:-1],
    sigmas1,
    label="type-1",
    c="lightcoral",
    lw=1,
)
ax.plot(
    np.linspace(start=1, stop=eps_stop, num=neps)[1:],
    sigmas2 + sigmas3,
    label="type-2",
    c="lightcoral",
    lw=1,
)
ax.axvline(1.0, color="grey", linestyle="--", lw=1)

ax.set_ylabel(r"$\tilde{\sigma}^{xxy}$/C ", fontsize=10)
ax.yaxis.set_label_coords(-0.15, 0.5)
ax.xaxis.set_label_coords(0.5, -0.15)
ax.set_xlabel(
    r"Velocity ratio $v_\mathrm{v} / v_\mathrm{c} $", fontsize=10, labelpad=-2
)
ax.tick_params(axis="both", labelsize=10)
eps_stop = 4
plt.xlim(0, eps_stop)
plt.xticks(range(0, eps_stop + 1, 1))
plt.ylim(-7, 10)
plt.savefig("integrate_kp_withoutdelta.png", dpi=1000, bbox_inches="tight")
plt.savefig("integrate_kp_withoutdelta.pdf", bbox_inches="tight", transparent=True)
# %%
