import sys
import json
from math import pi

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.colors import SymLogNorm

with open("config", "r") as json_f:
    json_dict = json.load(json_f)
nms = json_dict["nms"]
alpha = json_dict["alpha"]
beta = json_dict["beta"]
gamma = json_dict["gamma"]

# frac: (a, a, b)
# cart = (a, a, b) @ rlv =

rlv = np.array(
    [
        [0.78590898800056, 0.453744765780669, 0.305330301010679],
        [-0.78590898800056, 0.453744765780669, 0.305330301010679],
        [0.0, -0.907489531561339, 0.305330301010679],
    ]
)
a = np.array([0, 0, 0] @ rlv)
b = np.array([1, 1, 0] @ rlv)
c = np.array([0, 0, 1] @ rlv)
d = np.array([1, 1, 1] @ rlv)

theta = (
    np.arccos(np.dot(b - a, c - a) / np.linalg.norm(b - a) / np.linalg.norm(c - a))
    * 180
    / pi
)
print(theta)
print(np.cos((theta - 90) * pi / 180))
print(np.linalg.norm(c - a))
print(np.linalg.norm(b - a))


def sym_log_transform(I, linthresh=1):
    """对数据进行对称对数变换."""
    return np.sign(I) * np.log10(np.maximum(np.abs(I), linthresh))


for i, nm in enumerate(nms):
    I = np.loadtxt("I_surface_%d%d%d_%s.dat" % (alpha, beta, gamma, nm))
    I = I[:, ::-1]
    print(I.shape)

    fig, ax = plt.subplots(1, 1)

    # 使用 SymLogNorm 显示正负双对数坐标
    linthresh = 1  # 线性范围的阈值
    vmin, vmax = -1e2, 1e2
    norm = SymLogNorm(linthresh=linthresh, vmin=vmin, vmax=vmax, base=10)

    im = ax.imshow(I, norm=norm, cmap="seismic")
    im.set_transform(
        mtransforms.Affine2D()
        .scale(
            np.linalg.norm(b - a) / np.sqrt(2) * np.cos((theta - 90) * pi / 180),
            np.linalg.norm(c - a),
        )
        .skew_deg(0, theta - 90)
        + ax.transData
    )
    plt.xlim(0, 1500)
    plt.ylim(0, 3000)
    plt.colorbar(im)
    plt.title("$I_{%s}^{%s%s%s}$" % (nm, alpha, beta, gamma))
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.axis("off")
    plt.savefig(
        "cart_I_surface_%s%s%s_%s_log.png" % (alpha, beta, gamma, nm),
        dpi=1000,
        transparent=True,
    )
    # plt.show()
    plt.close()

    fig, ax = plt.subplots(1, 1)
    im = ax.imshow(I, vmin=-2, vmax=2, cmap="RdBu")
    im.set_transform(
        mtransforms.Affine2D()
        .scale(
            np.linalg.norm(b - a) / np.sqrt(2) * np.cos((theta - 90) * pi / 180),
            np.linalg.norm(c - a),
        )
        .skew_deg(0, theta - 90)
        + ax.transData
    )
    plt.xlim(0, 1500)
    plt.ylim(0, 3000)
    plt.colorbar(im)
    plt.title("$I_{%s}^{%s%s%s}$" % (nm, alpha, beta, gamma))
    plt.savefig(
        "cart_I_surface_%s%s%s_%s_colorbar.png" % (alpha, beta, gamma, nm), dpi=500
    )
    # plt.show()
    plt.close()

    fig, ax = plt.subplots(1, 1)
    im = ax.imshow(I, vmin=-2, vmax=2, cmap="RdBu")
    im.set_transform(
        mtransforms.Affine2D()
        .scale(
            np.linalg.norm(b - a) / np.sqrt(2) * np.cos((theta - 90) * pi / 180),
            np.linalg.norm(c - a),
        )
        .skew_deg(0, theta - 90)
        + ax.transData
    )
    plt.xlim(0, 1500)
    plt.ylim(0, 3000)
    plt.colorbar(im)
    plt.title("$I_{%s}^{%s%s%s}$" % (nm, alpha, beta, gamma))
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    plt.axis("off")
    plt.savefig(
        "cart_I_surface_%s%s%s_%s.png" % (alpha, beta, gamma, nm),
        dpi=1000,
        transparent=True,
    )
    # plt.show()
    plt.close()
