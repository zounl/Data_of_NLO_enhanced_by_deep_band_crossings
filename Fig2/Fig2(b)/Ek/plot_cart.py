import sys
import json
from math import pi

import numpy as np
from matplotlib import pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.colors import LogNorm


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

degeneracy = np.loadtxt("EigenValue_3.dat") - np.loadtxt("EigenValue_2.dat")
degeneracy = degeneracy[:, ::-1]
print(degeneracy.shape)

fig, ax = plt.subplots(1, 1)
im = ax.imshow(degeneracy, cmap="Blues", norm=LogNorm(vmax=1e-4, vmin=1e-1))
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
plt.title("$E_{cv}$")
plt.savefig("cart_degeneracy_colorbar.png", dpi=500)
# plt.show()
plt.close()

fig, ax = plt.subplots(1, 1)
im = ax.imshow(degeneracy, cmap="Blues", norm=LogNorm(vmax=1e-4, vmin=1e-1))
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
plt.title("$E_{cv}$")
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
plt.axis("off")
plt.savefig("cart_degeneracy.png", dpi=1000, transparent=True)
# plt.show()
plt.close()
