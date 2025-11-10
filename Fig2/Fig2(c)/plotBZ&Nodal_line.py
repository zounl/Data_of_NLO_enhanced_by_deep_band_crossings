# %%

from cmath import cos, pi, sin
from turtle import color
import numpy as np
from mayavi import mlab


filepath = "./data/"

rlv = np.array(
    [
        [0.78590898800056, 0.453744765780669, 0.305330301010679],
        [-0.78590898800056, 0.453744765780669, 0.305330301010679],
        [0.0, -0.907489531561339, 0.305330301010679],
    ]
).T


degeneracy = np.loadtxt(filepath + "EigenValue_3.dat") - np.loadtxt(
    filepath + "EigenValue_2.dat"
)


kx = np.loadtxt(filepath + "ksurface_x.txt")
ky = np.loadtxt(filepath + "ksurface_y.txt")
kz = np.loadtxt(filepath + "ksurface_z.txt")


epsilon = 1e-4

C3 = np.array(
    [
        [np.cos(2 * pi / 3), -np.sin(2 * pi / 3), 0],
        [np.sin(2 * pi / 3), np.cos(2 * pi / 3), 0],
        [0, 0, 1],
    ]
)

lx = [[] for i in range(9)]
ly = [[] for i in range(9)]
lz = [[] for i in range(9)]

imax, jmax = kx.shape
for i in range(imax):
    for j in range(jmax):
        if degeneracy[i, j] < epsilon:
            k = np.dot(rlv, [kx[i, j], ky[i, j], kz[i, j]])
            if k[1] > 0.085:
                for li in range(3):
                    lx[li].append(k[0])
                    ly[li].append(k[1])
                    lz[li].append(k[2])
                    k = np.dot(C3, k)
            elif k[1] < -0.085:
                for li in range(3):
                    lx[3 + li].append(k[0])
                    ly[3 + li].append(k[1])
                    lz[3 + li].append(k[2])
                    k = np.dot(C3, k)
            else:
                for li in range(3):
                    lx[6 + li].append(k[0])
                    ly[6 + li].append(k[1])
                    lz[6 + li].append(k[2])
                    k = np.dot(C3, k)


# %%
def get_distance(xs, ys, zs):
    pmax = len(xs)
    d = np.zeros((pmax, pmax), float)
    for i in range(pmax):
        for j in range(i):
            d[i, j] = np.linalg.norm([xs[i] - xs[j], ys[i] - ys[j], zs[i] - zs[j]])
            d[j, i] = d[i, j]
    return d


def plot_scatter_to_line(xs, ys, zs, **kwargs):
    d = get_distance(xs, ys, zs)
    nump = len(xs)
    pflist = list(range(nump))
    ips = 0
    while len(pflist) != 0:
        ps = pflist[ips]
        del pflist[ips]
        if len(pflist) == 0:
            pf = 0
        else:
            ipf = np.argmin(d[ps, pflist])
            pf = pflist[ipf]
        mlab.plot3d([xs[ps], xs[pf]], [ys[ps], ys[pf]], [zs[ps], zs[pf]], **kwargs)
        ips = ipf


def get_brillouin_zone_3d(cell):
    """
    Generate the Brillouin Zone of a given cell. The BZ is the Wigner-Seitz cell
    of the reciprocal lattice, which can be constructed by Voronoi decomposition
    to the reciprocal lattice.  A Voronoi diagram is a subdivision of the space
    into the nearest neighborhoods of a given set of points.

    https://en.wikipedia.org/wiki/Wigner%E2%80%93Seitz_cell
    https://docs.scipy.org/doc/scipy/reference/tutorial/spatial.html#voronoi-diagrams
    """

    cell = np.asarray(cell, dtype=float)
    assert cell.shape == (3, 3)

    px, py, pz = np.tensordot(cell, np.mgrid[-1:2, -1:2, -1:2], axes=[0, 0])
    points = np.c_[px.ravel(), py.ravel(), pz.ravel()]

    from scipy.spatial import Voronoi

    vor = Voronoi(points)

    bz_facets = []
    bz_ridges = []
    bz_vertices = []

    for pid, rid in zip(vor.ridge_points, vor.ridge_vertices):
        if pid[0] == 13 or pid[1] == 13:
            bz_ridges.append(vor.vertices[np.r_[rid, [rid[0]]]])
            bz_facets.append(vor.vertices[rid])
            bz_vertices += rid

    bz_vertices = list(set(bz_vertices))

    return vor.vertices[bz_vertices], bz_ridges, bz_facets


fig = mlab.figure(bgcolor=(1, 1, 1), size=(800, 800))

for li in range(9):
    if li % 3 == 0:
        plot_scatter_to_line(
            lx[li], ly[li], lz[li], tube_radius=0.005, color=(0.0, 0.0, 1.0)
        )
    # else:
    #     plot_scatter_to_line(lx[li],ly[li],lz[li],tube_radius=0.005,color=(0., 0., 1.),opacity = 0.05)

kz = np.dot(rlv, [0.5, 0.5, 0.5])[2]
kL = np.dot(C3, np.dot(C3, np.dot(rlv, [0.5, 0.0, 0.0])))
mlab.plot3d(
    [0.0, 0.0], [0.0, 0.0], [0.0, 2 * kz], tube_radius=0.006, color=(1.0, 0.0, 0.0)
)
mlab.points3d(0.0, 0.0, 0.0, scale_factor=0.03, color=(0.0, 0.0, 0.0), resolution=50)
mlab.points3d(0.0, 0.0, kz, scale_factor=0.03, color=(0.0, 0.0, 0.0), resolution=50)
mlab.points3d(0.0, 0.0, 2 * kz, scale_factor=0.03, color=(0.0, 0.0, 0.0), resolution=50)
mlab.points3d(
    kL[0], kL[1], kL[2], scale_factor=0.03, color=(0.0, 0.0, 0.0), resolution=50
)


v, e, f = get_brillouin_zone_3d(rlv.T)

for xx in e:
    mlab.plot3d(xx[:, 0], xx[:, 1], xx[:, 2], tube_radius=None, color=(0.6, 0.6, 0.6))

mlab.show()
