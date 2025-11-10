import numpy as np
from mayavi import mlab


def get_distance(xs, ys, zs):
    pmax = len(xs)
    d = np.zeros((pmax, pmax), float)
    for i in range(pmax):
        for j in range(i):
            d[i, j] = np.linalg.norm([xs[i] - xs[j], ys[i] - ys[j], zs[i] - zs[j]])
            d[j, i] = d[i, j]
    return d


def plot_scatter_to_line(xs, ys, zs):
    d = get_distance(xs, ys, zs)
    nump = len(xs)
    pflist = list(range(nump))
    ips = 0
    plist = [0]
    while len(pflist) != 0:
        ps = pflist[ips]
        del pflist[ips]
        if len(pflist) == 0:
            pf = 0
        else:
            ipf = np.argmin(d[ps, pflist])
            pf = pflist[ipf]
        plist.append(pf)
        ips = ipf
    return plist


# Define the meshgrid for x and y
x, y = np.mgrid[-1:1:500j, -1:1:500j]

# Compute the corresponding z values for each surface
z1 = np.sqrt(y**2 + x**2)
z2 = -np.sqrt(y**2 + x**2)
c = y / (x**2 + y**2)
z3 = -2 + 0.25 * y
zcut = 0.5 + 0.25 * y

# Only plot z1 where z1 < 2
z1[z1 >= 1] = np.nan

# Only plot z2 where z2 > -2
z2[z2 <= -1] = np.nan

# Create the figure and plot the surfaces
fig = mlab.figure(bgcolor=(1, 1, 1), fgcolor=(1, 1, 1), size=(800, 800))
mesh = mlab.mesh(
    x,
    y,
    z1,
    scalars=c,
    colormap="coolwarm",
    vmin=-3,
    vmax=3,
    representation="surface",
    opacity=0.9,
)
mlab.mesh(x, y, z2, scalars=c, colormap="coolwarm", vmin=-3, vmax=3, opacity=0.9)
mlab.surf(x, y, z3, color=(230 / 255, 230 / 255, 250 / 255))

# Adjust material properties for the surface
# mesh.actor.property.specular = 0  # Set specular reflection to 0, reducing the smoothness
# mesh.actor.property.specular_power = 0.5  # Set specular power to 0, further reducing the smoothness

# Create boolean mask for the condition z2 = zcut
mask = np.isclose(z1, zcut, rtol=5.0e-4)

# Apply the mask to extract the corresponding values
x0 = x[mask]
y0 = y[mask]
z0 = z1[mask]
plist = plot_scatter_to_line(x0, y0, z0)
plist.append(0)
mlab.plot3d(x0[plist], y0[plist], z0[plist], color=(0, 0, 0), tube_radius=0.015)
z0 = z3[mask]
c0 = c[mask]
mlab.plot3d(
    x0[plist],
    y0[plist],
    z0[plist],
    c0[plist],
    tube_radius=0.015,
    colormap="coolwarm",
    vmin=-3,
    vmax=3,
)
# mlab.outline(color=(0.5, 0.5, 0.5), extent=(-1, 1, -1, 1, -2.5, 1))


mlab.show()
