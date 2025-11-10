# %%
from math import sqrt
from matplotlib.pyplot import legend
import numpy as np
from mayavi import mlab
from tvtk.tools import visual


def Arrow_From_A_to_B(x1, y1, z1, x2, y2, z2):
    ar1 = visual.arrow(x=x1, y=y1, z=z1, color=(0, 0, 0), resolution=100)
    ar1.length_cone = 0.3
    ar1.radius_cone = 0.1

    arrow_length = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
    ar1.actor.scale = [arrow_length, arrow_length, arrow_length]
    ar1.pos = ar1.pos / arrow_length
    ar1.axis = [x2 - x1, y2 - y1, z2 - z1]
    return ar1


def plot_dotted_line(xs, ys, zs, **kwargs):
    for i in range(0, len(xs), 2):
        mlab.plot3d(xs[i : i + 2], ys[i : i + 2], zs[i : i + 2], **kwargs)


fig = mlab.figure(bgcolor=(1, 1, 1), fgcolor=(0.3, 0.3, 0.3), size=(800, 800))

o = np.array([0, 0, 0])
a1 = np.array([2.11533318, 1.22128818, 3.62985997])
a2 = np.array([-2.11533318, 1.22128818, 3.62985997])
a3 = np.array([0.0, -2.44257636, 3.62985997])
lv = np.array([a1, a2, a3]).T

# primitive cell
# line_radius = 0.03
# mlab.plot3d([o[0],a1[0]],[o[1],a1[1]],[o[2],a1[2]],tube_radius=line_radius,color=(0.6, 0.6, 0.6))
# mlab.plot3d([o[0],a2[0]],[o[1],a2[1]],[o[2],a2[2]],tube_radius=line_radius,color=(0.6, 0.6, 0.6))
# mlab.plot3d([o[0],a3[0]],[o[1],a3[1]],[o[2],a3[2]],tube_radius=line_radius,color=(0.6, 0.6, 0.6))
#
# mlab.plot3d([a1[0],(a1+a2)[0]],[a1[1],(a1+a2)[1]],[a1[2],(a1+a2)[2]],tube_radius=line_radius,color=(0.6, 0.6, 0.6))
# mlab.plot3d([a1[0],(a1+a3)[0]],[a1[1],(a1+a3)[1]],[a1[2],(a1+a3)[2]],tube_radius=line_radius,color=(0.6, 0.6, 0.6))
# mlab.plot3d([a2[0],(a1+a2)[0]],[a2[1],(a1+a2)[1]],[a2[2],(a1+a2)[2]],tube_radius=line_radius,color=(0.6, 0.6, 0.6))
# mlab.plot3d([a2[0],(a3+a2)[0]],[a2[1],(a3+a2)[1]],[a2[2],(a3+a2)[2]],tube_radius=line_radius,color=(0.6, 0.6, 0.6))
# mlab.plot3d([a3[0],(a3+a2)[0]],[a3[1],(a3+a2)[1]],[a3[2],(a3+a2)[2]],tube_radius=line_radius,color=(0.6, 0.6, 0.6))
# mlab.plot3d([a3[0],(a3+a1)[0]],[a3[1],(a3+a1)[1]],[a3[2],(a3+a1)[2]],tube_radius=line_radius,color=(0.6, 0.6, 0.6))
#
# mlab.plot3d([(a1+a2)[0],(a1+a2+a3)[0]],[(a1+a2)[1],(a1+a2+a3)[1]],[(a1+a2)[2],(a1+a2+a3)[2]],tube_radius=line_radius,color=(0.6, 0.6, 0.6))
# mlab.plot3d([(a3+a2)[0],(a1+a2+a3)[0]],[(a3+a2)[1],(a1+a2+a3)[1]],[(a3+a2)[2],(a1+a2+a3)[2]],tube_radius=line_radius,color=(0.6, 0.6, 0.6))
# mlab.plot3d([(a1+a3)[0],(a1+a2+a3)[0]],[(a1+a3)[1],(a1+a2+a3)[1]],[(a1+a3)[2],(a1+a2+a3)[2]],tube_radius=line_radius,color=(0.6, 0.6, 0.6))
#
# mlab.plot3d([(a1+a2)[0],(a1+a2+a3)[0]],[(a1+a2)[1],(a1+a2+a3)[1]],[(a1+a2)[2],(a1+a2+a3)[2]],tube_radius=line_radius,color=(0.6, 0.6, 0.6))


# mirror plane

# y,z = np.mgrid[-1.5:1.5:100j, 0:10.88957991:100j]
# yy,zz = np.mgrid[-3:3:100j, 0:10.88957991:100j]
# x1 = 0*y
# x2 = sqrt(3)*y
# x3 = -sqrt(3)*y

# mlab.mesh(x1,yy,zz,opacity=0.4,color=(1,1,1))
# mlab.mesh(x2,y,z,opacity=0.4,color=(1,1,1))
# mlab.mesh(x3,y,z,opacity=0.4,color=(1,1,1))
# mt = mlab.text3d(0,3,10.88957991,'Mirror',line_width = 0.2,scale = 0.5,orient_to_camera = False,orientation = [90.,0.,-90.])

# C3axis

N = 30
plot_dotted_line(
    np.linspace(o[0], (a1 + a2 + a3)[0], num=N),
    np.linspace(o[1], (a1 + a2 + a3)[1], num=N),
    np.linspace(o[2], (a1 + a2 + a3)[2], num=N),
    tube_radius=None,
    color=(0.6, 0.6, 0.6),
)

# mt = mlab.text3d(
#     0,
#     -0.5,
#     11.88957991,
#     "C3",
#     line_width=0.2,
#     scale=0.5,
#     orient_to_camera=False,
#     orientation=[90.0, 0.0, -90.0],
# )
visual.set_viewer(fig)
Arrow_From_A_to_B(0, 0, 0.0, -2, 0, 0)
Arrow_From_A_to_B(0, 0, 0.0, 0, 2, 0)
Arrow_From_A_to_B(0, 0, 0.0, 0, 0, 2)


Ge_size = 0.8
Te_size = 0.8
Ge_color = (1.0, 0, 0.0)
Te_color = (0.0, 0.0, 1.0)

ks = np.array(
    [
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 1, 0],
        [1, 0, 1],
        [0, 1, 1],
        [1, 1, 1],
        [2, 0, 0],
        [0, 2, 0],
        [0, 0, 2],
        [1, 1, -1],
        [1, -1, 1],
        [-1, 1, 1],
    ]
).T
points = np.dot(lv, ks)
mlab.points3d(
    points[0, :],
    points[1, :],
    points[2, :],
    scale_factor=Ge_size,
    color=Ge_color,
    resolution=100,
)


ks2 = np.array(
    [
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [-1, 0, 0],
        [0, -1, 0],
        [0, 0, -1],
        [1, -1, 0],
        [0, 1, -1],
        [-1, 0, 1],
        [-1, 1, 0],
        [0, -1, 1],
        [1, 0, -1],
    ]
) + np.array([0.53088])
points = np.dot(lv, ks2.T)
mlab.points3d(
    points[0, :],
    points[1, :],
    points[2, :],
    scale_factor=Te_size,
    color=Te_color,
    resolution=100,
)

# conventional cell
N = 20
ks = np.array(
    [
        [0, 0, 0],
        [1, 1, 1],
        [2, 0, 0],
        [0, 2, 0],
        [0, 0, 2],
        [1, 1, -1],
        [1, -1, 1],
        [-1, 1, 1],
    ]
).T
points = np.dot(lv, ks)

lines = [
    [0, 5],
    [5, 2],
    [2, 6],
    [6, 0],
    [7, 3],
    [3, 1],
    [1, 4],
    [4, 7],
    [0, 7],
    [6, 4],
    [2, 1],
    [5, 3],
]
for line in lines:
    start = points[:, line[0]]
    stop = points[:, line[1]]
    plot_dotted_line(
        np.linspace(start[0], stop[0], num=N),
        np.linspace(start[1], stop[1], num=N),
        np.linspace(start[2], stop[2], num=N),
        tube_radius=None,
        color=(1, 0.6, 0.6),
    )

# ks = np.array(
#     [
#         [0, 0, 0],
#         [1, 1, 1],
#         [2, 0, 0],
#         [0, 2, 0],
#         [0, 0, 2],
#         [1, 1, -1],
#         [1, -1, 1],
#         [-1, 1, 1],
#     ]
# ) + np.array([0.03088])
# points = np.dot(lv, ks.T)
# lines = [
#     [0, 5],
#     [5, 2],
#     [2, 6],
#     [6, 0],
#     [7, 3],
#     [3, 1],
#     [1, 4],
#     [4, 7],
#     [0, 7],
#     [6, 4],
#     [2, 1],
#     [5, 3],
# ]
# for line in lines:
#     start = points[:, line[0]]
#     stop = points[:, line[1]]
#     plot_dotted_line(
#         np.linspace(start[0], stop[0], num=N),
#         np.linspace(start[1], stop[1], num=N),
#         np.linspace(start[2], stop[2], num=N),
#         tube_radius=None,
#         color=(0.6, 0.6, 1),
#     )


# mirror plane边框

# xs = np.array([0 for i in range(20)])
# ys = np.linspace(-3,3,20)
# zs = np.array([0 for i in range(20)])
# plot_dotted_line(xs,ys,zs,tube_radius=None,color=(0.3, 0.3, 0.3))
# zs = np.array([10.88957991 for i in range(20)])
# plot_dotted_line(xs,ys,zs,tube_radius=None,color=(0.3, 0.3, 0.3))
#
# xs = np.array([0 for i in range(30)])
# ys = np.array([-3 for i in range(30)])
# zs = np.linspace(0,10.88957991,30)
# plot_dotted_line(xs,ys,zs,tube_radius=None,color=(0.3, 0.3, 0.3))
# ys = np.array([3 for i in range(30)])
# plot_dotted_line(xs,ys,zs,tube_radius=None,color=(0.3, 0.3, 0.3))


# other mirror
# ys = np.linspace(-1.5,1.5,20)
# xs = sqrt(3)*ys
# zs = np.array([0 for i in range(20)])
# plot_dotted_line(xs,ys,zs,tube_radius=None,color=(0.3, 0.3, 0.3))
# plot_dotted_line(-xs,ys,zs,tube_radius=None,color=(0.3, 0.3, 0.3))
# zs = np.array([10.88957991 for i in range(20)])
# plot_dotted_line(xs,ys,zs,tube_radius=None,color=(0.3, 0.3, 0.3))
# plot_dotted_line(-xs,ys,zs,tube_radius=None,color=(0.3, 0.3, 0.3))
#
# xs = np.array([-1.5*sqrt(3) for i in range(30)])
# ys = np.array([-1.5 for i in range(30)])
# zs = np.linspace(0,10.88957991,30)
# plot_dotted_line(xs,ys,zs,tube_radius=None,color=(0.3, 0.3, 0.3))
# plot_dotted_line(-xs,ys,zs,tube_radius=None,color=(0.3, 0.3, 0.3))
# plot_dotted_line(xs,-ys,zs,tube_radius=None,color=(0.3, 0.3, 0.3))
# plot_dotted_line(-xs,-ys,zs,tube_radius=None,color=(0.3, 0.3, 0.3))

visual.set_viewer(fig)
mlab.show()
# %%
# legend

# fig = mlab.figure(
#    bgcolor=(1, 1, 1),
#    fgcolor=(0.3,0.3,0.3),
#    size=(800, 800)
# )
#
# Ge_size = 0.8
# Te_size = 0.8
# Ge_color = (1., 0, 0.)
# Te_color = (0.,0.,1.)
# mlab.points3d(0,-4,10.88957991,scale_factor=Ge_size,color=Ge_color,resolution=100)
# mlab.points3d(0,-4,8.88957991,scale_factor=Te_size,color=Te_color,resolution=100)
#
# mlab.show()
# %%
