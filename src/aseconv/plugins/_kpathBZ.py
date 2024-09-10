"""kpath on Brillouinzone drawing functions"""

from __future__ import annotations
import numpy as np


def _bzone_3d(cell: ase.cell.Cell) -> scipy.spatial.Voronoi:
    """Draw Brillouinzone 3D

    Args:
        cell

    """
    from scipy.spatial import Voronoi

    px, py, pz = np.tensordot(cell, np.mgrid[-1:2, -1:2, -1:2], axes=[0, 0])
    points = np.c_[px.ravel(), py.ravel(), pz.ravel()]
    vor = Voronoi(points)

    bz_facets = []
    bz_ridges = []
    bz_vertices = []
    for pid, rid in zip(vor.ridge_points, vor.ridge_vertices):
        # 14th points out of 27 is (0,0,0)
        if pid[0] != 13 and pid[1] != 13:
            continue
        # print(rid,vor.vertices[rid])
        bz_ridges.append(vor.vertices[np.r_[rid, [rid[0]]]])
        bz_facets.append(vor.vertices[rid])
        bz_vertices += rid

    bz_vertices = list(set(bz_vertices))
    return vor.vertices[bz_vertices], bz_ridges, bz_facets


def draw_brillouinzone(show, bzfile, atom, kp):
    import matplotlib as mpl
    import sys, os

    # if ('ss' in args.kfile):
    # 	mpl.use('Agg')

    import matplotlib.pyplot as plt
    from matplotlib.text import Annotation
    from matplotlib.patches import FancyArrowPatch
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
    from mpl_toolkits.mplot3d.proj3d import proj_transform

    class Arrow3D(FancyArrowPatch):
        def __init__(self, xs, ys, zs, *args, **kwargs):
            FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
            self._verts3d = xs, ys, zs

        def do_3d_projection(self, renderer=None):
            xs3d, ys3d, zs3d = self._verts3d
            xs, ys, zs = proj_transform(xs3d, ys3d, zs3d, self.axes.M)
            self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
            return np.min(zs)

        def draw(self, renderer):
            xs3d, ys3d, zs3d = self._verts3d
            xs, ys, zs = proj_transform(xs3d, ys3d, zs3d, self.axes.M)  # renderer.M)
            self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
            FancyArrowPatch.draw(self, renderer)

    class Annotation3D(Annotation):
        def __init__(self, s, xyz, *args, **kwargs):
            Annotation.__init__(self, s, xy=(0, 0), *args, **kwargs)
            self._verts3d = xyz

        def draw(self, renderer):
            xs3d, ys3d, zs3d = self._verts3d
            xs, ys, zs = proj_transform(xs3d, ys3d, zs3d, self.axes.M)  # renderer.M)
            self.xy = (xs, ys)
            Annotation.draw(self, renderer)

    # TODO: custom figure size
    fig = plt.figure(figsize=(6, 4), dpi=300)  # , constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")
    afontsize = 10
    cell = np.array(atom.cell.reciprocal())
    v, e, f = _bzone_3d(cell)
    # Brillouin  zone
    for xx in e:
        verts = [xx]
        ax.add_collection3d(
            Poly3DCollection(verts, facecolors="w", linewidths=0, alpha=0.4)
        )
        ax.add_collection3d(
            Line3DCollection(
                verts, colors="k", linewidths=0.8, linestyles="-", alpha=0.5
            )
        )

    kpc = kp["point_coords"]
    kpp = kp["path"]
    Tcell = cell.T
    kpoints = []

    # K-path
    for p in kpp:
        v1 = np.dot(Tcell, kpc[p[0]])
        v2 = np.dot(Tcell, kpc[p[1]])
        if p[0] not in kpoints:
            kpoints.append(p[0])
        if p[1] not in kpoints:
            kpoints.append(p[1])
        ax.add_collection3d(
            Line3DCollection([[v1, v2]], colors="r", linewidths=1.5, alpha=0.8)
        )

    # K points
    for k in kpoints:
        v = kpc[k]
        sg = k
        if k == "G" or k == "GAMMA":
            sg = "\Gamma"
        for g in ["DELTA", "SIGMA", "LAMBDA"]:
            sg = sg.replace(g, "\\" + g[0] + g[1:].lower())

        s = r"$\mathrm{\mathsf{%s}}$" % sg
        p = np.dot(Tcell, v)
        # ls="{}:{:.2f},{:.2f},{:.2f}".format(s,v[0],v[1],v[2])
        sap = []
        for vv in v:
            vs = "%.2g" % vv
            vs = vs.replace("0.", ".")
            sap.append(vs)
        ls = r"$\mathrm{\mathsf{%s}}$:(%s, %s, %s)" % (sg, sap[0], sap[1], sap[2])
        # ls=re.sub('0(?=[.])', '',ls)
        ax.scatter(p[0], p[1], p[2], marker="o", color="k", s=4, label=ls)
        ax.add_artist(
            Annotation3D(
                s,
                xyz=p,
                fontsize=afontsize,
                xytext=(4, 1),
                textcoords="offset points",
                ha="right",
                va="bottom",
            )
        )

    # Reciprocal lattices
    for id, v in enumerate(cell, 1):
        a = Arrow3D(
            [0, v[0]],
            [0, v[1]],
            [0, v[2]],
            mutation_scale=5,
            lw=0.5,
            arrowstyle="-|>",
            color="b",
            linestyle="--",
        )
        ax.add_artist(a)
        s = r"$\mathrm{\mathsf{b_%d}}$" % id
        ax.add_artist(
            Annotation3D(
                s,
                xyz=v,
                fontsize=afontsize,
                xytext=(2, 3),
                textcoords="offset points",
                ha="right",
                va="bottom",
                color="b",
            )
        )

    blen = np.linalg.norm(cell, axis=1).max()
    b1 = b2 = b3 = blen * 0.8
    ax.set_xlim(-b1, b1)
    ax.set_ylim(-b2, b2)
    ax.set_zlim(-b3, b3)
    ttl = kp.get("spacegroup_international", "")
    if ttl != "":
        ttl += f"({kp['spacegroup_number']})"
    bext = kp.get("bravais_lattice_extended", None)
    if bext != None:
        ttl += "\n[%s]" % bext

    ax.legend(
        title=ttl,
        title_fontsize="x-small",
        handletextpad=-2.0,
        markerscale=0,
        fontsize=afontsize * 0.6,
        loc="center left",
        framealpha=0.5,
    )
    zoom = 1
    ax.set_box_aspect(None, zoom=zoom)
    # plt.tight_layout()
    if show:
        plt.show()

    ax.set_axis_off()
    ofile = "{}_D{}_A{:.1f}_E{:.1f}.png".format(str(bzfile), zoom, ax.azim, ax.elev)

    fig.savefig(ofile, dpi=300, transparent=True)  # pad_inches = 0)
    plt.close()
    # TODO check convert
    try:
        os.system("convert -trim %s %s" % (ofile, ofile))
    except:
        print(sys.exc_info())
        print(" - convert trim failed, please install ImageMagick...")

    print(" - Writing '{}'...".format(ofile))
