"""Illustration of how S2 assigns a shared vertex (from vgrid/dggs/s2.py)."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from vgrid.conversion.latlon2dggs import latlon2s2
from vgrid.conversion.dggs2geo.s22geo import s22geo
from vgrid.dggs import s2

OUT = Path(__file__).with_name("s2_shared_vertex.png")

FACE_COLORS = {
    0: "#4C78A8",
    1: "#F58518",
    2: "#54A24B",
    3: "#E45756",
    4: "#B279A2",
    5: "#72B7B2",
}
FACE_NAMES = {
    0: "+X  face 0",
    1: "+Y  face 1",
    2: "+Z  face 2",
    3: "−X  face 3",
    4: "−Y  face 4",
    5: "−Z  face 5",
}


def cube_face_verts(face, s=1.02):
    """Slightly inflated cube-face quads for a 3D drawing."""
    if face == 0:  # +X
        return [(s, -1, -1), (s, 1, -1), (s, 1, 1), (s, -1, 1)]
    if face == 1:  # +Y
        return [(-1, s, -1), (1, s, -1), (1, s, 1), (-1, s, 1)]
    if face == 2:  # +Z
        return [(-1, -1, s), (1, -1, s), (1, 1, s), (-1, 1, s)]
    if face == 3:  # -X
        return [(-s, -1, -1), (-s, -1, 1), (-s, 1, 1), (-s, 1, -1)]
    if face == 4:  # -Y
        return [(-1, -s, -1), (-1, -s, 1), (1, -s, 1), (1, -s, -1)]
    return [(-1, -1, -s), (1, -1, -s), (1, 1, -s), (-1, 1, -s)]


def draw_cube(ax):
    order = [5, 3, 4, 0, 1, 2]  # back faces first
    for face in order:
        verts = cube_face_verts(face)
        poly = Poly3DCollection(
            [verts],
            facecolors=FACE_COLORS[face],
            edgecolors="#222222",
            linewidths=1.1,
            alpha=0.88,
        )
        ax.add_collection3d(poly)

    # Face labels on visible faces
    ax.text(1.35, 0, 0, "0\n+X", ha="center", va="center", fontsize=9, fontweight="bold", color=FACE_COLORS[0])
    ax.text(0, 1.35, 0, "1\n+Y", ha="center", va="center", fontsize=9, fontweight="bold", color=FACE_COLORS[1])
    ax.text(0, 0, 1.28, "2  +Z", ha="center", va="center", fontsize=9, fontweight="bold", color=FACE_COLORS[2])
    ax.text(-1.35, 0, 0, "3\n−X", ha="center", va="center", fontsize=9, fontweight="bold", color=FACE_COLORS[3])

    # Cube corner (+1,+1,+1): three faces meet → Z wins → face 2
    ax.scatter([1], [1], [1], s=70, c="#111111", depthshade=False, zorder=10)
    ax.text(1.12, 1.12, 1.18, "cube corner\n→ face 2 (+Z)", fontsize=7.5, color="#111111")

    # Edge between +X and +Y
    ax.plot([1, 1], [1, 1], [-1, 1], color="#111111", lw=2.2)
    ax.text(1.15, 1.15, -0.15, "edge |x|=|y|\n→ Y wins → face 1", fontsize=7.5, color="#111111")

    ax.set_xlim(-1.7, 1.7)
    ax.set_ylim(-1.7, 1.7)
    ax.set_zlim(-1.7, 1.7)
    ax.set_axis_off()
    ax.view_init(elev=18, azim=35)
    ax.set_title("1. Pick one cube face\nlargest_abs_component: later axis wins", fontsize=11, pad=6)


def draw_half_open_grid(ax):
    # 2×2 cells in (i, j). Vertex at the origin. Winner is NE: [i0, i0+size) × [j0, j0+size)
    cells = [
        # (x0, y0, label, color, owner)
        (-1, -1, "SW\nexcluded", "#cfcfcf", False),
        (0, -1, "SE\nexcluded", "#cfcfcf", False),
        (-1, 0, "NW\nexcluded", "#cfcfcf", False),
        (0, 0, "NE  WINNER\nthis cell’s\nmin-i, min-j\ncorner", "#54A24B", True),
    ]
    for x0, y0, label, color, owner in cells:
        rect = Rectangle(
            (x0, y0),
            1,
            1,
            facecolor=color,
            edgecolor="#222222",
            linewidth=1.6 if owner else 1.0,
            alpha=0.92,
            zorder=1,
        )
        ax.add_patch(rect)
        ax.text(
            x0 + 0.5,
            y0 + 0.5,
            label,
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold" if owner else "normal",
            color="#0d3b0d" if owner else "#555555",
            zorder=2,
        )

    # Half-open edges of the winner: included (solid green) vs excluded (dashed)
    ax.plot([0, 1], [0, 0], color="#1b6b1b", lw=3.2, solid_capstyle="butt", zorder=3)  # min-j included
    ax.plot([0, 0], [0, 1], color="#1b6b1b", lw=3.2, solid_capstyle="butt", zorder=3)  # min-i included
    ax.plot([1, 1], [0, 1], color="#888888", lw=2.2, ls=(0, (4, 3)), zorder=3)  # max-i excluded
    ax.plot([0, 1], [1, 1], color="#888888", lw=2.2, ls=(0, (4, 3)), zorder=3)  # max-j excluded

    ax.plot(0, 0, "o", ms=14, color="#111111", zorder=5)
    ax.annotate(
        "shared vertex\nfloor(s), floor(t)\n→ this (i, j)",
        xy=(0, 0),
        xytext=(-0.92, -0.55),
        fontsize=8,
        color="#111111",
        arrowprops=dict(arrowstyle="->", color="#111111", lw=1.1),
        ha="left",
        va="center",
    )

    ax.annotate("[ included", xy=(0.5, 0.0), xytext=(0.5, -0.22), ha="center", fontsize=8, color="#1b6b1b")
    ax.annotate("[ included", xy=(0.0, 0.5), xytext=(-0.22, 0.5), ha="center", va="center", rotation=90, fontsize=8, color="#1b6b1b")
    ax.annotate(") excluded", xy=(1.0, 0.5), xytext=(1.18, 0.5), ha="center", va="center", rotation=90, fontsize=8, color="#666666")
    ax.annotate(") excluded", xy=(0.5, 1.0), xytext=(0.5, 1.16), ha="center", fontsize=8, color="#666666")

    ax.annotate("", xy=(1.45, -1.05), xytext=(-1.05, -1.05), arrowprops=dict(arrowstyle="->", lw=1.2, color="#333"))
    ax.annotate("", xy=(-1.05, 1.45), xytext=(-1.05, -1.05), arrowprops=dict(arrowstyle="->", lw=1.2, color="#333"))
    ax.text(1.48, -1.05, "i  (s)", fontsize=9, va="center")
    ax.text(-1.05, 1.50, "j  (t)", fontsize=9, ha="center")

    ax.set_xlim(-1.55, 1.65)
    ax.set_ylim(-1.45, 1.55)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title("2. On that face: half-open floor\ncell owns [i, i+size) × [j, j+size)", fontsize=11, pad=6)


def draw_cell_191(ax):
    token = "191"
    cid = s2.CellId.from_token(token)
    level = cid.level()
    cell = s2.Cell(cid)
    xyz = [cell.get_vertex(i) for i in range(4)]

    palette = ["#4C78A8", "#F58518", "#54A24B", "#E45756"]
    assigned_order = []
    rows = []
    for i, p in enumerate(xyz, 1):
        ll = s2.LatLng.from_point(p)
        lat, lon = ll.lat().degrees, ll.lng().degrees
        assigned = latlon2s2(lat, lon, level)
        rows.append((i, lat, lon, assigned))
        if assigned not in assigned_order:
            assigned_order.append(assigned)
    color = {t: palette[i] for i, t in enumerate(assigned_order)}

    # Fill cells
    input_poly = s22geo(token)
    xs, ys = input_poly.exterior.xy
    ax.fill(xs, ys, facecolor="#eeeeee", edgecolor="#222222", lw=1.6, zorder=1)

    labeled = set()
    for t in assigned_order:
        poly = s22geo(t)
        c = color[t]
        xs, ys = poly.exterior.xy
        ax.fill(xs, ys, facecolor=c, edgecolor=c, lw=1.2, alpha=0.22, zorder=2)
        if t in labeled:
            continue
        labeled.add(t)
        ax.annotate(
            t,
            (poly.centroid.x, poly.centroid.y),
            ha="center",
            va="center",
            fontsize=9,
            color=c,
            fontweight="bold",
            zorder=4,
        )

    corners = ["SW  min-i, min-j\nOWNED by 191", "SE", "NE", "NW"]
    for (i, lat, lon, assigned), corner in zip(rows, corners):
        c = color[assigned]
        ax.plot(lon, lat, "o", color=c, ms=10, zorder=6)
        extra = f"  → {assigned}"
        ax.annotate(
            f"v{i}  {corner}{extra if i > 1 else extra}",
            (lon, lat),
            textcoords="offset points",
            xytext=(8, 8) if i != 1 else (8, -22),
            fontsize=8,
            color=c,
            fontweight="bold",
            zorder=7,
        )

    ax.set_aspect("equal")
    ax.set_xlabel("lon")
    ax.set_ylabel("lat")
    ax.set_title(
        "3. Cell 191 (level 4): v1 is the min-i, min-j corner, so it stays on 191",
        fontsize=11,
        pad=6,
    )


def main():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titlesize": 11,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )
    fig = plt.figure(figsize=(13.2, 8.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.15], hspace=0.32, wspace=0.18)

    ax_cube = fig.add_subplot(gs[0, 0], projection="3d")
    draw_cube(ax_cube)

    ax_grid = fig.add_subplot(gs[0, 1])
    draw_half_open_grid(ax_grid)

    ax_map = fig.add_subplot(gs[1, :])
    draw_cell_191(ax_map)

    fig.suptitle(
        "How S2 chooses a cell for a shared vertex\n"
        "xyz_to_face_uv  →  uv_to_st  →  floor(MAX_SIZE · s), floor(MAX_SIZE · t)  →  parent(level)",
        fontsize=13,
        fontweight="bold",
        y=0.98,
    )
    fig.savefig(OUT, dpi=160, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
