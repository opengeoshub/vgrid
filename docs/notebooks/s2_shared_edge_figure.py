"""Illustration of how S2 assigns a point on a shared edge (from vgrid/dggs/s2.py)."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from vgrid.conversion.latlon2dggs import latlon2s2
from vgrid.conversion.dggs2geo.s22geo import s22geo
from vgrid.dggs import s2

OUT = Path(__file__).with_name("s2_shared_edge.png")

OWN = "#54A24B"
LOSE = "#cfcfcf"
OWN_EDGE = "#1b6b1b"
LOSE_EDGE = "#c44e52"


def draw_half_open_edges(ax):
    """One cell: min-i / min-j edges included, max-i / max-j excluded."""
    ax.add_patch(
        Rectangle((0, 0), 1, 1, facecolor="#e8f5e6", edgecolor="none", zorder=0)
    )
    # Neighbors (faded)
    ax.add_patch(Rectangle((1, 0), 1, 1, facecolor="#f4f4f4", edgecolor="#bbbbbb", lw=1, zorder=0))
    ax.add_patch(Rectangle((0, 1), 1, 1, facecolor="#f4f4f4", edgecolor="#bbbbbb", lw=1, zorder=0))
    ax.add_patch(Rectangle((1, 1), 1, 1, facecolor="#f4f4f4", edgecolor="#bbbbbb", lw=1, zorder=0))
    ax.add_patch(Rectangle((-1, 0), 1, 1, facecolor="#f4f4f4", edgecolor="#bbbbbb", lw=1, zorder=0))
    ax.add_patch(Rectangle((0, -1), 1, 1, facecolor="#f4f4f4", edgecolor="#bbbbbb", lw=1, zorder=0))

    ax.text(0.5, 0.5, "this cell\n[i, i+size)\n[j, j+size)", ha="center", va="center", fontsize=9, color="#1b4d1b", fontweight="bold")
    ax.text(1.5, 0.5, "+i neighbor", ha="center", va="center", fontsize=8, color="#666666")
    ax.text(0.5, 1.5, "+j neighbor", ha="center", va="center", fontsize=8, color="#666666")

    # Included edges (solid green): min-i (left) and min-j (bottom)
    ax.plot([0, 0], [0, 1], color=OWN_EDGE, lw=4.0, solid_capstyle="butt", zorder=3)
    ax.plot([0, 1], [0, 0], color=OWN_EDGE, lw=4.0, solid_capstyle="butt", zorder=3)
    # Excluded edges (dashed red): max-i (right) and max-j (top)
    ax.plot([1, 1], [0, 1], color=LOSE_EDGE, lw=3.2, ls=(0, (5, 3)), zorder=3)
    ax.plot([0, 1], [1, 1], color=LOSE_EDGE, lw=3.2, ls=(0, (5, 3)), zorder=3)

    # Edge midpoints (squares) — interior of the edge, not a vertex
    mids = [
        (0.5, 0.0, OWN_EDGE, "min-j  (south)\nfloor → this cell"),
        (0.0, 0.5, OWN_EDGE, "min-i  (west)\nfloor → this cell"),
        (1.0, 0.5, LOSE_EDGE, "max-i  (east)\nfloor → +i neighbor"),
        (0.5, 1.0, LOSE_EDGE, "max-j  (north)\nfloor → +j neighbor"),
    ]
    offsets = [(0.08, -0.38), (-0.95, 0.12), (0.12, 0.12), (0.08, 0.28)]
    for (x, y, c, text), (dx, dy) in zip(mids, offsets):
        ax.plot(x, y, "s", ms=11, color=c, markeredgecolor="#111111", markeredgewidth=0.6, zorder=5)
        ax.annotate(text, (x, y), xytext=(x + dx, y + dy), fontsize=7.5, color=c, fontweight="bold")

    # Vertices as small circles, for contrast
    for x, y in [(0, 0), (1, 0), (1, 1), (0, 1)]:
        ax.plot(x, y, "o", ms=6, color="#888888", zorder=4)

    ax.annotate("", xy=(2.15, -1.15), xytext=(-1.15, -1.15), arrowprops=dict(arrowstyle="->", lw=1.1, color="#333"))
    ax.annotate("", xy=(-1.15, 2.15), xytext=(-1.15, -1.15), arrowprops=dict(arrowstyle="->", lw=1.1, color="#333"))
    ax.text(2.18, -1.15, "i  (s)", fontsize=9, va="center")
    ax.text(-1.15, 2.22, "j  (t)", fontsize=9, ha="center")

    ax.set_xlim(-1.35, 2.35)
    ax.set_ylim(-1.55, 2.35)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.set_title(
        "A point in the interior of an edge is shared by 2 cells.\n"
        "floor keeps the min-i and min-j sides; the opposite sides go to the +i / +j neighbor.",
        fontsize=10,
        pad=8,
    )


def draw_cell_191_edges(ax):
    token = "191"
    cid = s2.CellId.from_token(token)
    level = cid.level()
    cell = s2.Cell(cid)
    xyz = [cell.get_vertex(i) for i in range(4)]

    palette = ["#4C78A8", "#F58518", "#E45756", "#54A24B"]
    assigned_order = []
    verts = []
    mids = []
    for i, p in enumerate(xyz, 1):
        ll = s2.LatLng.from_point(p)
        verts.append((i, ll.lat().degrees, ll.lng().degrees, latlon2s2(ll.lat().degrees, ll.lng().degrees, level)))
    for i in range(4):
        mid = (xyz[i] + xyz[(i + 1) % 4]).normalize()
        ll = s2.LatLng.from_point(mid)
        lat, lon = ll.lat().degrees, ll.lng().degrees
        assigned = latlon2s2(lat, lon, level)
        name = f"v{i + 1}_v{(i + 1) % 4 + 1}"
        mids.append((name, lat, lon, assigned))
        if assigned not in assigned_order:
            assigned_order.append(assigned)
    for _, _, _, a in verts:
        if a not in assigned_order:
            assigned_order.append(a)
    color = {t: palette[i % len(palette)] for i, t in enumerate(assigned_order)}

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

    edge_notes = {
        "v1_v2": "south  min-j  → 191",
        "v2_v3": "east  max-i  → 18f",
        "v3_v4": "north  max-j  → 19b",
        "v4_v1": "west  min-i  → 191",
    }
    mid_offsets = {
        "v1_v2": (8, -18),
        "v2_v3": (8, 8),
        "v3_v4": (8, 8),
        "v4_v1": (-88, 8),
    }
    for name, lat, lon, assigned in mids:
        c = color[assigned]
        ax.plot(lon, lat, "s", color=c, ms=10, markeredgecolor="#111111", markeredgewidth=0.4, zorder=6)
        ax.annotate(
            f"{name}\n{edge_notes[name]}",
            (lon, lat),
            textcoords="offset points",
            xytext=mid_offsets[name],
            fontsize=8,
            color=c,
            fontweight="bold",
            zorder=7,
        )

    for i, lat, lon, assigned in verts:
        ax.plot(lon, lat, "o", color="#888888", ms=5, zorder=5)

    handles = [
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#4C78A8", markeredgecolor="#111", markersize=8, label="edge midpoint"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#888888", markersize=6, label="vertex (see other figure)"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=8, framealpha=0.95)
    ax.set_aspect("equal")
    ax.set_xlabel("lon")
    ax.set_ylabel("lat")
    ax.set_title(
        "Cell 191: south and west midpoints stay on 191; east and north go to the neighbor",
        fontsize=11,
        pad=6,
    )


def main():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )
    fig = plt.figure(figsize=(13.2, 6.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.15], wspace=0.22)

    ax0 = fig.add_subplot(gs[0, 0])
    draw_half_open_edges(ax0)

    ax1 = fig.add_subplot(gs[0, 1])
    draw_cell_191_edges(ax1)

    fig.suptitle(
        "How S2 chooses a cell for a shared-edge point\n"
        "Same floor rule as a vertex: the point goes to the cell on the +i or +j side of the grid line",
        fontsize=13,
        fontweight="bold",
        y=1.02,
    )
    fig.savefig(OUT, dpi=160, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
