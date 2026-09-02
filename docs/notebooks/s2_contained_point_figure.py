"""Illustration of how S2 assigns an interior (contained) point."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

from vgrid.conversion.latlon2dggs import latlon2s2
from vgrid.conversion.dggs2geo.s22geo import s22geo
from vgrid.dggs import s2

OUT = Path(__file__).with_name("s2_contained_point.png")


def draw_interior_st(ax):
    ax.add_patch(Rectangle((0, 0), 1, 1, facecolor="#e8f5e6", edgecolor="none", zorder=0))
    ax.add_patch(Rectangle((1, 0), 1, 1, facecolor="#f4f4f4", edgecolor="#bbbbbb", lw=1, zorder=0))
    ax.add_patch(Rectangle((0, 1), 1, 1, facecolor="#f4f4f4", edgecolor="#bbbbbb", lw=1, zorder=0))
    ax.add_patch(Rectangle((1, 1), 1, 1, facecolor="#f4f4f4", edgecolor="#bbbbbb", lw=1, zorder=0))
    ax.add_patch(Rectangle((-1, 0), 1, 1, facecolor="#f4f4f4", edgecolor="#bbbbbb", lw=1, zorder=0))
    ax.add_patch(Rectangle((0, -1), 1, 1, facecolor="#f4f4f4", edgecolor="#bbbbbb", lw=1, zorder=0))

    ax.text(0.5, 0.42, "this cell\n[i, i+size) × [j, j+size)", ha="center", va="center", fontsize=9, color="#1b4d1b", fontweight="bold")

    ax.plot([0, 0], [0, 1], color="#1b6b1b", lw=3.5, zorder=3)
    ax.plot([0, 1], [0, 0], color="#1b6b1b", lw=3.5, zorder=3)
    ax.plot([1, 1], [0, 1], color="#c44e52", lw=2.4, ls=(0, (5, 3)), zorder=3)
    ax.plot([0, 1], [1, 1], color="#c44e52", lw=2.4, ls=(0, (5, 3)), zorder=3)

    # Interior samples
    pts = [(0.5, 0.5, "center"), (0.28, 0.62, ""), (0.72, 0.3, ""), (0.22, 0.22, "")]
    for x, y, label in pts:
        ax.plot(x, y, "*", ms=14, color="#4C78A8", markeredgecolor="#1a3a5c", zorder=5)
        if label:
            ax.annotate(label, (x, y), xytext=(8, 6), textcoords="offset points", fontsize=8, color="#4C78A8", fontweight="bold")

    ax.annotate(
        "strictly inside\n→ unique cell\n(no tie)",
        xy=(0.5, 0.5),
        xytext=(1.12, -0.55),
        fontsize=8,
        color="#4C78A8",
        fontweight="bold",
        arrowprops=dict(arrowstyle="->", color="#4C78A8", lw=1.1),
    )

    ax.set_xlim(-1.2, 2.2)
    ax.set_ylim(-1.25, 2.15)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.annotate("", xy=(2.05, -1.05), xytext=(-1.05, -1.05), arrowprops=dict(arrowstyle="->", lw=1.1, color="#333"))
    ax.annotate("", xy=(-1.05, 2.05), xytext=(-1.05, -1.05), arrowprops=dict(arrowstyle="->", lw=1.1, color="#333"))
    ax.text(2.08, -1.05, "i  (s)", fontsize=9, va="center")
    ax.text(-1.05, 2.12, "j  (t)", fontsize=9, ha="center")
    ax.set_title(
        "Interior: (s, t) not on a grid line.\nfloor lands in this cell’s half-open square — no neighbor can match.",
        fontsize=10,
        pad=8,
    )


def draw_cell_191_interior(ax):
    token = "191"
    cid = s2.CellId.from_token(token)
    level = cid.level()
    cell = s2.Cell(cid)
    center = cell.get_center()

    neighbors = [s2.CellId.to_token(n) for n in cid.get_edge_neighbors()]
    palette = {token: "#4C78A8"}
    nb_colors = ["#F58518", "#E45756", "#54A24B", "#B279A2"]
    for t, c in zip(neighbors, nb_colors):
        palette[t] = c

    for t, c in palette.items():
        poly = s22geo(t)
        xs, ys = poly.exterior.xy
        ax.fill(xs, ys, facecolor=c, edgecolor=c, lw=1.2, alpha=0.18 if t != token else 0.35, zorder=2)
        ax.annotate(t, (poly.centroid.x, poly.centroid.y), ha="center", va="center", fontsize=9, color=c, fontweight="bold", zorder=4)

    samples = [("center", center)]
    for k in range(4):
        v = cell.get_vertex(k)
        p = s2.Point(
            0.35 * v[0] + 0.65 * center[0],
            0.35 * v[1] + 0.65 * center[1],
            0.35 * v[2] + 0.65 * center[2],
        ).normalize()
        samples.append((f"in{k + 1}", p))

    all_191 = True
    for name, p in samples:
        ll = s2.LatLng.from_point(p)
        lat, lon = ll.lat().degrees, ll.lng().degrees
        assigned = latlon2s2(lat, lon, level)
        contained = cell.contains(p)
        if assigned != token:
            all_191 = False
        ax.plot(lon, lat, "*", ms=13, color="#4C78A8", markeredgecolor="#1a3a5c", zorder=6)
        if name == "center":
            ax.annotate(
                f"center → {assigned}\nCell.contains = {contained}",
                (lon, lat),
                textcoords="offset points",
                xytext=(10, -18),
                fontsize=8,
                color="#4C78A8",
                fontweight="bold",
                zorder=7,
            )

    ax.set_aspect("equal")
    ax.set_xlabel("lon")
    ax.set_ylabel("lat")
    ax.set_title(
        f"Cell 191: center and interior samples all → 191  (contains={all_191})",
        fontsize=11,
        pad=6,
    )
    ax.legend(
        handles=[Line2D([0], [0], marker="*", color="w", markerfacecolor="#4C78A8", markeredgecolor="#1a3a5c", markersize=12, label="interior point")],
        loc="lower right",
        fontsize=8,
    )


def main():
    plt.rcParams.update({"font.family": "DejaVu Sans", "figure.facecolor": "white", "savefig.facecolor": "white"})
    fig = plt.figure(figsize=(13.2, 6.4))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.15], wspace=0.22)
    draw_interior_st(fig.add_subplot(gs[0, 0]))
    draw_cell_191_interior(fig.add_subplot(gs[0, 1]))
    fig.suptitle(
        "How S2 assigns a contained (interior) point\n"
        "Same pipeline: xyz_to_face_uv → uv_to_st → floor → parent(level). No extra rule.",
        fontsize=13,
        fontweight="bold",
        y=1.02,
    )
    fig.savefig(OUT, dpi=160, bbox_inches="tight")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
