import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull


def plot_convex_hull(points):
    hull = ConvexHull(points)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    # Plot vertices
    ax.plot(points[:, 0], points[:, 1], points[:, 2], "o")

    # Plot faces
    for simplex in hull.simplices:
        s = hull.points[simplex]
        s = np.vstack((s, s[0]))
        ax.plot(s[:, 0], s[:, 1], s[:, 2], "k-")
    faces = Poly3DCollection(
        [hull.points[simplex] for simplex in hull.simplices], alpha=0.5
    )
    faces.set_facecolor("c")
    ax.add_collection3d(faces)
    plt.show()


def main():
    # Generate sample 3D points for testing
    # Create points on a sphere surface
    np.random.seed(42)
    n_points = 20
    theta = np.random.uniform(0, 2 * np.pi, n_points)
    phi = np.random.uniform(0, np.pi, n_points)
    r = 1.0

    x = r * np.sin(phi) * np.cos(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(phi)

    points = np.column_stack((x, y, z))

    print(f"Testing plot_convex_hull with {n_points} points")
    print(f"Points shape: {points.shape}")
    print(
        f"Points range: x=[{x.min():.2f}, {x.max():.2f}], "
        f"y=[{y.min():.2f}, {y.max():.2f}], "
        f"z=[{z.min():.2f}, {z.max():.2f}]"
    )

    plot_convex_hull(points)


if __name__ == "__main__":
    main()
