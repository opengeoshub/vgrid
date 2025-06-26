from shapely.geometry import shape, Polygon
import argparse
import json
from vgrid.utils import qtm
from vgrid.generator.settings import geodesic_dggs_to_feature
from shapely.ops import unary_union
from tqdm import tqdm

p90_n180, p90_n90, p90_p0, p90_p90, p90_p180 = (
    (90.0, -180.0),
    (90.0, -90.0),
    (90.0, 0.0),
    (90.0, 90.0),
    (90.0, 180.0),
)
p0_n180, p0_n90, p0_p0, p0_p90, p0_p180 = (
    (0.0, -180.0),
    (0.0, -90.0),
    (0.0, 0.0),
    (0.0, 90.0),
    (0.0, 180.0),
)
n90_n180, n90_n90, n90_p0, n90_p90, n90_p180 = (
    (-90.0, -180.0),
    (-90.0, -90.0),
    (-90.0, 0.0),
    (-90.0, 90.0),
    (-90.0, 180.0),
)


def generate_grid(resolution):
    levelFacets = {}
    QTMID = {}

    for lvl in tqdm(range(resolution), desc="Generating QTM DGGS"):
        levelFacets[lvl] = []
        QTMID[lvl] = []
        qtm_features = []  # Store GeoJSON features separately

        if lvl == 0:
            initial_facets = [
                [p0_n180, p0_n90, p90_n90, p90_n180, p0_n180, True],
                [p0_n90, p0_p0, p90_p0, p90_n90, p0_n90, True],
                [p0_p0, p0_p90, p90_p90, p90_p0, p0_p0, True],
                [p0_p90, p0_p180, p90_p180, p90_p90, p0_p90, True],
                [n90_n180, n90_n90, p0_n90, p0_n180, n90_n180, False],
                [n90_n90, n90_p0, p0_p0, p0_n90, n90_n90, False],
                [n90_p0, n90_p90, p0_p90, p0_p0, n90_p0, False],
                [n90_p90, n90_p180, p0_p180, p0_p90, n90_p90, False],
            ]

            for i, facet in enumerate(initial_facets):
                facet_geom = qtm.constructGeometry(facet)
                QTMID[0].append(str(i + 1))
                levelFacets[0].append(facet)
                qtm_id = QTMID[0][i]
                num_edges = 3
                qtm_feature = geodesic_dggs_to_feature(
                    "qtm", qtm_id, resolution, facet_geom, num_edges
                )
                qtm_features.append(qtm_feature)

        else:
            for i, pf in enumerate(levelFacets[lvl - 1]):
                subdivided_facets = qtm.divideFacet(pf)
                for j, subfacet in enumerate(subdivided_facets):
                    new_id = QTMID[lvl - 1][i] + str(j)
                    QTMID[lvl].append(new_id)
                    levelFacets[lvl].append(subfacet)
                    if lvl == resolution - 1:
                        subfacet_geom = qtm.constructGeometry(subfacet)
                        qtm_id = new_id
                        num_edges = 3
                        qtm_feature = geodesic_dggs_to_feature(
                            "qtm", qtm_id, resolution, subfacet_geom, num_edges
                        )
                        qtm_features.append(qtm_feature)
    return {"type": "FeatureCollection", "features": qtm_features}


def generate_grid_within_bbox(resolution, bbox):
    """Generates a Dutton QTM grid at a specific resolution within a bounding box and saves it as GeoJSON."""
    levelFacets = {}
    QTMID = {}
    qtm_features = []

    # Convert bbox to Polygon
    bbox_poly = Polygon(
        [
            (bbox[0], bbox[1]),  # min_lon, min_lat
            (bbox[2], bbox[1]),  # max_lon, min_lat
            (bbox[2], bbox[3]),  # max_lon, max_lat
            (bbox[0], bbox[3]),  # min_lon, max_lat
            (bbox[0], bbox[1]),  # Close the polygon
        ]
    )

    for lvl in tqdm(range(resolution), desc="Generating QTM DGGS"):
        levelFacets[lvl] = []
        QTMID[lvl] = []

        if lvl == 0:
            initial_facets = [
                [p0_n180, p0_n90, p90_n90, p90_n180, p0_n180, True],
                [p0_n90, p0_p0, p90_p0, p90_n90, p0_n90, True],
                [p0_p0, p0_p90, p90_p90, p90_p0, p0_p0, True],
                [p0_p90, p0_p180, p90_p180, p90_p90, p0_p90, True],
                [n90_n180, n90_n90, p0_n90, p0_n180, n90_n180, False],
                [n90_n90, n90_p0, p0_p0, p0_n90, n90_n90, False],
                [n90_p0, n90_p90, p0_p90, p0_p0, n90_p0, False],
                [n90_p90, n90_p180, p0_p180, p0_p90, n90_p90, False],
            ]

            for i, facet in enumerate(initial_facets):
                QTMID[0].append(str(i + 1))
                facet_geom = qtm.constructGeometry(facet)
                levelFacets[0].append(facet)
                if shape(facet_geom).intersects(bbox_poly) and resolution == 1:
                    qtm_id = QTMID[0][i]
                    num_edges = 3
                    qtm_feature = geodesic_dggs_to_feature(
                        "qtm", qtm_id, resolution, facet_geom, num_edges
                    )
                    qtm_features.append(qtm_feature)
        else:
            for i, pf in enumerate(levelFacets[lvl - 1]):
                subdivided_facets = qtm.divideFacet(pf)
                for j, subfacet in enumerate(subdivided_facets):
                    subfacet_geom = qtm.constructGeometry(subfacet)
                    if shape(subfacet_geom).intersects(
                        bbox_poly
                    ):  # Only keep intersecting facets
                        new_id = QTMID[lvl - 1][i] + str(j)
                        QTMID[lvl].append(new_id)
                        levelFacets[lvl].append(subfacet)
                        if (
                            lvl == resolution - 1
                        ):  # Only store final resolution in GeoJSON
                            qtm_id = new_id
                            num_edges = 3
                            qtm_feature = geodesic_dggs_to_feature(
                                "qtm", qtm_id, resolution, subfacet_geom, num_edges
                            )
                            qtm_features.append(qtm_feature)
    return {"type": "FeatureCollection", "features": qtm_features}


def generate_grid_resample(resolution, geojson_features):
    """Generates a Dutton QTM grid at a specific resolution within geojson_features and returns it as GeoJSON."""
    levelFacets = {}
    QTMID = {}
    qtm_features = []

    # Step 1: Union all input GeoJSON geometries
    geometries = [
        shape(feature["geometry"]) for feature in geojson_features["features"]
    ]
    unified_geom = unary_union(geometries)

    for lvl in tqdm(range(resolution), desc="Generating QTM DGGS"):
        levelFacets[lvl] = []
        QTMID[lvl] = []

        if lvl == 0:
            initial_facets = [
                [p0_n180, p0_n90, p90_n90, p90_n180, p0_n180, True],
                [p0_n90, p0_p0, p90_p0, p90_n90, p0_n90, True],
                [p0_p0, p0_p90, p90_p90, p90_p0, p0_p0, True],
                [p0_p90, p0_p180, p90_p180, p90_p90, p0_p90, True],
                [n90_n180, n90_n90, p0_n90, p0_n180, n90_n180, False],
                [n90_n90, n90_p0, p0_p0, p0_n90, n90_n90, False],
                [n90_p0, n90_p90, p0_p90, p0_p0, n90_p0, False],
                [n90_p90, n90_p180, p0_p180, p0_p90, n90_p90, False],
            ]

            for i, facet in enumerate(initial_facets):
                QTMID[0].append(str(i + 1))
                facet_geom = qtm.constructGeometry(facet)
                levelFacets[0].append(facet)

                if shape(facet_geom).intersects(unified_geom) and resolution == 1:
                    qtm_id = QTMID[0][i]
                    num_edges = 3
                    qtm_feature = geodesic_dggs_to_feature(
                        "qtm", qtm_id, resolution, facet_geom, num_edges
                    )
                    qtm_features.append(qtm_feature)

        else:
            for i, pf in enumerate(levelFacets[lvl - 1]):
                subdivided_facets = qtm.divideFacet(pf)
                for j, subfacet in enumerate(subdivided_facets):
                    subfacet_geom = qtm.constructGeometry(subfacet)
                    if shape(subfacet_geom).intersects(unified_geom):
                        new_id = QTMID[lvl - 1][i] + str(j)
                        QTMID[lvl].append(new_id)
                        levelFacets[lvl].append(subfacet)

                        if lvl == resolution - 1:
                            qtm_id = new_id
                            num_edges = 3
                            qtm_feature = geodesic_dggs_to_feature(
                                "qtm", qtm_id, resolution, subfacet_geom, num_edges
                            )
                            qtm_features.append(qtm_feature)

    return {"type": "FeatureCollection", "features": qtm_features}


def main():
    parser = argparse.ArgumentParser(description="Generate QTM DGGS.")
    parser.add_argument(
        "-r", "--resolution", required=True, type=int, help="Resolution [1..24]."
    )
    parser.add_argument(
        "-b",
        "--bbox",
        type=float,
        nargs=4,
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)",
    )

    args = parser.parse_args()
    resolution = args.resolution
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]

    if resolution < 1 or resolution > 24:
        print("Please select a resolution in [1..24] range and try again ")
        return

    if bbox == [-180, -90, 180, 90]:
        geojson_features = generate_grid(resolution)

    else:
        geojson_features = generate_grid_within_bbox(resolution, bbox)

    if geojson_features:
        geojson_path = f"qtm_grid_{resolution}.geojson"
        with open(geojson_path, "w") as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()
