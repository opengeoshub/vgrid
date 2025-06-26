#!/usr/bin/env python3

import pandas as pd
import geopandas as gpd
from shapely.geometry import (
    Point,
    MultiPoint,
    Polygon,
    LineString,
    MultiLineString,
    MultiPolygon,
)
from vgrid.conversion.vector2dggs.vector2h3 import (
    vector2h3,
    point2h3,
    polyline2h3,
    polygon2h3,
)

from pyproj import Geod
import tempfile
import os

# Create geod object for testing
geod = Geod(ellps="WGS84")


def test_all_input_types():
    """Test vector2h3 with all supported input types."""
    print("Testing all input types...")

    # Test data
    test_data = {
        "point": Point(106.68, 10.77),
        "multipoint": MultiPoint([(106.68, 10.77), (106.69, 10.78), (106.67, 10.76)]),
        "linestring": LineString([(106.68, 10.77), (106.69, 10.78)]),
        "multilinestring": MultiLineString(
            [[(106.68, 10.77), (106.69, 10.78)], [(106.67, 10.76), (106.70, 10.79)]]
        ),
        "polygon": Polygon(
            [
                (106.67, 10.76),
                (106.69, 10.76),
                (106.69, 10.78),
                (106.67, 10.78),
                (106.67, 10.76),
            ]
        ),
        "multipolygon": MultiPolygon(
            [
                Polygon(
                    [
                        (106.67, 10.76),
                        (106.69, 10.76),
                        (106.69, 10.78),
                        (106.67, 10.78),
                        (106.67, 10.76),
                    ]
                ),
                Polygon(
                    [
                        (106.70, 10.79),
                        (106.72, 10.79),
                        (106.72, 10.81),
                        (106.70, 10.81),
                        (106.70, 10.79),
                    ]
                ),
            ]
        ),
    }

    # Test 1: Shapely Geometry Objects
    print("\n1. Testing Shapely Geometry Objects:")
    for geom_type, geom in test_data.items():
        try:
            result = vector2h3(geom, resolution=9)
            print(f"  ✓ {geom_type}: Generated {len(result['features'])} features")
        except Exception as e:
            print(f"  ✗ {geom_type}: {e}")

    # Test 2: List of Shapely Geometries
    print("\n2. Testing List of Shapely Geometries:")
    try:
        geom_list = [test_data["point"], test_data["linestring"], test_data["polygon"]]
        result = vector2h3(geom_list, resolution=9)
        print(f"  ✓ List of geometries: Generated {len(result['features'])} features")
    except Exception as e:
        print(f"  ✗ List of geometries: {e}")

    # Test 3: DataFrame with Geometry Column
    print("\n3. Testing DataFrame with Geometry Column:")
    try:
        df_data = {
            "id": [1, 2, 3],
            "name": ["point1", "line1", "polygon1"],
            "geometry": [
                test_data["point"],
                test_data["linestring"],
                test_data["polygon"],
            ],
        }
        df = pd.DataFrame(df_data)
        result = vector2h3(df, resolution=9)
        print(f"  ✓ DataFrame: Generated {len(result['features'])} features")
    except Exception as e:
        print(f"  ✗ DataFrame: {e}")

    # Test 4: GeoDataFrame
    print("\n4. Testing GeoDataFrame:")
    try:
        gdf_data = {
            "id": [1, 2, 3],
            "name": ["point1", "line1", "polygon1"],
            "geometry": [
                test_data["point"],
                test_data["linestring"],
                test_data["polygon"],
            ],
        }
        gdf = gpd.GeoDataFrame(gdf_data, crs="EPSG:4326")
        result = vector2h3(gdf, resolution=9)
        print(f"  ✓ GeoDataFrame: Generated {len(result['features'])} features")
    except Exception as e:
        print(f"  ✗ GeoDataFrame: {e}")

    # Test 5: GeoJSON Dictionary
    print("\n5. Testing GeoJSON Dictionary:")
    try:
        geojson_data = {
            "type": "FeatureCollection",
            "features": [
                {
                    "type": "Feature",
                    "geometry": {"type": "Point", "coordinates": [106.68, 10.77]},
                    "properties": {"name": "test_point", "type": "city"},
                },
                {
                    "type": "Feature",
                    "geometry": {
                        "type": "LineString",
                        "coordinates": [[106.68, 10.77], [106.69, 10.78]],
                    },
                    "properties": {"name": "test_line", "type": "road"},
                },
                {
                    "type": "Feature",
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [
                            [
                                [106.67, 10.76],
                                [106.69, 10.76],
                                [106.69, 10.78],
                                [106.67, 10.78],
                                [106.67, 10.76],
                            ]
                        ],
                    },
                    "properties": {"name": "test_polygon", "type": "area"},
                },
            ],
        }
        result = vector2h3(geojson_data, resolution=9)
        print(f"  ✓ GeoJSON: Generated {len(result['features'])} features")
        # Check that properties are preserved
        if result["features"]:
            print(f"    Properties preserved: {result['features'][0]['properties']}")
    except Exception as e:
        print(f"  ✗ GeoJSON: {e}")

    print("✓ All input types test completed!")


def test_all_output_formats():
    """Test vector2h3 with all supported output formats."""
    print("\nTesting all output formats...")

    # Test geometry
    test_geom = Polygon(
        [
            (106.67, 10.76),
            (106.69, 10.76),
            (106.69, 10.78),
            (106.67, 10.78),
            (106.67, 10.76),
        ]
    )

    # Create temporary directory for file outputs
    with tempfile.TemporaryDirectory() as temp_dir:
        # Test 1: GeoJSON (default)
        print("\n1. Testing GeoJSON output:")
        try:
            result = vector2h3(test_geom, resolution=9, output_format="geojson")
            print(f"  ✓ GeoJSON: Generated {len(result['features'])} features")
        except Exception as e:
            print(f"  ✗ GeoJSON: {e}")

        # Test 2: Shapefile
        print("\n2. Testing Shapefile output:")
        try:
            output_path = os.path.join(temp_dir, "test_output.shp")
            result = vector2h3(
                test_geom,
                resolution=9,
                output_format="shapefile",
                output_path=output_path,
            )
            print(f"  ✓ Shapefile: Saved to {result}")
        except Exception as e:
            print(f"  ✗ Shapefile: {e}")

        # Test 3: GeoPackage
        print("\n3. Testing GeoPackage output:")
        try:
            output_path = os.path.join(temp_dir, "test_output.gpkg")
            result = vector2h3(
                test_geom, resolution=9, output_format="gpkg", output_path=output_path
            )
            print(f"  ✓ GeoPackage: Saved to {result}")
        except Exception as e:
            print(f"  ✗ GeoPackage: {e}")

        # Test 4: CSV
        print("\n4. Testing CSV output:")
        try:
            output_path = os.path.join(temp_dir, "test_output.csv")
            result = vector2h3(
                test_geom, resolution=9, output_format="csv", output_path=output_path
            )
            print(f"  ✓ CSV: Saved to {result}")
        except Exception as e:
            print(f"  ✗ CSV: {e}")

        # Test 5: Parquet
        print("\n5. Testing Parquet output:")
        try:
            output_path = os.path.join(temp_dir, "test_output.parquet")
            result = vector2h3(
                test_geom,
                resolution=9,
                output_format="parquet",
                output_path=output_path,
            )
            print(f"  ✓ Parquet: Saved to {result}")
        except Exception as e:
            print(f"  ✗ Parquet: {e}")

        # Test 6: WKT
        print("\n6. Testing WKT output:")
        try:
            output_path = os.path.join(temp_dir, "test_output.wkt")
            result = vector2h3(
                test_geom, resolution=9, output_format="wkt", output_path=output_path
            )
            print(f"  ✓ WKT: Saved to {result}")
        except Exception as e:
            print(f"  ✗ WKT: {e}")

        # Test 7: WKB
        print("\n7. Testing WKB output:")
        try:
            output_path = os.path.join(temp_dir, "test_output.wkb")
            result = vector2h3(
                test_geom, resolution=9, output_format="wkb", output_path=output_path
            )
            print(f"  ✓ WKB: Saved to {result}")
        except Exception as e:
            print(f"  ✗ WKB: {e}")

    print("✓ All output formats test completed!")


def test_predicates():
    """Test vector2h3 with different spatial predicates."""
    print("\nTesting spatial predicates...")

    # Test polygon
    test_polygon = Polygon(
        [
            (106.67, 10.76),
            (106.69, 10.76),
            (106.69, 10.78),
            (106.67, 10.78),
            (106.67, 10.76),
        ]
    )

    predicates = {
        None: "intersects (default)",
        "intersects": "intersects",
        "within": "within",
        "centroid_within": "centroid_within",
        "largest_overlap": "intersection >= 50%",
    }

    for predicate, description in predicates.items():
        try:
            result = vector2h3(test_polygon, resolution=9, predicate=predicate)
            print(f"  ✓ {description}: Generated {len(result['features'])} features")
        except Exception as e:
            print(f"  ✗ {description}: {e}")

    print("✓ Spatial predicates test completed!")


def test_resolutions():
    """Test vector2h3 with different H3 resolutions."""
    print("\nTesting different H3 resolutions...")

    test_geom = Point(106.68, 10.77)

    for resolution in [0, 5, 9, 12, 15]:
        try:
            result = vector2h3(test_geom, resolution=resolution)
            print(
                f"  ✓ Resolution {resolution}: Generated {len(result['features'])} features"
            )
        except Exception as e:
            print(f"  ✗ Resolution {resolution}: {e}")

    print("✓ H3 resolutions test completed!")


def test_compact_mode():
    """Test vector2h3 with compact mode."""
    print("\nTesting compact mode...")

    test_polygon = Polygon(
        [
            (106.67, 10.76),
            (106.69, 10.76),
            (106.69, 10.78),
            (106.67, 10.78),
            (106.67, 10.76),
        ]
    )

    try:
        # Without compact mode
        result_normal = vector2h3(test_polygon, resolution=9, compact=False)
        print(f"  ✓ Normal mode: Generated {len(result_normal['features'])} features")

        # With compact mode
        result_compact = vector2h3(test_polygon, resolution=9, compact=True)
        print(f"  ✓ Compact mode: Generated {len(result_compact['features'])} features")

        if len(result_compact["features"]) <= len(result_normal["features"]):
            print("    ✓ Compact mode reduced features as expected")
        else:
            print("    ⚠ Compact mode did not reduce features")

    except Exception as e:
        print(f"  ✗ Compact mode: {e}")

    print("✓ Compact mode test completed!")


def test_individual_functions():
    """Test individual geometry functions directly."""
    print("\nTesting individual geometry functions...")

    # Test point2h3
    point = Point(106.68, 10.77)
    point_properties = {"name": "test_point", "type": "city"}
    try:
        result = point2h3(9, point, point_properties)
        print(f"  ✓ point2h3: Generated {len(result['features'])} features")
        print(f"    Properties: {result['features'][0]['properties']}")
    except Exception as e:
        print(f"  ✗ point2h3: {e}")

    # Test polyline2h3
    line = LineString([(106.68, 10.77), (106.69, 10.78)])
    line_properties = {"name": "test_line", "type": "road"}
    try:
        result = polyline2h3(9, line, line_properties)
        print(f"  ✓ polyline2h3: Generated {len(result['features'])} features")
        if result["features"]:
            print(f"    Properties: {result['features'][0]['properties']}")
    except Exception as e:
        print(f"  ✗ polyline2h3: {e}")

    # Test polygon2h3
    polygon = Polygon(
        [
            (106.67, 10.76),
            (106.69, 10.76),
            (106.69, 10.78),
            (106.67, 10.78),
            (106.67, 10.76),
        ]
    )
    polygon_properties = {"name": "test_polygon", "type": "area"}
    try:
        result = polygon2h3(9, polygon, polygon_properties)
        print(f"  ✓ polygon2h3: Generated {len(result['features'])} features")
        if result["features"]:
            print(f"    Properties: {result['features'][0]['properties']}")
    except Exception as e:
        print(f"  ✗ polygon2h3: {e}")

    print("✓ Individual functions test completed!")


def test_error_handling():
    """Test error handling for invalid inputs."""
    print("\nTesting error handling...")

    # Test invalid resolution
    try:
        vector2h3(Point(106.68, 10.77), resolution=20)
        print("  ✗ Should have failed for invalid resolution 20")
    except ValueError as e:
        print(f"  ✓ Correctly caught invalid resolution: {e}")
    except Exception as e:
        print(f"  ✗ Unexpected error for invalid resolution: {e}")

    # Test invalid input type
    try:
        vector2h3("invalid_input", resolution=9)
        print("  ✗ Should have failed for invalid input type")
    except ValueError as e:
        print(f"  ✓ Correctly caught invalid input type: {e}")
    except Exception as e:
        print(f"  ✗ Unexpected error for invalid input type: {e}")

    # Test invalid output format
    try:
        vector2h3(
            Point(106.68, 10.77), resolution=9, output_format="invalid_format"
        )
        print("  ✗ Should have failed for invalid output format")
    except ValueError as e:
        print(f"  ✓ Correctly caught invalid output format: {e}")
    except Exception as e:
        print(f"  ✗ Unexpected error for invalid output format: {e}")

    print("✓ Error handling test completed!")


if __name__ == "__main__":
    print("🧪 Comprehensive vector2h3 Testing Suite")
    print("=" * 50)

    test_all_input_types()
    test_all_output_formats()
    test_predicates()
    test_resolutions()
    test_compact_mode()
    test_individual_functions()
    test_error_handling()

    print("\n" + "=" * 50)
    print("🎉 All tests completed!")
