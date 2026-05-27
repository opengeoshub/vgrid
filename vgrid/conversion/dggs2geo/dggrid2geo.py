"""
DGGRID to Geometry Module

This module provides functionality to convert DGGRID cell IDs to Shapely Polygons and GeoJSON FeatureCollection.

Key Functions:
    dggrid2geo: Convert DGGRID cell IDs to Shapely Polygons
    dggrid2geojson: Convert DGGRID cell IDs to GeoJSON FeatureCollection
    dggrid2geo_cli: Command-line interface for polygon conversion
    dggrid2geojson_cli: Command-line interface for GeoJSON conversion
"""

import json
import argparse

from dggrid4py import dggs_types
from dggrid4py.dggrid_runner import output_address_types
from vgrid.utils.io import (
    validate_dggrid_type,
    validate_dggrid_resolution,
    create_dggrid_instance,
)


def dggrid2geo(
    dggrid_instance,
    dggs_type,
    dggrid_ids,
    resolution=None,
    input_address_type="SEQNUM",
    split_antimeridian=False,
    aggregate=False,
    options=None,
):
    """
    Convert DGGRID cell IDs to Shapely geometry objects.

    Accepts a single dggrid_id (string/int) or a list of dggrid_ids. For each valid DGGRID cell ID,
    creates a Shapely Polygon representing the grid cell boundaries. Skips invalid or
    error-prone cells.

    Parameters
    ----------
    dggrid_instance : object
        DGGRID instance for processing.
    dggs_type : str
        DGGRID DGGS type (e.g., "ISEA7H", "ISEA4T").
    dggrid_ids : str, int, or list of str/int
        DGGRID cell ID(s) to convert. Can be a single string/int or a list of strings/ints.
        Example format: "783229476878"
    resolution : int, optional
        Resolution level for the DGGS. If None, will be validated based on dggs_type.
    input_address_type : str, default "SEQNUM"
        Input address type for the cell IDs.
    split_antimeridian : bool, optional
        When True, apply antimeridian fixing to the resulting polygons.
        Defaults to False when None or omitted.
    options : dict, optional
        Options to pass to grid_cell_polygons_from_cellids. 
        For example: {"densification": 2} to add densification points.
        Defaults to None.

    Returns
    -------
    geopandas.GeoDataFrame
        A GeoDataFrame containing polygon geometries for each valid DGGRID cell.
        Each row includes:
        - geometry: Polygon representing the cell boundaries
        - dggrid_{dggs_type.lower()}: The original cell ID
        - resolution: The resolution level

    Examples
    --------
    >>> dggrid2geo(instance, "ISEA7H", "783229476878", 13)
    <geopandas.GeoDataFrame object at ...>

    >>> dggrid2geo(instance, "ISEA7H", ["783229476878", "783229476879"], 13)
    <geopandas.GeoDataFrame object at ...>
    """
    dggs_type = validate_dggrid_type(dggs_type)
    resolution = validate_dggrid_resolution(dggs_type, resolution)
    if isinstance(dggrid_ids, (str, int)):
        dggrid_ids = [dggrid_ids]

    # Convert from input_address_type to SEQNUM if needed
    if input_address_type and input_address_type != "SEQNUM":
        address_type_transform = dggrid_instance.address_transform(
            dggrid_ids,
            dggs_type=dggs_type,
            resolution=resolution,
            mixed_aperture_level=None,
            input_address_type=input_address_type,
            output_address_type="SEQNUM",
        )
        # Extract all SEQNUM values, not just the first one
        dggrid_ids = address_type_transform["SEQNUM"].tolist()

    kwargs = {
        "split_dateline": split_antimeridian,  # to prevent polygon divided into 2 polygons
    }
    if options:
        kwargs.update(options)
    dggrid_cells = dggrid_instance.grid_cell_polygons_from_cellids(
        dggrid_ids,
        dggs_type,
        resolution,
        **kwargs,
    )

    # Convert global_id back to input_address_type if needed
    if input_address_type and input_address_type != "SEQNUM":
        # Get the SEQNUM values from global_id column
        seqnum_values = dggrid_cells["global_id"].tolist()

        # Transform back to input_address_type
        reverse_transform = dggrid_instance.address_transform(
            seqnum_values,
            dggs_type=dggs_type,
            resolution=resolution,
            mixed_aperture_level=None,
            input_address_type="SEQNUM",
            output_address_type=input_address_type,
        )

        # Replace global_id values with the original input_address_type values
        dggrid_cells["global_id"] = reverse_transform[input_address_type].values

    # Rename global_id column to dggrid_{dggs_type.lower()}
    dggrid_cells = dggrid_cells.rename(
        columns={"global_id": f"dggrid_{dggs_type.lower()}"}
    )
    # Add resolution property
    dggrid_cells["resolution"] = resolution

    # Apply antimeridian fixing if requested
    if split_antimeridian:
        if aggregate:
            dggrid_cells = dggrid_cells.dissolve(by=f"dggrid_{dggs_type.lower()}", as_index=False)  
    return dggrid_cells


def dggrid2geo_cli():
    """
    Command-line interface for dggrid2geo supporting multiple DGGRID cell IDs.
    """
    parser = argparse.ArgumentParser(
        description="Convert DGGRID cell ID(s) to Shapely Polygons. \
                                     Usage: dggrid2geo <cell_ids> <dggs_type> <res> [input_address_type]. \
                                     Ex: dggrid2geo 783229476878 ISEA7H 13 SEQNUM"
    )
    parser.add_argument(
        "dggs_type",
        choices=dggs_types,
        help="Select a DGGS type from the available options.",
    )

    parser.add_argument("dggrid_ids", nargs="+", help="Input DGGRID cell ID(s)")

    parser.add_argument("resolution", type=int, help="resolution")
    parser.add_argument(
        "input_address_type",
        choices=output_address_types,
        default="SEQNUM",
        nargs="?",  # This makes the argument optional
        help="Select an input address type from the available options.",
    )

    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,
        help="Apply antimeridian fixing to the resulting polygons",
    )
    parser.add_argument(
        "-aggregate",
        "--aggregate",
        action="store_true",
        help="Aggregate the resulting polygons",
    )
    parser.add_argument(
        "-options",
        "--options",
        type=str,
        default=None,
        help="JSON string of options to pass to grid_cell_polygons_from_cellids. "
             "Example: '{\"densification\": 2}'",
    )
    args = parser.parse_args()
    dggrid_instance = create_dggrid_instance()
    
    # Parse options JSON if provided
    options = None
    if args.options:
        try:
            options = json.loads(args.options)
        except json.JSONDecodeError as e:
            print(f"Error: Invalid JSON in options: {str(e)}")
            return
    
    polys = dggrid2geo(
        dggrid_instance,
        args.dggs_type,
        args.dggrid_ids,
        args.resolution,
        args.input_address_type,
        split_antimeridian=args.split_antimeridian,
        aggregate=args.aggregate,
        options=options,
    )
    return polys


def dggrid2geojson(
    dggrid_instance,
    dggs_type,
    dggrid_ids,
    resolution,
    input_address_type="SEQNUM",
    split_antimeridian=False,
    aggregate=False,
    options=None,
):
    """
    Convert DGGRID cell IDs to GeoJSON FeatureCollection.

    Accepts a single dggrid_id (string/int) or a list of dggrid_ids. For each valid DGGRID cell ID,
    creates a GeoJSON feature with polygon geometry representing the grid cell boundaries.
    Skips invalid or error-prone cells.

    Parameters
    ----------
    dggrid_instance : object
        DGGRID instance for processing.
    dggs_type : str
        DGGRID DGGS type (e.g., "ISEA7H", "ISEA4T").
    dggrid_ids : str, int, or list of str/int
        DGGRID cell ID(s) to convert. Can be a single string/int or a list of strings/ints.
        Example format: "783229476878"
    resolution : int
        Resolution level for the DGGS.
    input_address_type : str, default "SEQNUM"
        Input address type for the cell IDs.
    split_antimeridian : bool, optional
        When True, apply antimeridian fixing to the resulting polygons.
        Defaults to False when None or omitted.
    options : dict, optional
        Options to pass to grid_cell_polygons_from_cellids. 
        For example: {"densification": 2} to add densification points.
        Defaults to None.

    Returns
    -------
    dict
        A GeoJSON FeatureCollection containing polygon features for each valid DGGRID cell.
        Each feature includes:
        - geometry: Polygon representing the cell boundaries
        - properties: Contains the DGGRID cell ID, resolution level, and cell metadata

    Examples
    --------
    >>> dggrid2geojson(instance, "ISEA7H", "783229476878", 13)
    {'type': 'FeatureCollection', 'features': [...]}

    >>> dggrid2geojson(instance, "ISEA7H", ["783229476878", "783229476879"], 13)
    {'type': 'FeatureCollection', 'features': [...]}
    """
    if isinstance(dggrid_ids, (str, int)):
        dggrid_ids = [dggrid_ids]

    # Get the GeoDataFrame from dggrid2geo
    gdf = dggrid2geo(
        dggrid_instance,
        dggs_type,
        dggrid_ids,
        resolution,
        input_address_type,
        split_antimeridian=split_antimeridian,
        aggregate=aggregate,
        options=options,
    )
    # Normalize CRS for GeoJSON export to avoid non-URN CRS warnings.
    if gdf.crs is not None and gdf.crs != "EPSG:4326":
        try:
            gdf = gdf.to_crs("EPSG:4326")
        except Exception:
            # DGGRID can emit custom authalic CRS metadata; keep coordinates and
            # mark as EPSG:4326 for interoperable GeoJSON output.
            gdf = gdf.set_crs("EPSG:4326", allow_override=True)
    # Convert GeoDataFrame to GeoJSON dictionary
    geojson_dict = json.loads(gdf.to_json())

    return geojson_dict


def dggrid2geojson_cli():
    """
    Command-line interface for dggrid2geojson supporting multiple DGGRID cell IDs.
    """
    parser = argparse.ArgumentParser(
        description="Convert DGGRID cell ID(s) to GeoJSON. \
                                     Usage: dggrid2geojson <cell_ids> <dggs_type> <res> [input_address_type]. \
                                     Ex: dggrid2geojson 783229476878 ISEA7H 13 SEQNUM"
    )

    parser.add_argument(
        "dggs_type",
        choices=dggs_types,
        help="Select a DGGS type from the available options.",
    )
    parser.add_argument("dggrid_ids", nargs="+", help="Input DGGRID cell ID(s)")
    parser.add_argument("resolution", type=int, help="resolution")
    parser.add_argument(
        "input_address_type",
        choices=output_address_types,
        default="SEQNUM",
        nargs="?",  # This makes the argument optional
        help="Select an input address type from the available options.",
    )
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=True,
        help="Enable Antimeridian splitting",
    )
    parser.add_argument(
        "-aggregate",
        "--aggregate",
        action="store_true",
        help="Aggregate the resulting polygons",
    )
    parser.add_argument(
        "-options",
        "--options",
        type=str,
        default=None,
        help="JSON string of options to pass to grid_cell_polygons_from_cellids. "
             "Example: '{\"densification\": 2}'",
    )
    args = parser.parse_args()
    dggrid_instance = create_dggrid_instance()
    
    # Parse options JSON if provided
    options = None
    if args.options:
        try:
            options = json.loads(args.options)
        except json.JSONDecodeError as e:
            print(f"Error: Invalid JSON in options: {str(e)}")
            return
    
    geojson_data = json.dumps(
        dggrid2geojson(
            dggrid_instance,
            args.dggs_type,
            args.dggrid_ids,
            args.resolution,
            args.input_address_type,
            split_antimeridian=args.split_antimeridian,
            aggregate=args.aggregate,
            options=options,
        )
    )
    print(geojson_data)


if __name__ == "__main__":
    dggrid2geojson_cli()
