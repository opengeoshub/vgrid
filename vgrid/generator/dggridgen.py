"""
DGGRID Grid Generator Module

Generates DGGRID grids for multiple grid types with automatic cell generation and validation using the DGGRID library.

Key Functions:
- generate_grid(): Core grid generation function with DGGRID instance
- dggridgen(): User-facing function with multiple output formats
- dggridgen_cli(): Command-line interface for grid generation
"""

from shapely.geometry import box
import argparse
import json
import geopandas as gpd

from dggrid4py import dggs_types
from dggrid4py.dggrid_runner import output_address_types
from vgrid.utils.io import (
    convert_to_output_format,
    create_dggrid_instance,
    is_full_world_bbox,
    validate_bbox,
    validate_dggrid_type,
    validate_dggrid_resolution,
)
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS


def _ensure_active_geometry(gdf):
    """Ensure GeoDataFrame has active geometry for geopandas ops."""
    if not isinstance(gdf, gpd.GeoDataFrame):
        gdf = gpd.GeoDataFrame(gdf, geometry="geometry")
    if not isinstance(gdf.geometry, gpd.GeoSeries):
        gdf = gdf.set_geometry("geometry")
    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:4326")
    return gdf


def generate_grid(
    dggrid_instance,
    dggs_type,
    resolution,
    bbox,
    output_address_type="SEQNUM",
    split_antimeridian=False,
    aggregate=False,
    options=None,
):
    dggs_type = validate_dggrid_type(dggs_type)
    resolution = validate_dggrid_resolution(dggs_type, resolution)
    bbox = validate_bbox(bbox)
    # DGGRID clip polygons that span the antimeridian fail on ISEA3H/ISEA4H etc.
    if bbox and is_full_world_bbox(bbox):
        bbox = None

    ### considering using dggrid_instance.grid_cellids_for_extent('ISEA4T', 10, output_address_type='SEQNUM')
    if bbox:
        bounding_box = box(*bbox)
        # dggrid4py can fail for clip_geom+bbox when split_dateline=True.
        # For bounded extents, antimeridian splitting is unnecessary in practice.
        kwargs = {
            "split_dateline": False,
            "output_address_type": output_address_type,
        }
        if options:
            kwargs.update(options)
        dggrid_gdf = dggrid_instance.grid_cell_polygons_for_extent(
            dggs_type,
            resolution,
            clip_geom=bounding_box,
            **kwargs,
        )

    else:
        kwargs = {
            "split_dateline": split_antimeridian,
            "output_address_type": output_address_type,
        }
        if options:
            kwargs.update(options)
        dggrid_gdf = dggrid_instance.grid_cell_polygons_for_extent(
            dggs_type,
            resolution,
            **kwargs,
        )
    dggrid_gdf = _ensure_active_geometry(dggrid_gdf)

    # Apply antimeridian fixing if requested
    if split_antimeridian:
        if aggregate:
            dggrid_gdf = dggrid_gdf.dissolve(by="global_id", as_index=False)
            dggrid_gdf = _ensure_active_geometry(dggrid_gdf)

    return dggrid_gdf


def dggridgen(
    dggrid_instance,
    dggs_type,
    resolution,
    bbox=None,
    output_address_type="SEQNUM",
    output_format="gpd",
    split_antimeridian=False,
    aggregate=False,
    options=None,
):
    """
    Generate DGGRID grid for pure Python usage.

    Args:
        dggrid_instance: DGGRID instance for grid operations
        dggs_type (str): DGGS type from dggs_types
        resolution (int): Resolution level
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat].
            Defaults to None (whole world).
        output_address_type (str, optional): Address type for output. Defaults to "SEQNUM".
        output_format (str, optional): Output format handled entirely by convert_to_output_format
        split_antimeridian (bool, optional): When True, apply antimeridian fixing to the resulting polygons.
            Defaults to False when None or omitted.
        options (dict, optional): Options to pass to grid_cell_polygons_for_extent. 
            For example: {"densification": 2} to add densification points.
            Defaults to None.

    Returns:
        Delegated to convert_to_output_format
    """
    gdf = generate_grid(
        dggrid_instance,
        dggs_type,
        resolution,
        bbox,
        output_address_type,
        split_antimeridian=split_antimeridian,
        aggregate=aggregate,
        options=options,
    )
    output_name = f"dggrid_{dggs_type}_{resolution}"
    return convert_to_output_format(gdf, output_format, output_name)


def dggridgen_cli():
    parser = argparse.ArgumentParser(description="Generate DGGRID.")
    parser.add_argument(
        "-dggs",
        "--dggs_type",
        choices=dggs_types,
        help="Select a DGGS type from the available options.",
    )
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="resolution"
    )
    parser.add_argument(
        "-b",
        "--bbox",
        type=float,
        nargs=4,
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)",
    )
    parser.add_argument(
        "-a",
        "--output_address_type",
        choices=output_address_types,
        default=None,
        help="Select an output address type.",
    )
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        choices=OUTPUT_FORMATS,
        default="gpd",
        help="Select an output format.",
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
        help="JSON string of options to pass to grid_cell_polygons_for_extent. "
             "Example: '{\"densification\": 2}'",
    )
    args = parser.parse_args()

    dggrid_instance = create_dggrid_instance()

    resolution = args.resolution
    dggs_type = args.dggs_type
    bbox = args.bbox
    output_address_type = args.output_address_type
    
    # Parse options JSON if provided
    options = None
    if args.options:
        try:
            options = json.loads(args.options)
        except json.JSONDecodeError as e:
            print(f"Error: Invalid JSON in options: {str(e)}")
            return

    try:
        result = dggridgen(
            dggrid_instance,
            dggs_type,
            resolution,
            bbox,
            output_address_type,
            args.output_format,
            split_antimeridian=args.split_antimeridian,
            aggregate=args.aggregate,
            options=options,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    dggridgen_cli()
