from shapely.geometry import box
import argparse
import platform

if platform.system() == "Linux":
    from vgrid.utils.dggrid4py import DGGRIDv7, dggs_types
    from vgrid.utils.dggrid4py.dggrid_runner import output_address_types


def generate_grid(dggrid_instance, dggs_type, resolution, bbox, address_type):
    def address_transform(dggrid_seqnum, dggs_type, resolution, address_type):
        address_type_transform = dggrid_instance.address_transform(
            [dggrid_seqnum],
            dggs_type=dggs_type,
            resolution=resolution,
            mixed_aperture_level=None,
            input_address_type="SEQNUM",
            output_address_type=address_type,
        )
        return address_type_transform.loc[0, address_type]

    if bbox:
        bounding_box = box(*bbox)
        dggrid_gdf = dggrid_instance.grid_cell_polygons_for_extent(
            dggs_type,
            resolution,
            clip_geom=bounding_box,
            split_dateline=True,
            output_address_type=address_type,
        )
        try:
            if address_type != "SEQNUM":
                dggrid_gdf["name"] = dggrid_gdf["name"].astype(str)
                dggrid_gdf["name"] = dggrid_gdf["name"].apply(
                    lambda val: address_transform(
                        val, dggs_type, resolution, address_type
                    )
                )
                dggrid_gdf = dggrid_gdf.rename(columns={"name": address_type.lower()})
            else:
                dggrid_gdf = dggrid_gdf.rename(columns={"name": "seqnum"})
        except Exception:
            pass
        geojson_path = f"dggrid_{dggs_type}_{resolution}_{address_type}_bbox.geojson"
        dggrid_gdf.to_file(geojson_path, driver="GeoJSON")

    else:
        dggrid_gdf = dggrid_instance.grid_cell_polygons_for_extent(
            dggs_type, resolution, split_dateline=True, output_address_type=address_type
        )
        try:
            if address_type != "SEQNUM":
                dggrid_gdf["name"] = dggrid_gdf["name"].astype(str)
                dggrid_gdf["name"] = dggrid_gdf["name"].apply(
                    lambda val: address_transform(
                        val, dggs_type, resolution, address_type
                    )
                )
                dggrid_gdf = dggrid_gdf.rename(columns={"name": address_type.lower()})
            else:
                dggrid_gdf = dggrid_gdf.rename(columns={"name": "seqnum"})
        except Exception:
            pass

        geojson_path = f"dggrid_{dggs_type}_{resolution}_{address_type}.geojson"
        dggrid_gdf.to_file(geojson_path, driver="GeoJSON")

    print(f"GeoJSON saved as {geojson_path}")


def main():
    if platform.system() == "Linux":
        parser = argparse.ArgumentParser(description="Generate DGGRID.")
        parser.add_argument(
            "-t",
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
            "--address_type",
            choices=output_address_types,
            help="Select an output address type.",
        )
        args = parser.parse_args()

        dggrid_instance = DGGRIDv7(
            executable="/usr/local/bin/dggrid",
            working_dir=".",
            capture_logs=False,
            silent=True,
            tmp_geo_out_legacy=False,
            debug=False,
        )
        resolution = args.resolution
        dggs_type = args.dggs_type
        bbox = args.bbox
        address_type = args.address_type
        try:
            generate_grid(dggrid_instance, dggs_type, resolution, bbox, address_type)
        except Exception:
            print(
                "Please ensure that -a <address_type> are set appropriately, and there is an excutable DGGRID located at /usr/local/bin/dggrid. Please install DGGRID following instructions from https://github.com/sahrk/DGGRID/blob/master/INSTALL.md"
            )
    else:
        print(
            "dggrid only works on Linux with an excutable DGGRID at /usr/local/bin/dggrid. Please install DGGRID following instructions from https://github.com/sahrk/DGGRID/blob/master/INSTALL.md"
        )


if __name__ == "__main__":
    main()
