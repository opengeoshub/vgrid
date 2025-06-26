import argparse
import platform

if platform.system() == "Linux":
    from vgrid.utils.dggrid4py import DGGRIDv7, dggs_types


def dggrid_stats(dggrid_instance, dggs_type, resolution, output):
    dggrid_metrics = dggrid_instance.grid_stats_table(dggs_type, resolution)
    if output:
        dggrid_metrics.to_csv(output, index=False)
    else:
        print(dggrid_metrics)


def main():
    if platform.system() == "Linux":
        parser = argparse.ArgumentParser(description="Export or display DGGRID stats.")
        parser.add_argument(
            "-t",
            "--dggs_type",
            choices=dggs_types,
            help="Select a DGGS type from the available options.",
        )
        parser.add_argument(
            "-r", "--resolution", type=int, required=True, help="resolution"
        )
        parser.add_argument("-o", "--output", help="Output CSV file name.")
        args = parser.parse_args()

        dggrid_instance = DGGRIDv7(
            executable="/usr/local/bin/dggrid",
            working_dir=".",
            capture_logs=False,
            silent=True,
            tmp_geo_out_legacy=False,
            debug=False,
        )
        dggs_type = args.dggs_type
        resolution = args.resolution
        output = args.output
        try:
            dggrid_stats(dggrid_instance, dggs_type, resolution, output)
        except Exception:
            print(
                "Please ensure that the -r <resolution> are set appropriately, and there is an excutable DGGRID located at /usr/local/bin/dggrid. Please install DGGRID following instructions from https://github.com/sahrk/DGGRID/blob/master/INSTALL.md"
            )
    else:
        print(
            "dggrid only works on Linux with an excutable DGGRID at /usr/local/bin/dggrid. Please install DGGRID following instructions from https://github.com/sahrk/DGGRID/blob/master/INSTALL.md"
        )


if __name__ == "__main__":
    main()
