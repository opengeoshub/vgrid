"""
This module provides functions for generating comprehensive statistics across multiple DGGS types.
"""

import argparse
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import seaborn as sns

# Import all the individual inspect functions
from vgrid.stats.h3stats import h3inspect
from vgrid.stats.s2stats import s2inspect
from vgrid.stats.a5stats import a5inspect
from vgrid.stats.isea4tstats import isea4tinspect
from vgrid.stats.rhealpixstats import rhealpixinspect
from vgrid.stats.dggalstats import dggalinspect
from vgrid.stats.dggridstats import dggridinspect

# Import utilities
from vgrid.utils.constants import DGGS_INSPECT
from vgrid.utils.io import create_dggrid_instance
import warnings

warnings.filterwarnings(
    "ignore",
    message="driver ESRI Shapefile does not support open option DRIVER",
    category=RuntimeWarning,
)


def dggsinspect():
    """
    Multi-DGGS cell inspection using DGGS_INSPECT configuration.

    Returns:
        dict: Dictionary with DGGS types as keys and GeoDataFrames as values
    """

    # Define DGGS type configurations with their inspect functions
    # All inspect functions take a `resolution` parameter, so there is no need
    # to configure a custom parameter name here.
    dggs_configs = {
        # "h3": {"inspect_func": h3inspect, "cell_id_col": "h3"},
        # "s2": {"inspect_func": s2inspect, "cell_id_col": "s2"},
        # "a5": {"inspect_func": a5inspect, "cell_id_col": "a5"},
        # "isea4t": {
        #     "inspect_func": isea4tinspect,
        #     "cell_id_col": "isea4t",
        # },
        # "rhealpix": {
        #     "inspect_func": rhealpixinspect,
        #     "cell_id_col": "rhealpix",
        # },
        # "dggrid_isea7h": {
        #     "inspect_func": dggridinspect,
        #     "cell_id_col": "global_id",
        #     "dggs_type": "ISEA7H",
        # },
        # "dggrid_fuller7h": {
        #     "inspect_func": dggridinspect,
        #     "cell_id_col": "global_id",
        #     "dggs_type": "FULLER7H",
        # },
        # "dggrid_isea4d": {
        #     "inspect_func": dggridinspect,
        #     "cell_id_col": "global_id",
        #     "dggs_type": "ISEA4D",
        # },
        # "dggrid_fuller4d": {
        #     "inspect_func": dggridinspect,
        #     "cell_id_col": "global_id",
        #     "dggs_type": "FULLER4D",
        # },
        # "dggrid_isea4t": {
        #     "inspect_func": dggridinspect,
        #     "cell_id_col": "global_id",
        #     "dggs_type": "ISEA4T",
        # },
        # "dggrid_fuller4t": {
        #     "inspect_func": dggridinspect,
        #     "cell_id_col": "global_id",
        #     "dggs_type": "FULLER4T",
        # },
        "dggrid_igeo7": {
            "inspect_func": dggridinspect,
            "cell_id_col": "global_id",
            "dggs_type": "IGEO7",
        },
        # "dggal_ivea3h": {
        #     "inspect_func": dggalinspect,
        #     "cell_id_col": "dggal_ivea3h",
        #     "dggs_type": "ivea3h",
        # },
        #   "dggal_ivea4r": {
        #     "inspect_func": dggalinspect,
        #     "cell_id_col": "dggal_ivea4r",
        #     "dggs_type": "ivea4r",
        # },
        #   "dggal_ivea7h": {
        #     "inspect_func": dggalinspect,
        #     "cell_id_col": "dggal_ivea7h",
        #     "dggs_type": "ivea7h",
        # },
        # "dggal_ivea9r": {
        #     "inspect_func": dggalinspect,
        #     "cell_id_col": "dggal_ivea9r",
        #     "dggs_type": "ivea9r",
        # },
    }

    # Dictionary to store processed GeoDataFrames for return
    processed_gdfs = {}

    for dggs_type, config in dggs_configs.items():
        if dggs_type not in DGGS_INSPECT:
            print(f"Warning: {dggs_type} not found in DGGS_INSPECT configuration")
            continue

        inspect_config = DGGS_INSPECT[dggs_type]
        min_res = inspect_config["min_res"]
        max_res = inspect_config["max_res"]

        print(f"Processing {dggs_type} for resolutions {min_res}-{max_res}")

        dggs_type_gdfs = []

        for res in range(min_res, max_res + 1):
            try:
                # Call the specific inspect function with appropriate parameters
                if dggs_type.startswith("dggrid_"):
                    # Create dggrid instance once for all dggrid operations
                    dggrid_instance = create_dggrid_instance()

                    # For dggrid functions that need dggrid_instance and dggs_type parameters
                    gdf = config["inspect_func"](
                        dggrid_instance,
                        dggs_type=config["dggs_type"],
                        resolution=res,
                    )
                elif "dggs_type" in config:
                    # For dggal functions that need dggs_type parameter
                    gdf = config["inspect_func"](
                        dggs_type=config["dggs_type"], resolution=res
                    )
                else:
                    # For standard inspect functions that take a `resolution` parameter
                    gdf = config["inspect_func"](resolution=res)

                # Add dggs_type column
                gdf["dggs_type"] = dggs_type

                # Rename the cell ID column to a generic name
                cell_id_col = config["cell_id_col"]
                if cell_id_col in gdf.columns:
                    gdf = gdf.rename(columns={cell_id_col: "cell_id"})

                # Ensure all expected columns exist, add NaN for missing ones
                expected_columns = [
                    "dggs_type",
                    "resolution",
                    "cell_id",
                    "geometry",
                    "cell_area",
                    "cell_perimeter",
                    "crossed",
                    "norm_area",
                    "ipq",
                    "zsc",
                    "cvh",
                ]

                for col in expected_columns:
                    if col not in gdf.columns:
                        gdf[col] = None

                # Reorder columns to match expected format
                gdf = gdf[expected_columns]

                dggs_type_gdfs.append(gdf)

            except Exception as e:
                print(f"Error processing {dggs_type} at resolution {res}: {e}")
                continue

        # Process and save this DGGS type immediately after completing all resolutions
        if dggs_type_gdfs:
            # Combine GeoDataFrames for this DGGS type
            combined_gdf = pd.concat(dggs_type_gdfs, ignore_index=True)

            # Filter to keep only cells that do NOT cross the dateline
            #  combined_gdf = combined_gdf[~combined_gdf["crossed"]]

            # Ensure it's a GeoDataFrame
            if not isinstance(combined_gdf, gpd.GeoDataFrame):
                combined_gdf = gpd.GeoDataFrame(combined_gdf, geometry="geometry")

            # Save this DGGS type immediately as a geoparquet file
            output_file = f"{dggs_type}.parquet"
            print(f"✓ Completed {dggs_type}: {len(combined_gdf)} cells")
            print(f"  Saving to: {output_file}")
            combined_gdf.to_parquet(output_file, index=False)
            print(
                f"  ✓ Successfully saved {len(combined_gdf)} records to {output_file}"
            )

            # Store in processed_gdfs for return
            processed_gdfs[dggs_type] = combined_gdf
        else:
            print(f"Warning: No valid data generated for {dggs_type}")

    if not processed_gdfs:
        raise ValueError(
            "No valid DGGS data could be generated for the specified resolution range"
        )

    print("\n🎉 All DGGS types processed and saved!")
    print(f"Total DGGS types: {len(processed_gdfs)}")
    total_cells = sum(len(gdf) for gdf in processed_gdfs.values())
    print(f"Total cells across all types: {total_cells}")

    return processed_gdfs


def dggsinspect_cli():
    """
    Command-line interface for multi-DGGS cell inspection using DGGS_INSPECT configuration.
    """
    try:
        results = dggsinspect()
        return results
    except Exception as e:
        print(f"Error: {e}")
        return None


def dggsboxplot(
    parquet_file: str, y_column: str = "norm_area", y_lim: tuple = (0.5, 1.3)
) -> pd.DataFrame:
    """
    Create a seaborn boxplot from an existing DGGS inspection parquet file.

    Args:
        parquet_file (str): Path to the input parquet file containing DGGS inspection data
        y_column (str): Column name to plot on y-axis (default: "norm_area")
        y_lim (tuple): Y-axis limits as (min, max) tuple (default: (0.5, 1.3))

    Returns:
        pd.DataFrame: Summary statistics dataframe grouped by DGGS type
    """

    # Read the existing parquet file
    gdf = gpd.read_parquet(parquet_file)
    # Recalculate cvh by Cartersian to avoid cvh > 1
    # gdf["cvh"] = gdf["geometry"].area / gdf["geometry"].convex_hull.area
    # Convert to regular DataFrame (drop geometry column for plotting)
    df = pd.DataFrame(gdf.drop(columns=["geometry"]))

    # Filter to keep only cells that do NOT cross the dateline
    df = df[~df["crossed"]]

    # Convert dggs_type to uppercase for display and sort by dggs_type
    df["dggs_type"] = df["dggs_type"].str.upper()

    print(
        f"Loaded data with {len(df)} cells across DGGS types: {df['dggs_type'].unique()}"
    )
    print(f"Resolution range: {df['resolution'].min()}-{df['resolution'].max()}")

    # Create the boxplot
    plt.figure(figsize=(9, 9))

    # Use modern seaborn style
    plt.style.use("default")
    sns.set_style("whitegrid")

    # Define design of the outliers
    outlier_design = dict(
        marker="o",
        markerfacecolor="black",
        markersize=1,
        linestyle="none",
        markeredgecolor="black",
    )

    # Plot the boxplots
    chart = sns.boxplot(
        x="dggs_type",
        y=y_column,
        data=df,
        palette="viridis",
        saturation=0.9,
        showfliers=True,
        flierprops=outlier_design,
    )

    plt.xticks(
        rotation=90,
        horizontalalignment="center",
        fontweight="light",
        fontsize="xx-large",
    )

    plt.xlabel("", fontsize="x-large")

    plt.yticks(
        rotation=0,
        horizontalalignment="right",
        fontweight="light",
        fontsize="xx-large",
    )

    plt.ylabel(y_column, fontsize="xx-large")

    # Set min and max values for y-axis
    plt.ylim(y_lim)

    plt.tight_layout()

    # Save to current directory with predefined filename
    output_file = f"dggs_{y_column}_box.png"
    plt.savefig(output_file, bbox_inches="tight", dpi=300)

    plt.show()

    print("Boxplot created successfully!")
    print(
        f"Data contains {len(df)} cells across DGGS types: {df['dggs_type'].unique()}"
    )
    print(f"Resolution range: {df['resolution'].min()}-{df['resolution'].max()}")

    # Print some summary statistics
    print("\nSummary statistics by DGGS type:")
    summary = df.groupby("dggs_type")[y_column].agg(
        ["count", "mean", "std", "min", "max"]
    )
    print(summary)

    return summary


def dggsboxplot_cli():
    """
    Command-line interface for creating DGGS boxplots from inspection data.

    CLI options:
      --input: Input parquet file path (required)
      --y-column: Column name to plot on y-axis (default: norm_area)
      --y-min: Minimum y-axis value (default: 0.5)
      --y-max: Maximum y-axis value (default: 1.3)
    """

    parser = argparse.ArgumentParser(
        description="Create boxplots from DGGS inspection data"
    )
    parser.add_argument(
        "-input",
        "--input",
        type=str,
        required=True,
        help="Input parquet file path containing DGGS inspection data",
    )
    parser.add_argument(
        "-y",
        "--y-column",
        type=str,
        default="norm_area",
        help="Column name to plot on y-axis (default: norm_area)",
    )
    parser.add_argument(
        "-ymin",
        "--ymin",
        type=float,
        default=0.5,
        help="Minimum y-axis value (default: 0.5)",
    )
    parser.add_argument(
        "-ymax",
        "--ymax",
        type=float,
        default=1.3,
        help="Maximum y-axis value (default: 1.3)",
    )

    args = parser.parse_args()

    try:
        dggsboxplot(
            parquet_file=args.input,
            y_column=args.y_column,
            y_lim=(args.ymin, args.ymax),
        )
    except Exception as e:
        print(f"Error: {e}")
        return None


def dggsboxplot_folder(
    folder: str = ".", y_column: str = "norm_area", y_lim: tuple = (0.5, 1.3)
) -> pd.DataFrame:
    """
    Create a seaborn boxplot from DGGS inspection parquet files in a folder.

    Args:
        folder (str): Path to folder containing parquet files (default: current folder ".")
        y_column (str): Column name to plot on y-axis (default: "norm_area")
        y_lim (tuple): Y-axis limits as (min, max) tuple (default: (0.5, 1.3))

    Returns:
        pd.DataFrame: Summary statistics dataframe grouped by DGGS type
    """

    import os
    import glob

    all_gdfs = []

    # Find all parquet files in the specified folder
    parquet_pattern = os.path.join(folder, "*.parquet")
    parquet_files = glob.glob(parquet_pattern)

    if not parquet_files:
        raise ValueError(f"No parquet files found in folder: {folder}")

    print(f"Found {len(parquet_files)} parquet files in {folder}")

    # Read all parquet files
    for parquet_file in parquet_files:
        try:
            gdf = gpd.read_parquet(parquet_file)
            all_gdfs.append(gdf)
            print(f"Loaded {parquet_file} with {len(gdf)} cells")
        except Exception as e:
            print(f"Error reading {parquet_file}: {e}")
            continue

    if not all_gdfs:
        raise ValueError("No valid parquet files could be read")

    # Combine all GeoDataFrames
    combined_gdf = pd.concat(all_gdfs, ignore_index=True)

    # Convert to regular DataFrame (drop geometry column for plotting)
    df = pd.DataFrame(combined_gdf.drop(columns=["geometry"]))
    # Filter to keep only cells that do NOT cross the dateline
    df = df[~df["crossed"]]
    # Convert dggs_type to uppercase for display
    df["dggs_type"] = df["dggs_type"].str.upper()

    print(
        f"Combined data with {len(df)} cells across DGGS types: {df['dggs_type'].unique()}"
    )
    print(f"Resolution range: {df['resolution'].min()}-{df['resolution'].max()}")

    # Create the boxplot
    plt.figure(figsize=(9, 9))

    # Use modern seaborn style
    plt.style.use("default")
    sns.set_style("whitegrid")

    # Define design of the outliers
    outlier_design = dict(
        marker="o",
        markerfacecolor="black",
        markersize=1,
        linestyle="none",
        markeredgecolor="black",
    )

    # Plot the boxplots
    chart = sns.boxplot(
        x="dggs_type",
        y=y_column,
        data=df,
        palette="viridis",
        saturation=0.9,
        showfliers=True,
        flierprops=outlier_design,
    )

    plt.xticks(
        rotation=90,
        horizontalalignment="center",
        fontweight="light",
        fontsize="xx-large",
    )

    plt.xlabel("", fontsize="x-large")

    plt.yticks(
        rotation=0,
        horizontalalignment="right",
        fontweight="light",
        fontsize="xx-large",
    )

    plt.ylabel(y_column, fontsize="xx-large")

    # Set min and max values for y-axis
    plt.ylim(y_lim)

    plt.tight_layout()

    # Save to current directory with predefined filename
    output_file = "box_plot_area_folder.png"
    plt.savefig(output_file, bbox_inches="tight", dpi=300)

    plt.show()

    print("Boxplot created successfully!")
    print(
        f"Data contains {len(df)} cells across DGGS types: {df['dggs_type'].unique()}"
    )
    print(f"Resolution range: {df['resolution'].min()}-{df['resolution'].max()}")

    # Print some summary statistics
    print("\nSummary statistics by DGGS type:")
    summary = df.groupby("dggs_type")[y_column].agg(
        ["count", "mean", "std", "min", "max"]
    )
    print(summary)

    return summary


def dggsboxplot_folder_cli():
    """
    Command-line interface for creating DGGS boxplots from inspection data files in a folder.

    CLI options:
      --folder: Folder path containing parquet files (default: current folder)
      --y-column: Column name to plot on y-axis (default: norm_area)
      --y-min: Minimum y-axis value (default: 0.5)
      --y-max: Maximum y-axis value (default: 1.3)
    """

    parser = argparse.ArgumentParser(
        description="Create boxplots from DGGS inspection data files in a folder"
    )
    parser.add_argument(
        "-folder",
        "--folder",
        type=str,
        default=".",
        help="Folder path containing parquet files (default: current folder)",
    )
    parser.add_argument(
        "-y",
        "--y-column",
        type=str,
        default="norm_area",
        help="Column name to plot on y-axis (default: norm_area)",
    )
    parser.add_argument(
        "-ymin",
        "--ymin",
        type=float,
        default=0.5,
        help="Minimum y-axis value (default: 0.5)",
    )
    parser.add_argument(
        "-ymax",
        "--ymax",
        type=float,
        default=1.3,
        help="Maximum y-axis value (default: 1.3)",
    )

    args = parser.parse_args()

    try:
        dggsboxplot_folder(
            folder=args.folder, y_column=args.y_column, y_lim=(args.ymin, args.ymax)
        )
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    dggsinspect_cli()
