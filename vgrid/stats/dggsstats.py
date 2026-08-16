"""
This module provides functions for generating comprehensive statistics across multiple DGGS types.
"""

import argparse
import glob
import math
import os

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import TwoSlopeNorm

# Import all the individual inspect functions
from vgrid.stats.dggridstats import dggridinspect
from vgrid.stats.a5stats import a5inspect
from vgrid.stats.dggalstats import dggalinspect

# Import utilities
from vgrid.utils.constants import DGGS_INSPECT
from vgrid.utils.io import create_dggrid_instance
import warnings

warnings.filterwarnings(
    "ignore",
    message="driver ESRI Shapefile does not support open option DRIVER",
    category=RuntimeWarning,
)

_Y_COLUMN_LABELS = {
    "norm_area": "Normalized Area",
    "ipq": "Compactness (IPQ)",
}

_Y_LIM_DEFAULTS = {
    "norm_area": (0.5, 1.3),
    "ipq": (0.5, 1.0),
}

_DGGS_TYPE_DISPLAY = {
    "healpix": "HEALPix",
    "rhealpix": "rHEALPix",
}


def _y_column_label(y_column: str) -> str:
    return _Y_COLUMN_LABELS.get(y_column, y_column)


def _default_y_lim(y_column: str) -> tuple:
    return _Y_LIM_DEFAULTS.get(y_column, (0.5, 1.3))


def _format_dggs_type_label(dggs_type: str) -> str:
    key = str(dggs_type).strip().lower()
    return _DGGS_TYPE_DISPLAY.get(key, key.upper())


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
        # "dggrid_igeo7": {
        #     "inspect_func": dggridinspect,
        #     "cell_id_col": "global_id",
        #     "dggs_type": "IGEO7",
        # },
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
        "dggal_rhealpix": {
             "inspect_func": dggalinspect,
             "cell_id_col": "rhealpix",
             "dggs_type": "rhealpix",         
        },
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
    parquet_folder: str = ".", y_column: str = "norm_area", y_lim: tuple | None = None
) -> pd.DataFrame:
    """
    Create a seaborn boxplot from DGGS inspection parquet files in a folder.

    Args:
        parquet_folder (str): Folder containing ``*.parquet`` inspection files
            (default: current folder ``"."``)
        y_column (str): Column name to plot on y-axis (default: ``"norm_area"``)
        y_lim (tuple): Y-axis limits as (min, max). Defaults to (0.5, 1.3) for
            ``norm_area`` and (0.5, 1.0) for ``ipq``.

    Returns:
        pd.DataFrame: Summary statistics dataframe grouped by DGGS type
    """
    if y_lim is None:
        y_lim = _default_y_lim(y_column)

    parquet_files = sorted(glob.glob(os.path.join(parquet_folder, "*.parquet")))
    if not parquet_files:
        raise ValueError(f"No parquet files found in folder: {parquet_folder}")

    print(f"Found {len(parquet_files)} parquet files in {parquet_folder}")

    # Read attributes only: geometry is unused here and mixing CRS variants
    # (e.g. "WGS 84" and "WGS 84 (CRS84)") breaks GeoDataFrame concatenation.
    frames = []
    for parquet_file in parquet_files:
        try:
            frame = pd.read_parquet(parquet_file)
            frame = frame.drop(columns=["geometry"], errors="ignore")
            frames.append(pd.DataFrame(frame))
            print(f"Loaded {parquet_file} with {len(frame)} cells")
        except Exception as e:
            print(f"Error reading {parquet_file}: {e}")
            continue

    if not frames:
        raise ValueError(f"No valid parquet files could be read from: {parquet_folder}")

    df = pd.concat(frames, ignore_index=True)

    # Filter to keep only cells that do NOT cross the dateline
    df = df[~df["crossed"]]

    # Format DGGS type labels for display (e.g. HEALPix, rHEALPix)
    df["dggs_type"] = df["dggs_type"].map(_format_dggs_type_label)

    data_min = float(df[y_column].min())
    data_max = float(df[y_column].max())

    print(
        f"Loaded data with {len(df)} cells across DGGS types: {df['dggs_type'].unique()}"
    )
    print(f"Resolution range: {df['resolution'].min()}-{df['resolution'].max()}")
    print(f"{y_column} max: {data_max}")
    print(f"{y_column} min: {data_min}")

    # Wider figure so long type names (HEALPix, rHEALPix, ...) have room
    n_types = df["dggs_type"].nunique()
    fig_width = max(10, 1.4 * n_types + 2)
    plt.figure(figsize=(fig_width, 9))

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
        hue="dggs_type",
        data=df,
        order=sorted(df["dggs_type"].unique(), key=str.casefold),
        palette="viridis",
        saturation=0.9,
        showfliers=True,
        flierprops=outlier_design,
        legend=False,
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

    plt.ylabel(_y_column_label(y_column), fontsize="xx-large")

    # Set min and max values for y-axis
    plt.ylim(y_lim)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.22)

    # Save to current directory with predefined filename
    output_file = f"dggs_{y_column}_box.png"
    plt.savefig(output_file, bbox_inches="tight", dpi=300)

    plt.show()

    print("Boxplot created successfully!")
    print(
        f"Data contains {len(df)} cells across DGGS types: {df['dggs_type'].unique()}"
    )
    print(f"Resolution range: {df['resolution'].min()}-{df['resolution'].max()}")
    print(f"{y_column} max: {data_max}")
    print(f"{y_column} min: {data_min}")

    # Print some summary statistics
    print("\nSummary statistics by DGGS type:")
    summary = df.groupby("dggs_type")[y_column].agg(
        ["count", "mean", "std", "max", "min"]
    )
    summary["max_min"] = summary["max"] / summary["min"]
    print(summary)

    csv_file = f"dggs_{y_column}_box.csv"
    summary.to_csv(csv_file)
    print(f"Saved summary to: {csv_file}")

    return summary


def dggsboxplot_cli():
    """
    Command-line interface for creating DGGS boxplots from inspection data.

    CLI options:
      --folder: Folder containing parquet files (default: current folder)
      --y-column: Column name to plot on y-axis (default: norm_area)
      --y-min: Minimum y-axis value (default: column-specific)
      --y-max: Maximum y-axis value (default: column-specific)
    """

    parser = argparse.ArgumentParser(
        description="Create boxplots from DGGS inspection parquet files in a folder"
    )
    parser.add_argument(
        "-folder",
        "--folder",
        type=str,
        default=".",
        help="Folder containing parquet inspection files (default: current folder)",
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
        default=None,
        help="Minimum y-axis value (default: 0.5)",
    )
    parser.add_argument(
        "-ymax",
        "--ymax",
        type=float,
        default=None,
        help="Maximum y-axis value (default: 1.3 for norm_area, 1.0 for ipq)",
    )

    args = parser.parse_args()

    try:
        y_lim = None
        if args.ymin is not None or args.ymax is not None:
            default_min, default_max = _default_y_lim(args.y_column)
            y_lim = (
                default_min if args.ymin is None else args.ymin,
                default_max if args.ymax is None else args.ymax,
            )
        dggsboxplot(
            parquet_folder=args.folder,
            y_column=args.y_column,
            y_lim=y_lim,
        )
    except Exception as e:
        print(f"Error: {e}")
        return None


def dggsboxplot_folder(
    folder: str = ".", y_column: str = "norm_area", y_lim: tuple | None = None
) -> pd.DataFrame:
    """Alias for ``dggsboxplot`` (kept for backward compatibility)."""
    return dggsboxplot(parquet_folder=folder, y_column=y_column, y_lim=y_lim)


def dggsboxplot_folder_cli():
    """Alias CLI for ``dggsboxplot_cli`` (kept for backward compatibility)."""
    return dggsboxplot_cli()


_WORLD_COUNTRIES_URL = (
    "https://raw.githubusercontent.com/opengeoshub/vopendata/"
    "refs/heads/main/shape/world_countries.geojson"
)


def _panel_title_from_gdf(gdf: gpd.GeoDataFrame, parquet_path: str) -> str:
    if "dggs_type" in gdf.columns and gdf["dggs_type"].notna().any():
        return _format_dggs_type_label(gdf["dggs_type"].dropna().iloc[0])
    stem = os.path.splitext(os.path.basename(parquet_path))[0]
    for prefix in ("dggal_", "dggrid_"):
        if stem.lower().startswith(prefix):
            stem = stem[len(prefix) :]
            break
    return _format_dggs_type_label(stem)


def _load_norm_area_panel(parquet_path: str, y_column: str) -> gpd.GeoDataFrame:
    columns = ["resolution", "crossed", y_column, "geometry"]
    # dggs_type is optional but preferred for panel titles
    try:
        gdf = gpd.read_parquet(parquet_path, columns=[*columns, "dggs_type"])
    except Exception:
        gdf = gpd.read_parquet(parquet_path, columns=columns)

    if "resolution" not in gdf.columns or gdf["resolution"].isna().all():
        raise ValueError(f"No resolution column in {parquet_path}")
    resolution = int(gdf["resolution"].max())
    gdf = gdf[gdf["resolution"] == resolution]
    if gdf.empty:
        raise ValueError(f"No cells at resolution={resolution} in {parquet_path}")

    if "crossed" in gdf.columns:
        gdf = gdf[~gdf["crossed"].fillna(False)]
    if gdf.crs is None:
        gdf = gdf.set_crs("EPSG:4326")
    return gdf


def dggs_norm_area_distribution(
    parquet_folder: str = ".",
    y_column: str = "norm_area",
    y_lim: tuple | None = None,
    crs: str = "proj=moll",
):
    """
    Plot Mollweide maps of area distortion from DGGS inspection parquet files.

    Each panel uses the maximum resolution present in its parquet file.

    Style matches the per-DGGS ``*_norm_area`` maps: cell polygons, Mollweide CRS,
    ``RdYlBu_r``, and country outlines. Color limits and the axis label follow
    ``dggsboxplot``.

    Args:
        parquet_folder: Folder containing ``*.parquet`` inspection files.
        y_column: Column to map (default: ``"norm_area"``).
        y_lim: Color limits as (min, max). Defaults to the same limits as
            ``dggsboxplot`` for ``y_column``.
        crs: Plot CRS (default: Mollweide).
    """
    if y_lim is None:
        y_lim = _default_y_lim(y_column)
    vmin, vmax = y_lim
    vcenter = 1.0 if y_column == "norm_area" else (vmin + vmax) / 2.0
    norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)

    parquet_files = sorted(glob.glob(os.path.join(parquet_folder, "*.parquet")))
    if not parquet_files:
        raise ValueError(f"No parquet files found in folder: {parquet_folder}")

    print(f"Found {len(parquet_files)} parquet files in {parquet_folder}")

    panels = []
    for parquet_file in parquet_files:
        try:
            gdf = _load_norm_area_panel(parquet_file, y_column)
            title = _panel_title_from_gdf(gdf, parquet_file)
            res_used = int(gdf["resolution"].iloc[0])
            print(
                f"{os.path.basename(parquet_file)} res={res_used}: {len(gdf)} cells  "
                f"{y_column} min={gdf[y_column].min():.4f} max={gdf[y_column].max():.4f}"
            )
            panels.append((gdf, title))
        except Exception as e:
            print(f"skip {parquet_file}: {e}")
            continue

    if not panels:
        raise ValueError(f"No valid panels could be built from: {parquet_folder}")

    n = len(panels)
    ncols = min(3, n)
    nrows = math.ceil(n / ncols)

    world = gpd.read_file(_WORLD_COUNTRIES_URL).to_crs(crs)

    # Size the figure from the panel geometry so equal-aspect Mollweide maps
    # (2:1) fill their axes exactly and leave no vertical slack.
    panel_width = 5.0
    left, right, wspace = 0.01, 0.99, 0.05
    title_height = 0.35
    row_gap = 0.55
    top_margin = 0.7
    bottom_margin = 1.0

    fig_width = panel_width * ncols
    axes_width = (right - left) * fig_width / (ncols + (ncols - 1) * wspace)
    axes_height = axes_width / 2.0
    row_pitch = title_height + row_gap
    fig_height = (
        top_margin
        + title_height
        + nrows * axes_height
        + (nrows - 1) * row_pitch
        + bottom_margin
    )

    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height))
    fig.subplots_adjust(
        left=left,
        right=right,
        top=1 - (top_margin + title_height) / fig_height,
        bottom=bottom_margin / fig_height,
        wspace=wspace,
        hspace=row_pitch / axes_height,
    )
    axes_list = np.atleast_1d(axes).ravel()

    for ax, (gdf, title) in zip(axes_list, panels):
        gdf.to_crs(crs).plot(
            column=y_column,
            ax=ax,
            norm=norm,
            cmap="RdYlBu_r",
            linewidth=0.05,
            edgecolor="face",
            legend=False,
        )
        world.boundary.plot(ax=ax, color=None, edgecolor="black", linewidth=0.2)
        ax.set_title(title, fontsize=11, fontweight="bold", pad=4)
        ax.axis("off")
        ax.margins(0)
        ax.set_aspect("equal")
        ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)

    for ax in axes_list[len(panels) :]:
        ax.set_visible(False)

    sm = plt.cm.ScalarMappable(cmap="RdYlBu_r", norm=norm)
    sm.set_array([])
    # Dedicated colorbar axes inside the reserved bottom margin, so it never
    # steals space from (or overlaps) the map panels.
    cbar_ax = fig.add_axes([0.25, 0.5 / fig_height, 0.5, 0.18 / fig_height])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cbar.set_label(_y_column_label(y_column), fontsize=14)
    cbar.ax.tick_params(labelsize=12)
    cbar.set_ticks(np.linspace(vmin, vmax, 9))

    fig.suptitle(
        "DGGS Area Distortion Distribution",
        fontsize=15,
        fontweight="bold",
        y=1 - 0.28 / fig_height,
    )

    output_file = f"dggs_{y_column}_distribution.png"
    fig.savefig(output_file, dpi=300)
    plt.show()
    print(f"Saved: {os.path.abspath(output_file)}")
    return output_file


def dggs_norm_area_distribution_cli():
    """
    Command-line interface for DGGS normalized-area distribution maps.

    Each panel uses the maximum resolution present in its parquet file.

    CLI options:
      --folder: Folder containing parquet files (default: current folder)
      --y-column: Column to map (default: norm_area)
      --y-min / --y-max: Color limits (default: column-specific)
    """
    parser = argparse.ArgumentParser(
        description="Plot DGGS area distortion maps from inspection parquet files"
    )
    parser.add_argument(
        "-folder",
        "--folder",
        type=str,
        default=".",
        help="Folder containing parquet inspection files (default: current folder)",
    )
    parser.add_argument(
        "-y",
        "--y-column",
        type=str,
        default="norm_area",
        help="Column to map (default: norm_area)",
    )
    parser.add_argument(
        "-ymin",
        "--ymin",
        type=float,
        default=None,
        help="Minimum color value (default: column-specific)",
    )
    parser.add_argument(
        "-ymax",
        "--ymax",
        type=float,
        default=None,
        help="Maximum color value (default: column-specific)",
    )
    args = parser.parse_args()

    try:
        y_lim = None
        if args.ymin is not None or args.ymax is not None:
            default_min, default_max = _default_y_lim(args.y_column)
            y_lim = (
                default_min if args.ymin is None else args.ymin,
                default_max if args.ymax is None else args.ymax,
            )
        dggs_norm_area_distribution(
            parquet_folder=args.folder,
            y_column=args.y_column,
            y_lim=y_lim,
        )
    except Exception as e:
        print(f"Error: {e}")
        return None


if __name__ == "__main__":
    dggsinspect_cli()
