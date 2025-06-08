# Welcome to vgrid

[![image](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/opengeoshub/vgrid/blob/master)
[![image](https://studiolab.sagemaker.aws/studiolab.svg)](https://studiolab.sagemaker.aws/import/github/opengeoshub/vgrid/blob/master/docs/notebooks/00_conversion.ipynb)
[![image](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/opengeoshub/vgrid/HEAD)
[![image](https://img.shields.io/pypi/v/vgrid.svg)](https://pypi.python.org/pypi/vgrid)
[![image](https://static.pepy.tech/badge/vgrid)](https://pepy.tech/project/vgrid)
[![image](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<div align="center">
<img src="https://raw.githubusercontent.com/opengeoshub/vgrid/master/docs/assets/logo_rect.png" alt="logo">
</div>

**A Python package for geospatial analysis and interactive mapping in a Jupyter environment.**

-   GitHub repo: <https://github.com/opengeoshub/vgrid>
-   Documentation: <https://vgrid.org>
-   PyPI: <https://pypi.org/project/vgrid>
-   Conda-forge: <https://anaconda.org/conda-forge/vgrid>
-   vgrid tutorials on YouTube: <https://youtube.com/@giswqs>
-   Free software: [MIT license](https://opensource.org/licenses/MIT)

Join our Discord server 👇

[![](https://dcbadge.limes.pink/api/server/https://discord.gg/UgZecTUq5P)](https://discord.gg/UgZecTUq5P)

## Introduction

**vgrid** is a Python package for interactive mapping and geospatial analysis with minimal coding in a Jupyter environment. It is a spin-off project of the [geemap](https://geemap.org) Python package, which was designed specifically to work with [Google Earth Engine](https://earthengine.google.com) (GEE). However, not everyone in the geospatial community has access to the GEE cloud computing platform. vgrid is designed to fill this gap for non-GEE users. It is a free and open-source Python package that enables users to analyze and visualize geospatial data with minimal coding in a Jupyter environment, such as Google Colab, Jupyter Notebook, JupyterLab, and [marimo](https://github.com/marimo-team/marimo). vgrid is built upon several open-source packages, such as [folium](https://github.com/python-visualization/folium) and [ipyleaflet](https://github.com/jupyter-widgets/ipyleaflet) (for creating interactive maps), [WhiteboxTools](https://github.com/jblindsay/whitebox-tools) and [whiteboxgui](https://github.com/opengeos/whiteboxgui) (for analyzing geospatial data), and [ipywidgets](https://github.com/jupyter-widgets/ipywidgets) (for designing interactive graphical user interfaces [GUIs]). vgrid has a toolset with various interactive tools that allow users to load vector and raster data onto the map without coding. In addition, users can use the powerful analytical backend (i.e., WhiteboxTools) to perform geospatial analysis directly within the vgrid user interface without writing a single line of code. The WhiteboxTools library currently contains **500+** tools for advanced geospatial analysis, such as [GIS Analysis](https://jblindsay.github.io/wbt_book/available_tools/gis_analysis.html), [Geomorphometric Analysis](https://jblindsay.github.io/wbt_book/available_tools/geomorphometric_analysis.html), [Hydrological Analysis](https://jblindsay.github.io/wbt_book/available_tools/hydrological_analysis.html), [LiDAR Data Analysis](https://jblindsay.github.io/wbt_book/available_tools/lidar_tools.html), [Mathematical and Statistical Analysis](https://jblindsay.github.io/wbt_book/available_tools/mathand_stats_tools.html), and [Stream Network Analysis](https://jblindsay.github.io/wbt_book/available_tools/stream_network_analysis.html).

## Acknowledgments

This project is supported by Amazon Web Services ([AWS](https://aws.amazon.com)).

## Statement of Need

There is a plethora of Python packages for geospatial analysis, such as [geopandas](https://geopandas.org) for vector data analysis and [xarray](https://docs.xarray.dev) for raster data analysis. As listed at [pyviz.org](https://pyviz.org), there are also many options for plotting data on a map in Python, ranging from libraries focused specifically on maps like [ipyleaflet](https://ipyleaflet.readthedocs.io) and [folium](https://python-visualization.github.io/folium) to general-purpose plotting tools that also support geospatial data types, such as [hvPlot](https://hvplot.pyviz.org), [bokeh](http://bokeh.org), and [plotly](https://plotly.com/python). While these tools provide powerful capabilities, displaying geospatial data from different file formats on an interactive map and performing basic analyses can be challenging, especially for users with limited coding skills. Furthermore, many tools lack bi-directional communication between the frontend (browser) and the backend (Python), limiting their interactivity and usability for exploring map data.

vgrid addresses these challenges by leveraging the bidirectional communication provided by ipyleaflet, enabling users to load and visualize geospatial datasets with just one line of code. vgrid also provides an interactive graphical user interface (GUI) for loading geospatial datasets without any coding. It is designed for anyone who wants to analyze and visualize geospatial data interactively in a Jupyter environment, making it particularly accessible for novice users with limited programming skills. Advanced programmers can also benefit from vgrid for geospatial data analysis and building interactive web applications.

## Usage

Launch the interactive notebook tutorial for the **vgrid** Python package with Amazon SageMaker Studio Lab, Microsoft Planetary Computer, Google Colab, or Binder:

[![image](https://studiolab.sagemaker.aws/studiolab.svg)](https://studiolab.sagemaker.aws/import/github/opengeoshub/vgrid/blob/master/examples/notebooks/00_key_features.ipynb)
[![image](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/opengeoshub/vgrid/blob/master)
[![image](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/opengeoshub/vgrid/HEAD)

Check out this excellent article on Medium - [vgrid a new Python Package for Geospatial data science](https://link.medium.com/HRRKDcynYgb)

## Key Features

vgrid offers a wide range of features and capabilities that empower geospatial data scientists, researchers, and developers to unlock the potential of their data. Some of the key features include:

-   **Creating an interactive map with just one line of code:** vgrid makes it easy to create an interactive map by providing a simple API that allows you to load and visualize geospatial datasets with minimal coding.

-   **Switching between different mapping backends:** vgrid supports multiple mapping backends, including ipyleaflet, folium, kepler.gl, pydeck, and bokeh. You can switch between these backends to create maps with different visualization styles and capabilities.

-   **Changing basemaps interactively:** vgrid allows you to change basemaps interactively, providing a variety of options such as OpenStreetMap, Stamen Terrain, CartoDB Positron, and many more.

-   **Adding XYZ, WMS, and vector tile services:** You can easily add XYZ, WMS, and vector tile services to your map, allowing you to overlay additional geospatial data from various sources.

-   **Displaying vector data:** vgrid supports various vector data formats, including Shapefile, GeoJSON, GeoPackage, and any vector format supported by GeoPandas. You can load and display vector data on the map, enabling you to visualize and analyze spatial features.

-   **Displaying raster data:** vgrid allows you to load and display raster data, such as GeoTIFFs, on the map. This feature is useful for visualizing satellite imagery, digital elevation models, and other gridded datasets.

-   **Creating custom legends and colorbars:** vgrid provides tools for customizing legends and colorbars on the map, allowing you to represent data values with different colors and corresponding labels.

-   **Creating split-panel maps and linked maps:** With vgrid, you can create split-panel maps to compare different datasets side by side. You can also create linked maps that synchronize interactions between multiple maps, providing a coordinated view of different spatial data.

-   **Downloading and visualizing OpenStreetMap data:** vgrid allows you to download and visualize OpenStreetMap data, providing access to detailed street maps, buildings, and other points of interest.

-   **Creating and editing vector data interactively:** vgrid includes tools for creating and editing vector data interactively on the map. You can draw points, lines, and polygons, and modify them as needed.

-   **Searching for geospatial data:** vgrid provides functionality for searching and accessing geospatial data from sources such as SpatialTemporal Asset Catalogs (STAC), Microsoft Planetary Computer, AWS Open Data Registry, and OpenAerialMap.

-   **Inspecting pixel values interactively:** vgrid allows you to interactively inspect pixel values in raster datasets, helping you analyze and understand the data at a more granular level.

-   **Creating choropleth maps and heat maps:** vgrid supports the creation of choropleth maps, where colors represent different data values for specific geographic areas. You can also create heat maps to visualize data density.

-   **Displaying data from a PostGIS database:** vgrid provides tools for connecting to a PostGIS database and displaying spatial data stored in the database on the map.

-   **Creating time series animations:** vgrid enables the creation of time series animations from both vector and raster data, allowing you to visualize temporal changes in your geospatial datasets.

-   **Analyzing geospatial data with whitebox:** vgrid integrates with WhiteboxTools and whiteboxgui, providing a suite of geospatial analyses, such as hydrological analysis, terrain analysis, and LiDAR processing.

-   **Segmenting and classifying remote sensing imagery:** vgrid integrates the segment-geospatial package, which provides tools for segmenting and classifying remote sensing imagery using deep learning algorithms.

-   **Building interactive web apps:** vgrid supports the development of interactive web applications using frameworks like Voila, Streamlit, and Solara. This allows you to share your geospatial analyses and visualizations with others in a user-friendly web interface.

These features and capabilities make vgrid a powerful tool for geospatial data exploration, analysis, and visualization. Whether you are a beginner or an experienced geospatial data scientist, vgrid provides an accessible and efficient way to work with geospatial data in Python.

## Citations

If you find **vgrid** useful in your research, please consider citing the following paper to support my work. Thank you for your support.

-   Wu, Q. (2021). vgrid: A Python package for interactive mapping and geospatial analysis with minimal coding in a Jupyter environment. _Journal of Open Source Software_, 6(63), 3414. <https://doi.org/10.21105/joss.03414>

## Demo

![](https://assets.gishub.org/images/vgrid_demo.gif)

## YouTube Channel

I have created a [YouTube Channel](https://youtube.com/@giswqs) for sharing geospatial tutorials. You can subscribe to my channel for regular updates. Check out the following videos for 3D mapping with MapLibre and vgrid.

[![MapLibre tutorials](https://assets.gishub.org/images/maplibre-tutorials.png)](https://bit.ly/maplibre)
