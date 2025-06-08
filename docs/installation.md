# Installation

## Install from PyPI

**vgrid** is available on [PyPI](https://pypi.org/project/vgrid/). To install **vgrid**, run this command in your terminal:

```bash
pip install vgrid
```

## Install from conda-forge

**vgrid** is also available on [conda-forge](https://anaconda.org/conda-forge/vgrid). If you have
[Anaconda](https://www.anaconda.com/distribution/#download-section) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installed on your computer, you can install vgrid using the following command:

```bash
conda install vgrid -c conda-forge
```

The vgrid package has some optional dependencies (e.g., [geopandas](https://geopandas.org/) and [localtileserver](https://github.com/banesullivan/localtileserver)), which can be challenging to install on some computers, especially Windows. It is highly recommended that you create a fresh conda environment to install geopandas and vgrid. Follow the commands below to set up a conda env and install [geopandas](https://geopandas.org), [localtileserver](https://github.com/banesullivan/localtileserver), [keplergl](https://docs.kepler.gl/docs/keplergl-jupyter), [pydeck](https://deckgl.readthedocs.io/), and vgrid.

```bash
conda install -n base mamba -c conda-forge
mamba create -n geo vgrid geopandas localtileserver python -c conda-forge
```

Optionally, you can install some [Jupyter notebook extensions](https://github.com/ipython-contrib/jupyter_contrib_nbextensions), which can improve your productivity in the notebook environment. Some useful extensions include Table of Contents, Gist-it, Autopep8, Variable Inspector, etc. See this [post](https://towardsdatascience.com/jupyter-notebook-extensions-517fa69d2231) for more information.

```bash
conda install jupyter_contrib_nbextensions -c conda-forge
```

## Install from GitHub

To install the development version from GitHub using [Git](https://git-scm.com/), run the following command in your terminal:

```bash
pip install git+https://github.com/opengeoshub/vgrid
```


## Upgrade vgrid

If you have installed **vgrid** before and want to upgrade to the latest version, you can run the following command in your terminal:

```bash
pip install -U vgrid
```

If you use conda, you can update vgrid to the latest version by running the following command in your terminal:

```bash
conda update -c conda-forge vgrid
```

To install the development version from GitHub directly within Jupyter notebook without using Git, run the following code:

```python
import vgrid
vgrid.update_package()
```

