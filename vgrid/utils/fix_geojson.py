import json
from typing import Dict, Any
from vgrid.utils.antimeridian import fix_geojson

# Load GeoJSON file
with open("./vgrid/utils/polyhedra/equi7_t6_proj.geojson", "r") as file:
    geojson_data = json.load(file)

# Call the function with the loaded GeoJSON data
fixed_geojson = fix_geojson(
    geojson_data,
    force_north_pole=False,  # Set these based on your needs
    force_south_pole=False,
    fix_winding=True
)

# Optionally, save the fixed GeoJSON to a new file
with open("./vgrid/utils/polyhedra/equi7_t6_proj_fixed.geojson", "w") as output_file:
    json.dump(fixed_geojson, output_file, indent=2)
