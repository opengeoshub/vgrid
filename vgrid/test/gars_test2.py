from vgrid.geocode.edgarsgrid import EDGARSGrid

# from latlon
ggrid = EDGARSGrid.from_latlon(-89.55, -179.57, resolution=3)

# from GARS ID
ggrid = EDGARSGrid("D01AA23")

print(ggrid)
# # get bounding poly
# grid_poly = ggrid.polygon

# # get GARS ID
# gars_id = str(ggrid)

# # UTM CRS EPSG Code
# epsg_code = ggrid.utm_epsg

from vgrid.geocode.gedgarsgrid import GEDGARSGrid

# from latlon
ggrid = GEDGARSGrid.from_latlon(-89.55, -179.57, resolution=30)

# # from GARS ID
# ggrid = GEDGARSGrid("GD1A")

# # get bounding poly
# grid_poly = ggrid.polygon

# # get GARS ID
# gars_id = str(ggrid)
