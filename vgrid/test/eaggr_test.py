from vgrid.utils.eaggr.eaggr import Eaggr
from vgrid.utils.eaggr.enums.model import Model
from vgrid.utils.eaggr.enums.dggs_shape_type import DggsShapeType
from vgrid.utils.eaggr.enums.dggs_shape_location import DggsShapeLocation
from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
from vgrid.utils.eaggr.exceptions.eaggr_exception import EaggrException
from vgrid.utils.eaggr.shapes.lat_long_point import LatLongPoint
from vgrid.utils.eaggr.shapes.lat_long_linestring import LatLongLinestring
from vgrid.utils.eaggr.shapes.lat_long_polygon import LatLongPolygon
from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
from vgrid.utils.eaggr.shapes.dggs_linestring import DggsLinestring
from vgrid.utils.eaggr.shapes.dggs_polygon import DggsPolygon
from vgrid.utils.eaggr.enums.dggs_analysis_type import DggsAnalysisType
from vgrid.utils.eaggr.shapes.dggs_shape import DggsShape

dggs = Eaggr(Model.ISEA3H)
dggs = Eaggr(Model.ISEA4T)

latitude, longitude = 10.775275567242561, 106.70679737574993# 
resolution =  10
 # Create the lat/long points
lat_long_point = LatLongPoint(10.775275567242561, 106.70679737574993, resolution)
# Initialise the DGGS model
# Convert the first lat/long point
dggs_cell = dggs.convert_point_to_dggs_cell(lat_long_point)
print (dggs_cell.get_cell_id())