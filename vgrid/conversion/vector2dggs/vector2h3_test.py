from vgrid.conversion.vector2dggs.vector2h3 import vector2h3
import json

point = "./data/shape/point.geojson"
point_h3 = vector2h3(point, 5, topology=True)
# print(point_h3)

# Save the H3 conversion result as GeoJSON
with open("vector2h3.geojson", "w") as f:
    json.dump(point_h3, f, indent=2)

print("H3 conversion saved as vector2h3.geojson")
