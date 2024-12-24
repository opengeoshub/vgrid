
# Reference: https://github.com/pistell/Leaflet.DumbMGRS
# https://mgrs-mapper.com/app, https://military-history.fandom.com/wiki/Military_Grid_Reference_System
# https://codesandbox.io/s/how-to-toggle-react-leaflet-layer-control-and-rectangle-grid-f43xi?file=/src/App.js:1057-1309
# https://github.com/GeoSpark/gridoverlay/tree/d66ed86636c7ec3f02ca2e9298ac3086c2023f1d
# https://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html#//00700000001n000000
# https://storymaps.arcgis.com/stories/842edf2b4381438b9a4edefed124775b
# https://github.com/dnlbaldwin/React-Leaflet-MGRS-Graticule
# https://dnlbaldwin.github.io/React-Leaflet-MGRS-Graticule/
# https://earth-info.nga.mil/index.php?dir=coordsys&action=mgrs-100km-polyline-dloads
# https://mgrs-data.org/metadata/
# https://ufl.maps.arcgis.com/apps/dashboards/2539764e24e74bd78f265f49c7adc2d1
# https://earth-info.nga.mil/index.php?dir=coordsys&action=gars-20x20-dloads

import json
from shapely.geometry import box, mapping
from tqdm import tqdm

bands = ['C', 'D', 'E', 'F', 'G', 'H', 'J', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X']

def generate_grid(polar=True):
    features = []

    def export_polygon(lon, lat, width, height, gzd):
        rect = box(lon, lat, lon + width, lat + height)
        features.append({
            'type': 'Feature',
            'geometry': mapping(rect),
            'properties': {
                'gzd': gzd
            }
        })

    if polar:
        export_polygon(-180, -90, 180, 10, 'A')
        export_polygon(0, -90, 180, 10, 'B')

    lat = -80
    for b in tqdm(bands, desc="Generating Bands"):
        if b == 'X':
            height = 12
            lon = -180
            for i in range(1, 31):
                gzd = '{:02d}{}'.format(i, b)
                width = 6
                export_polygon(lon, lat, width, height, gzd)
                lon += width
            export_polygon(lon, lat, 9, height, '31X')
            lon += 9
            export_polygon(lon, lat, 12, height, '33X')
            lon += 12
            export_polygon(lon, lat, 12, height, '35X')
            lon += 12
            export_polygon(lon, lat, 9, height, '37X')
            lon += 9
            for i in range(38, 61):
                gzd = '{:02d}{}'.format(i, b)
                width = 6
                export_polygon(lon, lat, width, height, gzd)
                lon += width
        else:
            height = 8
            lon = -180
            for i in range(1, 61):
                gzd = '{:02d}{}'.format(i, b)
                if b == 'V' and i == 31:
                    width = 3
                elif b == 'V' and i == 32:
                    width = 9
                else:
                    width = 6
                export_polygon(lon, lat, width, height, gzd)
                lon += width
        lat += height

    if polar:
        export_polygon(-180, 84, 180, 6, 'Y')
        export_polygon(0, 84, 180, 6, 'Z')

    return features


def main():
    try:
        # Generate gzd grid
        features = generate_grid()

        # Save to GeoJSON
        geojson = {
            'type': 'FeatureCollection',
            'features': features
        }
        with open('gzd.geojson', 'w') as f:
            json.dump(geojson, f)
        
        print("GeoJSON saved to gzd.geojson")

    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == '__main__':
    main()
