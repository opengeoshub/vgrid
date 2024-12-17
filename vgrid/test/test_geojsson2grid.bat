@echo off
REM Test geojson2grid with specified parameters
geojson2h3 -r 11 -geojson ./data/shape/multipolygon.geojson
geojson2s2 -r 18 -geojson ./data/shape/multipolygon.geojson
geojson2rhealpix -r 11 -geojson ./data/shape/multipolygon.geojson
geojson2eaggrisea4t -r 17 -geojson ./data/shape/multipolygon.geojson
pause
