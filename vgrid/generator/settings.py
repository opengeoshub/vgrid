from pyproj import Geod
geod = Geod(ellps="WGS84")
from shapely.geometry import mapping

max_cells = 1_000_000
chunk_size = 10_000

def graticule_dggs_metrics(cell_polygon):
    min_lon, min_lat, max_lon, max_lat = cell_polygon.bounds
    center_lat = round((min_lat + max_lat) / 2,7)
    center_lon = round((min_lon + max_lon) / 2,7)
    cell_width = round(geod.line_length([min_lon, max_lon], [min_lat, min_lat]),3)
    cell_height = round(geod.line_length([min_lon, min_lon], [min_lat, max_lat]),3)
    cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),3)  # Area in square meters     
    return center_lat, center_lon, cell_width, cell_height, cell_area

def geodesic_dggs_metrics(cell_polygon, num_edges):
    cell_centroid = cell_polygon.centroid
    center_lat =  round(cell_centroid.y, 7)
    center_lon = round(cell_centroid.x, 7)    
    cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),3) # Area in square meters     
    cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
    avg_edge_len = round(cell_perimeter / num_edges,3)    
    return center_lat, center_lon, avg_edge_len, cell_area

# Convert Graticule DGGS cell to GeoJSON feature
def graticule_dggs_to_feature(dggs_name, cell_id, resolution, cell_polygon):
    center_lat,center_lon,cell_width,cell_height,cell_area =  graticule_dggs_metrics(cell_polygon)
    feature = {
                "type": "Feature",
                "geometry":mapping(cell_polygon),
                "properties": {
                    f"{dggs_name}": str(cell_id),
                    "resolution": resolution,
                    "center_lat": center_lat,
                    "center_lon": center_lon,
                    "cell_width": cell_width,
                    "cell_height": cell_height,
                    "cell_area": cell_area
                }
            }
    return feature

# Convert Geodesic DGGS cell to GeoJSON feature
def geodesic_dggs_to_feature(dggs_name, cell_id, resolution, cell_polygon,num_edges):
    center_lat,center_lon,avg_edge_len,cell_area =  geodesic_dggs_metrics(cell_polygon,num_edges)
    feature = {
                "type": "Feature",
                "geometry":mapping(cell_polygon),
                "properties": {
                    f"{dggs_name}": str(cell_id),
                    "resolution": resolution,
                    "center_lat": center_lat,
                    "center_lon": center_lon,
                    "avg_edge_len": avg_edge_len,
                    "cell_area": cell_area
                }
            }
    return feature

isea4t_res_accuracy_dict = {
    0: 25_503_281_086_204.43,
    1: 6_375_820_271_551.114,    
    2: 1_593_955_067_887.7715,
    3: 398_488_766_971.94995,
    4: 99_622_191_742.98041,
    5: 24905_547_935.752182,
    6: 6_226_386_983.930966,
    7: 1_556_596_745.9898202,
    8: 389_149_186.4903765,
    9: 97_287_296.6296727,
    10: 24_321_824.150339592,
    11: 6_080_456.0446634805,
    12: 1_520_114.0040872877,
    13: 380_028.5081004044,
    14: 95_007.11994651864,
    15: 23_751.787065212124,
    16: 5_937.9396877205645,
    17: 1_484.492000512607,
    18: 371.1159215456855,
    19: 92.78605896888773,    
    20: 23.189436159755584,
    21: 5.804437622405244,
    22: 1.4440308231349632,
    23: 0.36808628825008866,
    24: 0.0849429895961743,
    25: 0.028314329865391435,
    
    26: 7.08*10**-3, # accuracy returns 0.0, avg_edge_len =  0.11562
    27: 1.77*10**-3, # accuracy returns 0.0, avg_edge_len =  0.05781
    28: 4.42*10**-4, # accuracy returns 0.0, avg_edge_len =  0.0289
    29: 1.11*10**-4, # accuracy returns 0.0, avg_edge_len =  0.01445
    30: 2.77*10**-5, # accuracy returns 0.0, avg_edge_len = 0.00723
    31: 6.91*10**-6, # accuracy returns 0.0, avg_edge_len =  0.00361
    32: 1.73*10**-6, # accuracy returns 0.0, avg_edge_len =  0.00181
    33: 5.76*10**-7, # accuracy returns 0.0, avg_edge_len = 0.0009
    34: 1.92*10**-7, # accuracy returns 0.0, avg_edge_len = 0.00045
    35: 6.40*10**-8, # accuracy returns 0.0, avg_edge_len = 0.00023
    36: 2.13*10**-8, # accuracy returns 0, avg_edge_len = 0.00011
    37: 7.11*10**-9, # accuracy returns 0.0, avg_edge_len = 6*10**(-5)
    38: 2.37*10**-9, # accuracy returns 0.0, avg_edge_len = 3*10**(-5)
    39: 7.90*10**-10 # accuracy returns 0.0, avg_edge_len = 10**(-5)
    }    

isea3h_res_accuracy_dict = {
        0: 25_503_281_086_204.43,
        1: 17_002_187_390_802.953,
        2: 5_667_395_796_934.327,
        3: 1_889_131_932_311.4424,
        4: 629_710_644_103.8047,
        5: 209_903_548_034.5921,
        6: 69_967_849_344.8546,
        7: 23_322_616_448.284866,
        8: 7_774_205_482.77106,
        9: 2_591_401_827.5809155,
        10: 863_800_609.1842003,
        11: 287_933_536.4041716,
        12: 95_977_845.45861907,
        13: 31_992_615.152873024,
        14: 10_664_205.060395785,
        15: 3_554_735.0295700384,
        16: 1_184_911.6670852362,
        17: 394_970.54625696875,
        18: 131_656.84875232293,
        19: 43_885.62568888426, 
        20: 14628.541896294753,
        21: 4_876.180632098251,
        22: 1_625.3841059227952,
        23: 541.7947019742651,
        24: 180.58879588146658,
        25: 60.196265293822194,
        26: 20.074859874562527,
        27: 6.6821818482323785,
        
        28: 2.2368320593659234,
        29: 0.7361725765001773,
        30: 0.2548289687885229,
        31: 0.0849429895961743,
        32: 0.028314329865391435,
       
        33: 0.009438109955130478,  
        34: 0.0031460366517101594,  
        35: 0.0010486788839033865,      
        36: 0.0003495596279677955, 
        37: 0.0001165198769892652,   
        38: 0.0000388399589964217,
        39: 0.0000129466529988072,      
        40: 0.0000043155509996024
    }    
    
isea3h_accuracy_res_dict = {
            25_503_281_086_204.43: 0,
            17_002_187_390_802.953: 1,
            5_667_395_796_934.327: 2,
            1_889_131_932_311.4424: 3,
            629_710_644_103.8047: 4,
            209_903_548_034.5921: 5,
            69_967_849_344.8546: 6,
            23_322_616_448.284866: 7,
            7_774_205_482.77106: 8,
            2_591_401_827.5809155: 9,
            863_800_609.1842003: 10,
            287_933_536.4041716: 11,
            95_977_845.45861907: 12,
            31_992_615.152873024: 13,
            10_664_205.060395785: 14,
            3_554_735.0295700384: 15,
            1_184_911.6670852362: 16,
            394_970.54625696875: 17,
            131_656.84875232293: 18,
            43_885.62568888426: 19,
            14628.541896294753: 20,
            4_876.180632098251: 21,
            1_625.3841059227952: 22,
            541.7947019742651: 23,
            180.58879588146658: 24,
            60.196265293822194: 25,
            20.074859874562527: 26,
            6.6821818482323785: 27,
            
            2.2368320593659234: 28,
            0.7361725765001773: 29,
            0.2548289687885229: 30,
            0.0849429895961743: 31,
            0.028314329865391435: 32,
            
            0.0: 33, # isea3h2point._accuracy always returns 0.0 from res 33
            0.0: 34,
            0.0: 35,
            0.0: 36,
            0.0: 37,
            0.0: 38,
            0.0: 39,
            0.0: 40
        }

