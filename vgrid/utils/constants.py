"""
Constants module for vgrid.

This module contains constants and dictionaries used across the vgrid package
"""

import math

MAX_CELLS = 10_000_000
CHUNK_SIZE = 100_000
MIN_CELL_AREA = 0.1  # m^2

AUTHALIC_RADIUS = 6_371_007.180918473897976252  # m, ref: https://github.com/ecere/dggal/blob/7c496d4a8dff94821a38f33b4c37ad6abf459725/src/dggrs.ec#L29C24-L29C50
AUTHALIC_AREA = 4 * math.pi * AUTHALIC_RADIUS * AUTHALIC_RADIUS  # m^2
STANDARD_METERS_PER_PIXEL = 0.00028  # 0.28 mm/pixels -- following standard WMS 1.3.0 [OGC 06-042], SE and WMTS
METERS_PER_DEGREE = AUTHALIC_RADIUS * math.pi / 180  # m

ICOSA_EDGE_RADS = math.atan((2))
ICOSA_EDGE_M = ICOSA_EDGE_RADS * AUTHALIC_RADIUS

#######  Use for DGGS Binning
R_TRI = 0.88
R_QUAD = 0.7
R_PEN = 0.65
R_HEX = 0.62

#######  Use for DGGS Inspect
VMIN_HEX = 0.8
VMAX_HEX = 0.95
VCENTER_HEX = 0.9

VMIN_PEN = 0.7
VMAX_PEN = 0.9
VCENTER_PEN = 0.8

VMIN_QUAD = 0.6
VMAX_QUAD = 0.8
VCENTER_QUAD = 0.75

VMIN_TRI = 0.5
VMAX_TRI = 0.7
VCENTER_TRI = 0.65

FIX_ANTIMERIDIAN_CHOICES = [
    "shift",
    "shift_balanced",
    "shift_west",
    "shift_east",
    "split",
    "none",
]

OUTPUT_FORMATS = [
    None,
    "geojson_dict",
    "json_dict",
    "gpd",
    "geopandas",
    "gdf",
    "geodataframe",
    "csv",
    "geojson",
    "json",
    "shp",
    "shapefile",
    "gpkg",
    "geopackage",
    "parquet",
    "geoparquet",
]

STRUCTURED_FORMATS = [
    "gpd",
    "geopandas",
    "gdf",
    "geodataframe",
    "geojson_dict",
    "json_dict",
]

STATS_OPTIONS = [
    "count",
    "min",
    "max",
    "sum",
    "mean",
    "median",
    "std",
    "var",
    "range",
    "minority",
    "majority",
    "variety",
]

RASTER_STATS_OPTIONS = [
    "min",
    "max",
    "mean",
    "median",
    "sum",
    "count",
]

RASTER2DGGS_METHODS = [
    "binning",
    "nearest_neighbour",
]

# DGGAL types with their default resolution ranges (min_res, max_res) and class names.
# You can adjust these bounds manually as needed.
DGGAL_TYPES = {
    "gnosis": {
        "min_res": 0,
        "max_res": 28,
        "default_res": 16,
        "class_name": "GNOSISGlobalGrid",
    },
    "isea4r": {"min_res": 0, "max_res": 25, "default_res": 12, "class_name": "ISEA4R"},
    "isea9r": {"min_res": 0, "max_res": 16, "default_res": 10, "class_name": "ISEA9R"},
    "isea3h": {"min_res": 0, "max_res": 33, "default_res": 21, "class_name": "ISEA3H"},
    "isea7h": {"min_res": 0, "max_res": 19, "default_res": 11, "class_name": "ISEA7H"},
    "isea7h_z7": {
        "min_res": 0,
        "max_res": 19,
        "default_res": 11,
        "class_name": "ISEA7H_Z7",
    },
    "ivea4r": {"min_res": 0, "max_res": 25, "default_res": 12, "class_name": "IVEA4R"},
    "ivea9r": {"min_res": 0, "max_res": 16, "default_res": 10, "class_name": "IVEA9R"},
    "ivea3h": {"min_res": 0, "max_res": 33, "default_res": 21, "class_name": "IVEA3H"},
    "ivea7h": {"min_res": 0, "max_res": 19, "default_res": 11, "class_name": "IVEA7H"},
    "ivea7h_z7": {
        "min_res": 0,
        "max_res": 19,
        "default_res": 11,
        "class_name": "IVEA7H_Z7",
    },
    "rtea4r": {"min_res": 0, "max_res": 25, "default_res": 12, "class_name": "RTEA4R"},
    "rtea9r": {"min_res": 0, "max_res": 16, "default_res": 10, "class_name": "RTEA9R"},
    "rtea3h": {"min_res": 0, "max_res": 33, "default_res": 21, "class_name": "RTEA3H"},
    "rtea7h": {"min_res": 0, "max_res": 19, "default_res": 11, "class_name": "RTEA7H"},
    "rtea7h_z7": {
        "min_res": 0,
        "max_res": 19,
        "default_res": 11,
        "class_name": "RTEA7H_Z7",
    },
    "healpix": {
        "min_res": 0,
        "max_res": 26,
        "default_res": 18,
        "class_name": "HEALPix",
    },
    "rhealpix": {
        "min_res": 0,
        "max_res": 16,
        "default_res": 10,
        "class_name": "rHEALPix",
    },
}

DGGRID_TYPES = {
    "SUPERFUND": {
        "min_res": 0,
        "max_res": 17,
        "default_res": 9,
    },  # error when calling dggrid_instance.grid_stats_table
    "PLANETRISK": {"min_res": 0, "max_res": 22, "default_res": 13},
    "ISEA3H": {"min_res": 0, "max_res": 35, "default_res": 20},
    "ISEA4H": {"min_res": 0, "max_res": 29, "default_res": 16},
    "ISEA4T": {"min_res": 0, "max_res": 29, "default_res": 15},
    "ISEA4D": {"min_res": 0, "max_res": 29, "default_res": 16},
    "ISEA43H": {"min_res": 0, "max_res": 18, "default_res": 10},
    "ISEA7H": {"min_res": 0, "max_res": 21, "default_res": 11},
    "IGEO7": {
        "min_res": 0,
        "max_res": 20,
        "default_res": 12,
    },  # error when calling dggrid_instance.grid_stats_table
    "FULLER3H": {"min_res": 0, "max_res": 35, "default_res": 20},
    "FULLER4H": {"min_res": 0, "max_res": 30, "default_res": 16},
    "FULLER4T": {"min_res": 0, "max_res": 29, "default_res": 15},
    "FULLER4D": {"min_res": 0, "max_res": 30, "default_res": 16},
    "FULLER43H": {"min_res": 0, "max_res": 18, "default_res": 10},
    "FULLER7H": {"min_res": 0, "max_res": 21, "default_res": 11},
}

DGGS_TYPES = {
    "h3": {"min_res": 0, "max_res": 15, "default_res": 10},
    "s2": {"min_res": 0, "max_res": 30, "default_res": 16},
    "a5": {"min_res": 0, "max_res": 30, "default_res": 15},
    # "healpix": {"min_res": 0, "max_res": 29, "default_res": 10},  # to be checked
    "rhealpix": {"min_res": 0, "max_res": 15, "default_res": 10},
    "isea4t": {"min_res": 0, "max_res": 39, "default_res": 16},
    "isea3h": {"min_res": 0, "max_res": 40, "default_res": 20},
    "ease": {"min_res": 0, "max_res": 6, "default_res": 4},
    "qtm": {"min_res": 1, "max_res": 24, "default_res": 18},
    "olc": {"min_res": 2, "max_res": 15, "default_res": 10},
    "geohash": {"min_res": 1, "max_res": 12, "default_res": 7},
    "georef": {"min_res": 0, "max_res": 10, "default_res": 3},
    "mgrs": {"min_res": 0, "max_res": 5, "default_res": 3},
    "tilecode": {"min_res": 0, "max_res": 29, "default_res": 18},
    "quadkey": {"min_res": 0, "max_res": 29, "default_res": 18},
    "maidenhead": {"min_res": 1, "max_res": 4, "default_res": 4},
    "gars": {"min_res": 1, "max_res": 4, "default_res": 4},
    "digipin": {"min_res": 1, "max_res": 10, "default_res": 6},
}

DGGS_INSPECT = {
    "h3": {
        "min_res": 2,
        "max_res": 4,
    },
    "s2": {
        "min_res": 2,
        "max_res": 8,
    },
    "a5": {
        "min_res": 2,
        "max_res": 7,
    },
    "isea4t": {"min_res": 2, "max_res": 7},
    "rhealpix": {"min_res": 2, "max_res": 5},
    "dggrid_isea4t": {
        "min_res": 2,
        "max_res": 7,
    },  
    "dggrid_fuller4t": {
        "min_res": 2,
        "max_res": 7,
    },  
    "dggrid_isea4d": {
        "min_res": 2,
        "max_res": 7,
    },  
    "dggrid_fuller4d": {
        "min_res": 2,
        "max_res": 7,
    },  
    "dggrid_isea7h": {
        "min_res": 2,
        "max_res": 5,
    },  
    "dggrid_fuller7h": {
        "min_res": 2,
        "max_res": 5,
    },  
    "dggrid_igeo7": {
        "min_res": 2,
        "max_res": 5,
    },  
    "dggal_ivea3h": {"min_res": 2, "max_res": 9},
    "dggal_ivea4r": {"min_res": 2, "max_res": 7},
    "dggal_ivea7h": {"min_res": 2, "max_res": 5},
    "dggal_ivea9r": {"min_res": 2, "max_res": 5},
}


# ISEA4T Resolution to Accuracy mapping
ISEA4T_RES_ACCURACY_DICT = {
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
    26: 7.08 * 10**-3,  # accuracy returns 0.0, avg_edge_len =  0.11562
    27: 1.77 * 10**-3,  # accuracy returns 0.0, avg_edge_len =  0.05781
    28: 4.42 * 10**-4,  # accuracy returns 0.0, avg_edge_len =  0.0289
    29: 1.11 * 10**-4,  # accuracy returns 0.0, avg_edge_len =  0.01445
    30: 2.77 * 10**-5,  # accuracy returns 0.0, avg_edge_len = 0.00723
    31: 6.91 * 10**-6,  # accuracy returns 0.0, avg_edge_len =  0.00361
    32: 1.73 * 10**-6,  # accuracy returns 0.0, avg_edge_len =  0.00181
    33: 5.76 * 10**-7,  # accuracy returns 0.0, avg_edge_len = 0.0009
    34: 1.92 * 10**-7,  # accuracy returns 0.0, avg_edge_len = 0.00045
    35: 6.40 * 10**-8,  # accuracy returns 0.0, avg_edge_len = 0.00023
    36: 2.13 * 10**-8,  # accuracy returns 0, avg_edge_len = 0.00011
    37: 7.11 * 10**-9,  # accuracy returns 0.0, avg_edge_len = 6*10**(-5)
    38: 2.37 * 10**-9,  # accuracy returns 0.0, avg_edge_len = 3*10**(-5)
    39: 7.90 * 10**-10,  # accuracy returns 0.0, avg_edge_len = 10**(-5)
}

# ISEA3H Resolution to Accuracy mapping
ISEA3H_RES_ACCURACY_DICT = {
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
    40: 0.0000043155509996024,
}

# ISEA3H Accuracy to Resolution mapping
ISEA3H_ACCURACY_RES_DICT = {
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
    0.009438109955130478: 33,
    0.0031460366517101594: 34,
    0.0010486788839033865: 35,
    0.0003495596279677955: 36,
    0.0001165198769892652: 37,
    0.0000388399589964217: 38,
    0.0000129466529988072: 39,
    0.0000043155509996024: 40,
    # For resolutions 33-40, accuracy returns 0.0, so we map 0.0 to the highest resolution
    0.0: 40,
}
GEOREF_RESOLUTION_DEGREES = {
    0: 15.0,  # 15° x 15°
    1: 1.0,  # 1° x 1°
    2: 1 / 60,  # 1-minute
    3: 1 / 600,  # 0.1-minute
    4: 1 / 6000,  # 0.01-minute
    5: 1 / 60_000,  # 0.001-minute
    6: 1 / 600_000,  # 0.0001-minute
    7: 1 / 60_000_000,  # 0.00001-minute
    8: 1 / 600_000_000,  # 0.000001-minute
    9: 1 / 60_000_000_000,  # 0.0000001-minute
    10: 1 / 600_000_000_000,  # 0.00000001-minute
}

GARS_RESOLUTION_MINUTES = {
    1: 30,  # 30 minutes
    2: 15,  # 15 minutes
    3: 5,  # 5 minutes
    4: 1,  # 1 minute
}

INITIAL_GEOHASHES = [
    "b",
    "c",
    "f",
    "g",
    "u",
    "v",
    "y",
    "z",
    "8",
    "9",
    "d",
    "e",
    "s",
    "t",
    "w",
    "x",
    "0",
    "1",
    "2",
    "3",
    "p",
    "q",
    "r",
    "k",
    "m",
    "n",
    "h",
    "j",
    "4",
    "5",
    "6",
    "7",
]


ISEA4T_BASE_CELLS = [
    "00",
    "01",
    "02",
    "03",
    "04",
    "05",
    "06",
    "07",
    "08",
    "09",
    "10",
    "11",
    "12",
    "13",
    "14",
    "15",
    "16",
    "17",
    "18",
    "19",
]

ISEA3H_BASE_CELLS = [
    "00000,0",
    "01000,0",
    "02000,0",
    "03000,0",
    "04000,0",
    "05000,0",
    "06000,0",
    "07000,0",
    "08000,0",
    "09000,0",
    "10000,0",
    "11000,0",
    "12000,0",
    "13000,0",
    "14000,0",
    "15000,0",
    "16000,0",
    "17000,0",
    "18000,0",
    "19000,0",
]
