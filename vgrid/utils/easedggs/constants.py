'''
Module for generating EASE-DGGS grid_specs and levels_specs.

© Regents of the University of Minnesota. All rights reserved.
This software is released under an Apache 2.0 license. Further details about the Apache 2.0 license are available in the license.txt file.
'''

import numpy as np

grid_spec = {'ease': {'min_x': -17367530.445161372,
                    'min_y': -7314540.830638504,
                    'max_x': 17367530.445161372,
                    'max_y': 7314540.830638504},
            'geo': {'min_x': -180.0,
                    'min_y': -85.04456640737216,
                    'max_x': 180.0,
                    'max_y': 85.04456640737216}
            }

levels_specs = {0: {"refine_ratio": 4,
                        "n_row": 406,
                        "n_col": 964,
                        "x_length": 36032.22084058376,
                        "y_length": 36032.22084058376},
                1: {"refine_ratio": 3,
                        "n_row": 1624,
                        "n_col": 3856,
                        "x_length": 9008.05521014594,
                        "y_length": 9008.05521014594},
                2: {"refine_ratio": 3,
                        "n_row": 4872,
                        "n_col": 11568,
                        "x_length": 3002.6850700486466,
                        "y_length": 3002.6850700486466},
                3: {"refine_ratio": 10,
                        "n_row": 14616,
                        "n_col": 34704,
                        "x_length": 1000.8950233495489,
                        "y_length": 1000.895023349549},
                4: {"refine_ratio": 10,
                        "n_row": 146160,
                        "n_col": 347040,
                        "x_length": 100.08950233495489,
                        "y_length": 100.0895023349549},
                5: {"refine_ratio": 10,
                        "n_row": 1461600,
                        "n_col": 3470400,
                        "x_length": 10.008950233495488,
                        "y_length": 10.00895023349549},
                6: {"refine_ratio": 10,
                        "n_row": 14616000,
                        "n_col": 34704000,
                        "x_length": 1.0008950233495488,
                        "y_length": 1.000895023349549}
                }

ease_crs = 6933
geo_crs = 4326
mult_fac = [1, 4, 12, 36, 360, 3600, 36000]
cell_scale_factors = np.array([1.0 / mf for mf in mult_fac])

# def gen_pixel_scale_factors(**levels_specs):
#     '''
#     Generates the pixel scaling factors for the grid, using the
# characteristics from grid specifications.

#     Parameters
#     ----------
#     levels_specs : Dictionary
#        Dictionary with the grid charactersic, including refine_ratio.

#     Returns
#     -------
#     Numpy array with scale factors on the grid index scale.
#     '''
#     pixel_scale_factors =[1.]
#     scale_accum = 1

#     for item in levels_specs.items():
#         f = 1 / item[1]['refine_ratio']
# #     print(f)
#         scale_accum *= f
#         pixel_scale_factors.append(scale_accum)

#     return np.array(pixel_scale_factors)
