'''
This file configures logging for the library

© Regents of the University of Minnesota. All rights reserved.
This software is released under an Apache 2.0 license. Further details about the Apache 2.0 license are available in the license.txt file.
'''

import logging
import os

if "DEBUG" in os.environ and os.getenv("DEBUG") == "True":
    # print(f'OS DEBUG exists and is {os.getenv("DEBUG")}')
    numeric_level = getattr(logging, "DEBUG", None)
else:
    numeric_level = getattr(logging, "WARNING", None)

if not isinstance(numeric_level, int):
    raise ValueError('Invalid log level: %s' % numeric_level)

logging.basicConfig(
    level = numeric_level,
    format = '%(levelname)s:%(asctime)s %(message)s',
    encoding='utf-8'
)

logger = logging.getLogger('easedggs')
