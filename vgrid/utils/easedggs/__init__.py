"""
Intiialize src directory

© Regents of the University of Minnesota. All rights reserved.
This software is released under an Apache 2.0 license. Further details about the Apache 2.0 license are available in the license.txt file.
"""

import pyproj

# enabling PROJ Content Delivery Network to ensure consistent results between PowerShell and Anaconda/Command Prompt
pyproj.network.set_network_enabled(active=True)
