#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 16:59:47 2021
más info de icesat https://nsidc.org/data/icesat-2/

@author: Nico Donoso nicolasdonosocastro@gmail.com
basado en https://github.com/icesat2py/icepyx
para analizar el lago subglaciar CECs 
"""

import icepyx as ipx
import os
import shutil

short_name = 'ATL08'
#spatial_extent = [-70.278409, -33.272140, -70.065084, -33.041871] #Olivares ampliado
spatial_extent = [(-70.3,-33.1), (-70.0, -33.0), (-70.0, -33.3), (-70.3, -33.3), (-70.3,-33.1)] # polygon vertices (here equivalent to the bounding box, above)
"""
 spatial extent = a region of interest to search within. This can be entered as a bounding box, polygon vertex coordinate pairs, or a polygon geospatial file (currently shp, kml, and gpkg are supported).

    bounding box: Given in decimal degrees for the lower left longitude, lower left latitude, upper right longitude, and upper right latitude
    polygon vertices: Given as longitude, latitude coordinate pairs of decimal degrees with the last entry a repeat of the first.
    polygon file: A string containing the full file path and name.
"""

#date_range = ['2018-10-1','2019-10-01']
#date_range = ['2019-10-1','2020-10-01']
#date_range = ['2020-10-1','2021-9-01']
date_range = ['2020-10-1','2021-09-01']
region_a = ipx.Query(short_name, spatial_extent, date_range)

region_a.visualize_spatial_extent()
region_a.dataset_summary_info()
print(region_a.latest_version())

#build and view the parameters that will be submitted in our query
region_a.CMRparams

earthdata_uid = 'ndonoso'
email = 'nicolasdonosocastro@gmail.com'
print("password = Ni3141592654")
region_a.earthdata_login(earthdata_uid, email)

request_mode = 'async' #Prueba

region_a.order_granules(verbose=True,subset=False)
region_a.granules.orderIDs

#path = '/home/ndonso/ice_sat2/python_icesat2/olivares/data'
path = '/home/ndonso/ice_sat2/python_icesat2/olivares/data'
region_a.download_granules(path)
