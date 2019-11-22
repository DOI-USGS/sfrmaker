"""Script to set up test data from larger files
"""
import json

import numpy as np
from GISio import shp2df, df2shp
from flopy.utils.reference import SpatialReference

# basic grid parameters
name = 'map_test'
epsg = 5070
xul, yul = 520487.3, 1194668.3
nrow, ncol = 20, 20
dxy = 5280 * .3048
buf = 1e4
bounds = xul - buf, \
         yul - dxy * nrow - buf, \
         xul + dxy * ncol + buf, \
         yul + buf

# make version of preprocessed flowlines filtered to bounding box
df = shp2df('/Users/aleaf/Documents/MAP/repos/sfr_output/preprocessed/flowlines_gt20km/flowlines_gt20km_edited.shp',
            filter=bounds)
df2shp(df, 'data/{}_flowlines.shp'.format(name), epsg=epsg)

# make a spatial reference object defining the grid
sr = SpatialReference(delr=np.ones(ncol, dtype=float) * dxy,
                      delc=np.ones(nrow, dtype=float) * dxy,
                      xul=xul, yul=yul, epsg=epsg)
# export sr info to json file
model_info = sr.attribute_dict
model_info['nrow'] = sr.nrow
model_info['ncol'] = sr.ncol
model_info['delr'] = sr.delr[0]
model_info['delc'] = sr.delc[0]
model_info['epsg'] = sr.epsg

with open('data/{}_grid.json'.format(name), 'w') as output:
    json.dump(model_info, output, indent=4, sort_keys=True)
