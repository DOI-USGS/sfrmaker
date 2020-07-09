SFRmaker
===
SFRmaker is a python package for automating construction of stream flow routing networks from hydrography data. Hydrography are input from a polyline shapefile and intersected with a structured grid defined using a shapefile or a Flopy `StructuredGrid` instance. Attribute data are supplied via `.dbf` files (NHDPlus input option) or via specified fields in the hydrography shapefile. Line fragments representing intersections between the flowlines and model grid cells are converted to SFR reaches using the supplied attribute data. MODFLOW-NWT/2005 or MODFLOW-6 SFR package input can then be written, along with shapefiles for visualizing the SFR package dataset.


### Version 0.4
[![Build Status](https://travis-ci.com/aleaf/SFRmaker.svg?branch=master)](https://travis-ci.com/aleaf/SFRmaker)
[![Build status](https://ci.appveyor.com/api/projects/status/0jk596k6osooyx1p/branch/master?svg=true)](https://ci.appveyor.com/project/aleaf/sfrmaker/branch/master)
[![Coverage Status](https://codecov.io/github/aleaf/SFRmaker/coverage.svg?branch=master)](https://codecov.io/github/aleaf/SFRmaker/coverage.svg?branch=master)
[![PyPI version](https://badge.fury.io/py/sfrmaker.svg)](https://badge.fury.io/py/sfrmaker)


Getting Started
----------------------------------------------- 
See the [SFRmaker documentation](https://aleaf.github.io/sfrmaker/index.html)


Installation
-----------------------------------------------
See the [Installation Instructions](https://aleaf.github.io/sfrmaker/installation.html)
