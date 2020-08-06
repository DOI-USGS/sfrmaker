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


Disclaimer
----------

This software is preliminary or provisional and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software.