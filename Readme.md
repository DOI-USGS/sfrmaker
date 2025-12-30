SFRmaker
===
SFRmaker is a python package for automating construction of stream flow routing networks from hydrography data. Hydrography are input from a polyline shapefile and intersected with a structured grid defined using a shapefile or a Flopy `StructuredGrid` instance. Attribute data are supplied via `.dbf` files (NHDPlus input option) or via specified fields in the hydrography shapefile. Line fragments representing intersections between the flowlines and model grid cells are converted to SFR reaches using the supplied attribute data. MODFLOW-NWT/2005 or MODFLOW-6 SFR package input can then be written, along with shapefiles for visualizing the SFR package dataset.


### Version 0.13

![Tests](https://github.com/doi-usgs/sfrmaker/workflows/Tests/badge.svg)
[![Coverage Status](https://codecov.io/github/doi-usgs/SFRmaker/coverage.svg?branch=develop)](https://codecov.io/github/doi-usgs/SFRmaker/coverage.svg?branch=develop)
[![PyPI version](https://badge.fury.io/py/sfrmaker.svg)](https://badge.fury.io/py/sfrmaker)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/sfrmaker.svg)](https://anaconda.org/conda-forge/sfrmaker)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

Getting Started
----------------------------------------------- 
See the [SFRmaker documentation](https://doi-usgs.github.io/sfrmaker/latest/index.html)


Installation
-----------------------------------------------
See the [Installation Instructions](https://doi-usgs.github.io/sfrmaker/latest/installation.html)

How to cite
--------------
###### Citations for SFRmaker

Leaf, A.T., Fienen, M.N. and Reeves, H.W. (2021), SFRmaker and Linesink-Maker: Rapid Construction of Streamflow Routing Networks from Hydrography Data. Groundwater, 59: 761-771. https://doi.org/10.1111/gwat.13095

Leaf, A.T., 2023, Automated construction of Streamflow-Routing networks for MODFLOW—Application in the Mississippi Embayment region: U.S. Geological Survey Scientific Investigations Report 2023–5051, 28 p., https://doi.org/10.3133/sir20235051.

###### Software/Code citation for SFRmaker (IP-122355):
Leaf, A.T., Fienen, M.N. and Reeves, H.W., 2021, SFRmaker version 0.7.1: U.S. Geological Survey Software Release, 29 Jan. 2021, [https://doi.org/10.5066/P9U2T031](https://doi.org/10.5066/P9U2T031)

Disclaimer
----------

This software is preliminary or provisional and is subject to revision. It is
being provided to meet the need for timely best science. The software has not
received final approval by the U.S. Geological Survey (USGS). No warranty,
expressed or implied, is made by the USGS or the U.S. Government as to the
functionality of the software and related material nor shall the fact of release
constitute any such warranty. The software is provided on the condition that
neither the USGS nor the U.S. Government shall be held liable for any damages
resulting from the authorized or unauthorized use of the software. It is the responsibility of the user to check the accuracy of the results.

Any use of trade, firm, or product names is for descriptive purposes only and does not imply endorsement by the U.S. Government.