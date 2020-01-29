SFRmaker
===
SFRmaker is a python package for automating construction of stream flow routing networks from hydrography data. Hydrography are input from a polyline shapefile and intersected with a structured grid defined using a Flopy `SpatialReference` instance. Attribute data are supplied via `.dbf` files (NHDPlus input option) or via specified fields in the hydrography shapefile. Line fragments representing intersections between the flowlines and model grid cells are converted to SFR reaches using the supplied attribute data. MODFLOW-NWT/2005 or MODFLOW-6 SFR package input can then be written, along with shapefiles for visualizing the SFR package dataset.


### Version 0.2
[![Build Status](https://travis-ci.com/aleaf/SFRmaker.svg?branch=master)](https://travis-ci.com/aleaf/SFRmaker)
[![Build status](https://ci.appveyor.com/api/projects/status/0jk596k6osooyx1p/branch/master?svg=true)](https://ci.appveyor.com/project/aleaf/sfrmaker/branch/master)
[![Coverage Status](https://codecov.io/github/aleaf/SFRmaker/coverage.svg?branch=master)](https://codecov.io/github/aleaf/SFRmaker/coverage.svg?branch=master)
[![PyPI version](https://badge.fury.io/py/sfrmaker.svg)](https://badge.fury.io/py/sfrmaker)


Getting Started
----------------------------------------------- 

```python
import flopy
import sfrmaker
```
#### create an instance of the Lines class from NHDPlus data 
* alternatively, **`Lines`** can also be created from a shapefile or dataframe containing LineString features representing streams
* see input data requirements below for more details

```python
lns = sfrmaker.Lines.from_nhdplus_v2(NHDFlowlines='NHDFlowlines.shp',  
                            			PlusFlowlineVAA='PlusFlowlineVAA.dbf',  
                            			PlusFlow='PlusFlow.dbf',  
                            			elevslope='elevslope.dbf',  
                            			filter='data/grid.shp')
```

#### create an instance of `Lines` from a hydrography shapefile
* when creating `Lines` from a shapefile or dataframe, attribute field or column names can be supplied in lieu of the NHDPlus attribute tables (.dbf files).


```python
lns = Lines.from_shapefile(flowlines_file,
                           id_column='COMID',
                           routing_column='tocomid',
                           width1_column='width1',
                           width2_column='width2',
                           up_elevation_column='elevupsmo',
                           dn_elevation_column='elevdnsmo',
                           name_column='GNIS_NAME',
                           attr_length_units='feet',
                           attr_height_units='feet')
```
                     
#### create a flopy `SpatialReference` instance defining the model grid

```python

sr = flopy.utils.SpatialReference(delr=np.ones(160)*250,
                                  delc=np.ones(112)*250,
                                  lenuni=1,
                                  xll=682688, yll=5139052, rotation=0,
                                  proj_str='epsg:26715')
```

#### alternatively, the model grid can be defined with a `sfrmaker.StructuredGrid` instance
* can be created from a shapefile or `SpatialReference`.  
* an `active_area` polygon defines the area within the grid where SFR will be populated
* See example scripts in `Examples/` for more details.

```python

grd = StructuredGrid.from_shapefile(shapefile='Examples/data/badriver/grid.shp',
                                    icol='i',
                                    jcol='j',
                                    active_area='Examples/data/badriver/active_area.shp'.format(data_dir)
                                    )
```

#### intersect the lines with the model grid
* results in an **`sfrdata`** class instance

```python
sfr = lns.to_sfr(sr=sr)
```
or  

```python
sfr = lns.to_sfr(grid=grd)
```

#### write a sfr package file

```python
sfr.write_package('model.sfr')
```
#### write a MODFLOW 6 SFR package file:

```python
sfr.write_package('model.sfr6', version='mf6')
```
#### write shapefiles for visualizing the SFR package
```python
sfr.export_cells('sfr_cells.shp')
sfr.export_outlets('sfr_outlets.shp')
sfr.export_transient_variable('flow', 'sfr_inlets.shp') # inflows to SFR network
sfr.export_lines('sfr_lines.shp')
sfr.export_routing('sfr_routing.shp')
```

Installation
-----------------------------------------------

**Python versions:**

SFRmaker requires **Python** 3.6 (or higher)

**Dependencies:**  
pyyaml  
numpy  
pandas  
fiona  
rasterio  
rasterstats  
shapely  
pyproj  
rtree    
flopy  

### [Instructions for installing python and the required package dependencies](https://github.com/aleaf/SFRmaker/blob/master/Installing_dependencies.md)

Before using sfrmaker, the conda environment has to be activated:

```
conda activate sfrmaker
```
### Install SFRmaker from PyPi
with the `sfrmaker` environment activated:  

```
pip install sfrmaker
```


### Clone and install SFRmaker
see the [python install instructions](https://github.com/aleaf/SFRmaker/blob/master/Installing_dependencies.md) for more details on cloning a package.

```
git clone https://github.com/aleaf/SFRmaker.git

```

#### Install SFRmaker to the site_packages folder
from the root folder for the package (that contains `setup.py`):
  
```
python setup.py install
```
#### Install SFRmaker in current location (to current python path)
Instead of copying the source code to the python `site_packages` folder, this option creates a link so that python uses the source code in-situ. This is the best option if you want to modify the source code.

from the root folder for the package (that contains `setup.py`):


```  
pip install -e .
```
### running the test suite
from the root folder for the package (that contains `setup.py`):


```  
py.test
```
(requires the pytest package)

Input data requirements
-----------------------------------------------


#### 1) Hydrography data
##### NHDPlus v2 hydrography datasets    
 * Available at <http://www.horizon-systems.com/NHDPlus/NHDPlusV2_data.php>
 * Archives needed/relevant files:
 	* **NHDPlusV21\_XX\_YY\_NHDSnapshot_**.7z**   
 		* NHDFcode.dbf  
 		* NHDFlowline.dbf, .prj, .shp, .shx  
 	* **NHDPlusV21\_XX\_YY\_NHDPlusAttributes\_**.7z**  
 		* elevslope.dbf  
		* PlusFlow.dbf  
		* PlusFlowlineVAA.dbf
	* If your model domain encompasses multiple drainage areas, each type of NHDPlus file (e.g. NHDFlowline.shp, PlusFlow.dbf, etc.) can be supplied as a list. e.g.   
	
		```python
		NHDFlowlines=['<path to drainage area 1>/NHDFlowlines.shp',
		              '<path to drainage area 2>/NHDFlowlines.shp'...
		              ]
		
		```
	

	**Notes:**  

	* XX is Drainage Area ID (e.g., GL for Great Lakes) and YY is the Vector Processing Unit (VPU; e.g. 04) in the  above (see NHDPlus 	website for details).  


##### Other hydrography   
Any Polyline shapefile can be supplied in lieu of NHDPlus, but it must have the following columns, as shown in the second example:  

**flowlines\_file**: path to shapefile  
**id\_column**: unique identifier for each polyline  
**routing\_column**: downstream connection (ID), 0 if none  
**width1\_column**: channel width at start of line, in `attr\_length\_units` (optional)  
**width2\_column**: channel width at end of line, in `attr_length_units` (optional)  
**up\_elevation\_column**: streambed elevation at start of line, in `attr_height_units `  
**dn\_elevation\_column**: streambed elevation at end of line, in `attr_height_units `  
**name\_column**: stream name (optional)  
**attr\_length\_units**: channel width units  
**attr\_height\_units**: streambed elevation units  



#### 2) Model grid information
is supplied by creating a 	[`flopy.utils.SpatialReference`](https://github.com/modflowpy/flopy/blob/develop/flopy/utils/reference.py) instance or via a shapefile, as shown in the examples.


Running the example script
-----------------------------------------------
from the Examples folder:

```
python make_sfr.py
```

This creates an sfr package and adds it to the model in `Examples/data/badriver/tylerforks`.
Shapefiles for visualization of the sfr package are written to `Examples/temp`.