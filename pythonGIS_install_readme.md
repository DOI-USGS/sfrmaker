##Instructions for installing the python GIS tool chain on Windows

The family of packages described below add geospatial functionality to Python, via bindings to underlying open-source GIS libraries that are generally robust, fast, and available for free. These packages can provide a powerful alternative to ArcPy with several advantages, including speed, reduced dependency on input and output operations, and easier integration with python datastructions. However, because of their dependency on external libraries (mostly written in C or C++), they have historically been difficult to install on Windows.

Christoph Gohlke's Unofficial Windows Binaries site (<http://www.lfd.uci.edu/~gohlke/pythonlibs/>) has made windows installation a lot easier, by packaging the python modules together with their library dependencies into binary ("wheel") files that can be installed using Python's **pip** package management system.

A brief description of the packages and their dependencies:  

| python package| what it does | underlying library  | python dependencies | unofficial binary |  
|:--------------|:------------ | :------------------:| -------------------:| :-----------------|
| **fiona**     | reads and writes shapefiles | OGR | six, cligj |<http://www.lfd.uci.edu/~gohlke/pythonlibs/#fiona>  |
| **shapely**   | provides containers and operations for manipulating vector data (points, lines and polygons)|  GEOS | -- | <http://www.lfd.uci.edu/~gohlke/pythonlibs/#shapely> |
| **gdal** | python bindings for GDAL library |  GDAL | numpy |<http://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal> |
| **rasterio** | reads and writes rasters | GDAL | **gdal**, affine, cligj (and click), enum34, numpy | <http://www.lfd.uci.edu/~gohlke/pythonlibs/#rasterio>|
| **rtree**     | spatial indexing to speed up large intersection operations | libspatialindex | -- | <http://www.lfd.uci.edu/~gohlke/pythonlibs/#rtree>|
| **pyproj**     | cartographic transformations | PROJ.4 | -- |<http://www.lfd.uci.edu/~gohlke/pythonlibs/#pyproj>|
| **rasterstats**| summarizes raster data (e.g. zonal statistics)  | GDAL | **gdal**, **rasterio**, **fiona**, **shapely** and numpy| just install using pip (see below) |

###1) Installing Python
If you are not already a python user, I recommend installing the Anaconda Python Distribution (<https://store.continuum.io/cshop/anaconda/>), which is available for free and comes with many popular python packages pre-installed. Assuming you have a 64-bit system, I recommend installing the 64-bit version of Anaconda, as the 32-bit version has memory limitations.

###2) Installing pip
Anaconda comes with its own package management system, **conda**, which is especially good at managing dependencies between packages. However, many of the above packages are not available for Windows on **conda**. To minimize the potential for version conflicts, I recommend installing this whole GIS tool chain through the Unofficial Windows Binaries site, using **pip**. Pip can be installed using **conda**. After installing Anaconda, open a command window anywhere (e.g., by typing "cmd" in the Start Menu search box) and enter:  

```  
>conda install pip  
```
###3) Downloading the wheel files
Installer files (called "wheels") can be downloaded for each package at the links given in the above table (alternatively, one can go to the main page listed above and search for the package name). Select the appropriate file for your python version. For example, to install gdal version 1.11.3 on 64-bit Python 2.7, one would choose **GDAL‑1.11.3‑cp27‑none‑win_amd64.whl**.

###4) Installing the wheels using pip
Open a command window in your download folder and install each package using pip:  
e.g:

```
>pip install Fiona-1.6.2-cp27-none-win_amd64.whl
```
**Note:** **GDAL** must be installed prior to installing **rasterio**, and **gdal**, **rasterio**, **fiona** and **shapely** must be installed prior to installing **rasterstats**. After these other packages are installed, **rasterstats** can simply be installed from **PyPl** using pip.  
e.g:
(note that the actual *.whl filenames will vary depending on package and python versions):

```
>pip install Fiona-1.6.2-cp27-none-win_amd64.whl
>pip install Shapely-1.5.12-cp27-none-win_amd64.whl
>pip install GDAL-1.11.3-cp27-none-win_amd64.whl
>pip install rasterio-0.28.0-cp27-none-win_amd64.whl
>pip install rasterstats
```
###Common problems
* #####message **"requirement already satisfied"** when trying to install a package 
	**Cause:** some form of the package is already installed on your python distribution  
	**Solution:** try either uninstalling the existing package first

	```
>pip uninstall <packagename>
>pip install <wheel file>
```
or upgrading
  
	```
>pip install --upgrade <wheel file>
```  

* #####"side-by-side configuration" error message when trying to use rasterstats or rasterio 
	**Cause:** gdal may be installed incorrectly or there is a version conflict.  
	**Solution:** try ```pip uninstall gdal```, obtain a gdal binary wheel from <http://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal>, and install the wheel using pip.
