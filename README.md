SFRmaker
=
Python package to automate construction of stream flow routing networks from hydrography data.

example usage:  

```python
import flopy
import sfrmaker
```
#### create an instance of the lines class from NHDPlus data 
* alternatively, **`lines`** can also be create from a shapefile or dataframe containing LineString features representing streams

```python
lns = sfrmaker.lines.from_NHDPlus_v2(NHDFlowlines='NHDFlowlines.shp',  
                            			PlusFlowlineVAA='PlusFlowlineVAA.dbf',  
                            			PlusFlow='PlusFlow.dbf',  
                            			elevslope='elevslope.dbf',  
                            			filter='data/grid.shp')
```
#### create a flopy `SpatialReference` instance defining the model grid

```python
sr = flopy.utils.SpatialReference(delr=np.ones(160)*250,
                                  delc=np.ones(112)*250,
                                  lenuni=1,
                                  xll=682688, yll=5139052, rotation=0,
                                  proj4_str='+init=epsg:26715')
```

#### intersect the lines with the model grid
* results in an **`sfrdata`** class instance

```python
sfr = lns.to_sfr(sr=sr)
```

#### write the sfr package file

```python
sfr.write_package('model.sfr')
```
#### to write a MODFLOW 6 SFR package:

```python
sfr.write_package('model.sfr', version='mf6')
```