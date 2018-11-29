SFRmaker
===
Python package for automating construction of stream flow routing networks from hydrography data.

example usage:  

```python
import flopy
import sfrmaker
```
#### create an instance of the lines class from NHDPlus data 
* alternatively, **`lines`** can also be created from a shapefile or dataframe containing LineString features representing streams

```python
lns = sfrmaker.lines.from_NHDPlus_v2(NHDFlowlines='NHDFlowlines.shp',  
                            			PlusFlowlineVAA='PlusFlowlineVAA.dbf',  
                            			PlusFlow='PlusFlow.dbf',  
                            			elevslope='elevslope.dbf',  
                            			filter='data/grid.shp')
```
#### create an instance of `lines` from a hydrography shapefile
* when creating `lines` from a shapefile or dataframe, attribute field or column names can be supplied in lieu of the NHDPlus attribute tables (.dbf files).


```python
lns = lines.from_shapefile(flowlines_file,
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
                                  proj4_str='+init=epsg:26715')
```

#### intersect the lines with the model grid
* results in an **`sfrdata`** class instance

```python
sfr = lns.to_sfr(sr=sr)
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