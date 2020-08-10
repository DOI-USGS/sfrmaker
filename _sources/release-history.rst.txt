===============
Release History
===============

Version 0.5.0 (2020-08-10)
--------------------------
* added from_yaml method to construct an SFR package from a configuration file
* deprecated use of the Flopy SpatialReference object
* add option to write MF6 packagedata block to an external file
* add default writing of source hydrography line_ids to MF6 packagedata as an auxiliary variable
* use pyproj CRS module internally for more robust handling of coordinate reference systems
* add starting gage package unit number attribute to SFRData that can be set by the user

Version 0.4.0 (2020-4-25)
--------------------------
* add automated setup (``SFRdata.add_observations`` method) of sfr observation locations (gages for mf2005 or obs input for mf6), from (x, y) locations, line_ids in source hydrography, or at specified reach numbers
* some minor patches to the gis module to use the new CRS module in pyproj to parse epsg codes and length units
* added screening of inactive cellids (set gw cellid to None) in MODFLOW6 flopy SFR package instance creation
* added support for MODFLOW-2005 gage package setup
* added support for specified inflows in MODFLOW-2005 (add_to_segment_data method)
* added to_riv method to convert SFR segments to RIV package; DataPackage base class

Version 0.3.0 (2020-2-12)
--------------------------
* replace FloPy SpatialReference with support for FloPy modelgrid
* updated modflow-6 execs to v 6.1
* minor big fixes/improvements to get_upsegs function in routing module

Version 0.2.1 (2019-12-12)
--------------------------
* fixed bug that was causing cases with only one (intersected stream) segment to fail

Version 0.2.0 (2019-12-08)
--------------------------
* added support for MODFLOW6 observation setup
* added support for specified inflows in MODFLOW6
* added get_inflow_locations_from_parent_model to get specified stream inflow values from a parent model
* bug fixes related to MODFLOW6 support

Version 0.1 Initial release (2019-07-12)
----------------------------------------
* see prior GitHub commits for "prehistory" of the project dating back to ArcPy scripts in 2013
