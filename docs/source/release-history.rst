===============
Release History
===============

Version 0.8.0 (2021-04-19)
--------------------------
* Inflows: 
    * Allow multiple inflows to be specified along a single headwater to outlet path, via an optional *one_inflow_per_path argument* (default False; previously was hardcoded as True).
* Runoff and other stress period data (including inflows)
    * Add support for specifying runoff (or other stress period inputs), allowing for multiple flows to be specified to a line, via lines IDs that may not be in the SFR package, but which route to lines in the SFR package. Routing information is provided by the ``flowline_routing argument`` to :meth:`sfrmaker.flows.add_to_perioddata` or :meth:`sfrmaker.sfrdata.SFRData.add_to_perioddata`. Line IDs may be missing if the linework was culled (i.e. to higher order streams or streams with a minimum arbolate sum). Runoff generated within the catchments of these culled lines still needs to be routed to the first line in the stream network.
    * Add option to distribute specified flows evenly among reaches associated with a line.
    * Refactor :meth:`sfrmaker.sfrdata.SFRData.period_data` to be indexed by stress period and reach number, allowing incremental updating (e.g. via pandas.DataFrame.update()). Previously, each call to :meth:`sfrmaker.flows.add_to_perioddata` would reset :meth:`sfrmaker.sfrdata.SFRData.period_data`. These changes allow for the specification of runoff to the SFR package in addition to specified inflows.
* Observations: 
    * base unique observations on name and type; allowing multiple observation types (e.g. downstream-flow and stage) to be appended to the observations table via add_observations
* Bug fixes:
    * Fix issue with starting arbolate sums, that was causing artificially narrow estimated widths on the first segment of any streams originating from outside of the model, by computing starting arbolate sums for each reach from the ending asum minus the line length.

Version 0.7.1 (2021-01-29)
--------------------------
USGS software release associated with `Groundwater` publication

Version 0.7.0 (2021-01-15)
--------------------------
* major speed-up (and overhead reduction) to finding routing paths (by replacing recursion strategy with simple for loop)
* in preprocessing module, use 1st percentile elevations sampled from DEM to avoid outliers (bad pixels)
* in preprocessing module, add option to re-use output from zonal statistics
* bug fix: refactor calls to gisutils.df2shp to use crs instead of epsg, etc.

Version 0.6.2 (2020-11-12)
--------------------------
* write unconnected reaches to connectiondata, as required by MODFLOW-6 v6.2

Version 0.6.1 (2020-11-04)
--------------------------
* deprecate sfrmaker.gis.CRS class in favor of :class:`pyproj.crs.CRS`
* add :func:`sfrmaker.routing.get_previous_ids_in_subset` function that can find outlet locations if the specified line IDs for outlets aren't in a consolidated (one_reach_per_cell=True) sfr network
* some bug fixes to the ``add_outlets`` option
* some fixes to the input data for the MERAS example

Version 0.6.0 (2020-10-15)
--------------------------
* ``add_outlets`` argument to :meth:`sfrmaker.lines.Lines.to_sfr` to add outlet conditions (outseg=0) at the locations of specified line IDs
* add :mod:`sfrmaker.preprocesing` module for culling NHDPlus flowlines, handling divergences, incorporating widths from the North American River Width (NARWidth) database, and reproducible editing of flowlines.
* fix :func:`sfrmaker.utils.width_from_arbolate_sum`: minimum width wasn't being implemented
* small fix to observations module to treat observation names as strings, even if they are digits
* small fix to `meth`:`sfrmaker.sfrdata.SFRData.write_package` if no options are supplied, set default fileout and obs6 entries to same location as SFR package file (by just writing the file names; previously, the full path to the SFR package file was written)

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
