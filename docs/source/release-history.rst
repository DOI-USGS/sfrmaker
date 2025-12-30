===============
Release History
===============

Version 0.13.0 (2025-12-30)
---------------------------
New Features:
**************
* Add support for incorporating field measurements of streambed top elevation.

Bug Fixes:
*************
* fix issue with elevations.sample_reach_elevations where csvs read by geopandas were being mistaken for vector data
* refactor code to support numpy API change from in1d to isin
* fix issue with the routing.find_path function: if the starting ID == the end ID, return the end ID (which is the last ID in any other returned path), instead of raising an error.
* drop_NHDPlusIDs compared as strings in load_nhdplus_hr()
* fix issue with load_nhdplus_hr() where user-supplied list of NHDPlusIDs were not being cast to strings for comparison with other NHDPlusIDs

Version 0.12.1 (2025-05-09)
---------------------------
* Fix Flopy packagedata error when there are auxiliary columns in the packagedata.

Version 0.12.0 (2025-05-09)
---------------------------
New Features:
**************
* Add support for preprocessing NHDPlus HR flowlines and waterbodies to shapefiles (`preprocessing.preprocess_nhdplus_hr_flowlines()` and `preprocessing.preprocess_nhdplus_hr_waterbodies()` functions)
* Add support for readling multiple NHDPlus HR geodatabases (which can then be consolidated)
* Add demo for converty MODFLOW SFR2 to MODFLOW 6

Internal Refactor:
*******************
* Use strings universally for line and other feature identifiers, to avoid confusion with floats vs ints, and issues with long identifiers such as those in NHDPlus HR that can't be represented as integers in 32-bit contexts

Bug Fixes:
*************
* fix issue with intersections when there are invalid geometries, by employging shapely’s `make_valid()` method.
* fix issue with reading and writing shapefiles when there are multiple geometry types (e.g. Multipolygon and Polygon)
* fix issue with convering MODFLOW 2005 SFR2 packages to MODFLOW 6, where line_id auxiliary variable wasn't included in the MODFLOW 6 packagedata.


Version 0.11.3 (2024-10-18)
---------------------------
* fix issue where SFR reach slopes were not being updated when streambed top elevations were resampled from a DEM
* fix issue with utils.assign_layers where a 2D array of the model bottom elevations was being returned instead of a 3D array of all layer bottom elevations.

Version 0.11.2 (2024-01-15)
---------------------------
* fix issue with line_id dtypes, so that NHDPlus Hi-Res COMIDs (which require 64-bits as integers) work more robustly on Windows.
* fixes to support Flopy renaming 'rno' to 'ifno'
* fixes to SFRData.add_to_segment_data:
    * add one_inflow_per_path argument similar to add_to_perioddata function for MODFLOW 6 with default of False; previously variables such as inflows were limited to a single value per flowline routing path (in the interest of avoiding duplicates), but many applications may have multiple inflows along a path.
    * sum any multiple values so that there is only one value per SFR segment (required for SFR2 Package)
* fixes to SFRData.create_modflow_sfr2: if no model is specified but the SFR package is MODFLOW 2005-style and transient (for example, with inflows specified by stress period), add stress periods and model grid information to flopy model instance (that is created so that flopy can be used to write the SFR package).
* fixes to mf5to6.py::segment_data_to_period_data: use nseg column instead of dataframe index for segment numbers
* bump supported python versions to 3.11 and 3.12
* refactor project setup to use pyproject.toml

Version 0.11.1 (2023-08-14)
---------------------------
* Refactor the lines class to replace the `attr_length_units`` arg with separate `asum_units` and `width_units` args, and rename the `attr_height_units`` arg to `elevation_units`.
* replace references to `epsg` and `proj_str` input args with `crs`.
* add warning to :py:func:`~sfrmaker.utils.assign_layers` if more than 1% of reaches land in an inactive cell.

Version 0.11.0 (2023-07-25)
---------------------------
 * enhancements to ``preprocessing.swb_runoff_to_csv()``:
    * add option to limit runoff output to a polygon area
    * add option to include x, y coordinates in output in a specified CRS (were previously included by default in the same CRS as the NetCDF input)
    * add start\_ and end_datetime args to limit the timespan of runoff output
    * add more robust detection of coordinate reference system of NetCDF input
    * add additional plotting and GIS layers export for diagnostic/QC checking (done by default)
 * fixes to adapted to breaking changes with pandas 2.0 and revised flopy modelgrid/crs interface
 * bug fix in cull_flowlines() function in preprocessing.py

Version 0.10.1 (2023-02-06)
---------------------------
* fixes to adapt to breaking changes in numpy 1.24

Version 0.10.0 (2023-02-02)
---------------------------
* the Lines.df attribute is now a GeoDataFrame with an attached crs attribute representing the current coordinate reference system for the flowlines.
* more fixes related to passing CRS references from the constructor method for Lines and the model grid attribute, through to reprojection to the model grid CRS and creation of the SFR package.
* rename `filter` argument to methods that read shapefiles to `bbox_filter`, to avoid overwriting built-in `filter`.

Version 0.9.4 (2023-01-19)
--------------------------
* add ability drop NHDPlus High Resolution segments by ftype (e.g., 428) and NHDPlusID, similar to the "drop_fcodes" option in the ``from_nhdplus_hr()`` method.
* fix issue where the crs argument for pyproj CRS-style coordinate references wasn't being passed to the base grid class; get the crs from a Flopy-style grid if it exists.
* updates to the install instructions, concerning Mamba and IPython kernels

Version 0.9.3 (2022-12-19)
--------------------------
* fixes to ``SFRdata.assign_layers``:
    * only compare new cell bottoms to old cell bottoms if a new bottom array is returned
    * when using idomain/ibound to set reach layers, recast array to boolean array of inactive/active cells, otherwise the ``np.argmax`` approach for finding the lowest active layer might identify the wrong layer (if ibound or idomain has values > 1). 
    * bug in file naming for revised bottom arrays (all arrays were previously named using the bottom layer).
* fixes to ``Lines`` class: 
    * add check for last update to routing dictionary, to avoid recursion if a user is trying to manually update the routing column in the attached dataframe (df) using the dictionary.
    * refactor internal updating of dataframe routing column in Lines.to_sfr() to its own method
    * use correct expression of ibound/idomain ``> 0`` when making isfr array, to handle cases where values > 1 are used to denote active cells.
* fixes to ``preprocessing.py`` module:
    * add support for cases where ``cull_invalid=False``, ``cull_isolated=False`` and ``asum_thresh=0`` (i.e. where all flowlines within the model area are kept
    * ``preprocessing.preprocess_nhdplus``: handle cases where no main stem is classified at a divergence; include warning in log file
* revert documentation links to ``aleaf`` fork, until GitHub Actions are allowed in ``doi-usgs`` org.
* add demo of preprocessing module to documentation

Version 0.9.2 (2022-08-01)
--------------------------
* add support for python 3.10; drop support for 3.8
* Bug fixes:
  * fix issue with SFRData where MODFLOW 6 datasets weren't being updated with any changes made to SFRdata
  * fix issue with empty elevation lists when determining routing


Version 0.9.1 (2021-12-30)
--------------------------
* Added support for NHDPlus High Resolution
* Features added to preprocessing:
  * Add get_flowline_routing function to get consolidated routing table in CSV format, for all flowlines in a set of NHDPlus v2 drainage basins
  * Add swb_runoff_to_csv function to associate gridded runoff estimates from Soil Water Balance Code with NHDPlus COMIDs (catchments), so that it can be input to sfrmaker.flows.add_to_perioddata()
  * Add keep_comids argument to cull_flowlines to retain specified COMIDs regardless of whether they meet culling criteria (asum, etc.)
* Features added to utils.assign_layers:
  * add arguments for streambed top and bottom thickness column names; 
  * add idomain argument to only adjust layer bottoms at active model locations. 
  * In places where the lowest active layer is not the bottom layer, the lowest active layer bottom will be pushed down, and the (inactive) layers beneath will be given zero-thickness
* Add arguments to set min/max streambed slope
* Bug fixes:
  * flows.add_to_perioddata: when distributing flows to reaches, handle comids that are in reach_data but not in the input data
  * handle structured grids with missing cells
  * fix some issues with rotated grid support, including the rotation sign and how lines outside of the model are culled
  * fix issue with how arbolate sums are computed in lines.to_sfr
  * fix issue with days to years unit conversion
  * fix an issue with observations.locate_sites when reaches are specified by segment and reach
  * preprocessing.preprocess_nhdplus: fix issue with auto reprojection and filtering for NARWidth data
  * observations.add_observations: force observation name column to object dtype on read_csv, so that USGS site numbers with leading 0s are preserved, for example
  * SFRData.sample_reach_elevations: ignore nan elevations sampled from DEM (or from beyond DEM extent)
* Refactoring:
  * sfrdata.set_streambed_top_elevations_from_dem: deprecate dem_z_units arg; replaced with 'elevation_units'

Version 0.9.0 (2021-04-19)
--------------------------
This release is the same as 0.8.0 due to a mistake, and has been yanked from PyPI.

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
