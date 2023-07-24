Philosophy
==========

Motivation
----------
Increasingly pressing and often complex water resources challenges require groundwater models that more holistically consider the integrated water cycle, and that can be built and updated on faster timeframes. While more advanced boundary condition formulations such as the Streamflow Routing (SFR) Package for MODFLOW (Prudic et al. 2004; Niswonger and Prudic 2005; Langevin et al. 2017) have long been available, they are often underutilized, presumably due to the complexity of input requirements (Anderson et al. 2015). Preparation of SFR input requires mapping linear stream features to a finite-difference grid, and populating each resulting stream reach with attribute properties. This process can be arduous and error-prone, requiring many geoprocessing operations or even extensive hand-digitizing of features. The number and complexity of operations presents a fundamental challenge to `scientific reproducibility <https://doi.org/10.1126/science.1213847>`_ and `step-wise modeling <https://haitjema.com/stepwise.html>`_ (Haitjema 1995).


What SFRmaker does
------------------------------
SFRmaker translates hydrography ("flowlines" or LineString features representing streams) into input for the MODFLOW-2005 (Harbaugh 2005) or MODFLOW-6 (Langevin 2017) SFR package, using the open-source Python geospatial stack. NHDPlus data (McKay et al. 2012) can be read directly, or any hydrography with routing attribute information can be used. Additional auxillary data such as observation locations, specified inflows, and surface runoff can be supplied in general formats such as shapefiles, CSV or NetCDF files. Streambed elevations can be estimated from a DEM, by sampling and then monotonically smoothing minimum elevations in the vicinity of each line. Overly dense stream networks can be pruned, and some common issues resolved semi-automatically with a preprocessing module that produces an intermediate set of (grid-independent) flowlines for a project area.

In general, SFRmaker is designed for the basic use case of a structured model grid and rectangular stream channel, where stage is estimated as part of the SFR package solution (ICALC=1 in MODFLOW-2005 parlance, or the input available in MODFLOW-6). Other options, such as specifying stream depth, unsaturated flow parameters and lake connections, can be accomplished within a script by editing the ``reach_data`` and ``segment_data`` tables associated with an ``SFRData`` object (see Basic Usage section). SFRmaker produces an SFR package input file, along with optional input to the MODFLOW 6 observation utility or HYDMOD, and shapefiles for visualizing the SFR package.


What SFRmaker doesn't do
------------------------------
SFRmaker doesn't map streams. Input hydrography can be obtained from databases such as NHDPlus, developed via flow accumulation methods (for example, Gardner et al. 2018) or even digitized by hand in a GIS environment. Development of SFRmaker has been primarily driven by project needs. As a result, it does not provide comprehensive support for the features available in the MODFLOW SFR Package. Currently, only structured model grids (those with row, column indices) are supported. As noted above, advanced SFR input options including lake connections, unsaturated flow beneath streams and diversions must be specified manually by editing the reach or segment data tables. Non-rectangular channel geometries are not supported. Support for unstructured grids and other SFR features may be added in the future. Contributions are welcome through GitHub pull requests; see the Contributing page.


Read more
------------------------------
Further description of SFRmaker and the companion `Linesink-maker package <https://github.com/DOI-USGS/linesink-maker>`_ for `analytic element models <https://www.epa.gov/ceam/gflow-groundwater-flow-analytic-element-model>`_ can be found in `Leaf et al. 2021. <https://doi.org/10.1111/gwat.13095>`_

Cited :ref:`References`
_________________________