Philosophy
==========

Motivation
----------
As the world grapples with increasingly pressing and often complex water resources challenges, demand has grown for groundwater models that more realistically incorporate surface water. Often, groundwater models are motivated by questions regarding groundwater/surface water interactions and the availability or quality of in-stream flows. The Streamflow Routing (SFR) Package for MODFLOW (Prudic, 2004) allows for more realistic simulation of groundwater/surface water interactions and in-stream flows but remains underutilized, presumably due to the complexity of input requirements (Anderson and others, 2015). Preparation of SFR input requires mapping linear stream features to a finite-difference grid, and populating each finite difference cell with attributes describing the properties of stream reaches within that cell. This process can be arduous and error-prone, requiring many geoprocessing operations or even extensive hand-digitizing of features. The number and complexity of operations presents a fundamental challenge to scientific reproducibility (e.g ref: Peng 2011 Science article) and step-wise modeling (Haitjema, 1995).


What SFRmaker does
------------------------------
SFRmaker translates hydrography data into input for the MODFLOW-2005 (Harbaugh, 2005) or MODFLOW-6 (Langevin, 2017) SFR package, using the open-source python geospatial stack. NHDPlus data (McKay and others, 2012) can be read directly, or any hydrography with routing information (in an attribute field) can be used. In general, SFRmaker is designed for the basic use case of a structured model grid and rectangular stream channel, where stage is estimated as part of the SFR package solution (ICALC=1 in MODFLOW-2005 parlance, or the input available in MODFLOW-6). Other options, such as specifying stream depth, unsaturated flow parameters and lake connections, can be accomplished by editing ``reach_data`` and ``segment_data`` tables associated with an ``SFRData`` object (see Basic Usage section). SFRmaker produces an SFR package input file as output, as well as shapefiles for visualizing the SFR package.


What SFRmaker doesn't do
------------------------------
Development of Sfrmaker has been primarily driven by project needs. As a result, it does not provide comprehensive support for the features available in the MODFLOW SFR Package. Currently, only structured model grids (those with row, column indices) are supported. As noted above, advanced SFR input options including lake connections, unsaturated flow beneath streams and diversions must be specified manually by editing the reach or segment data tables. Non-rectangular channel geometries are not supported. Support for unstructured grids and other SFR features may be added in the future. Contributions are welcome through GitHub pull requests; see the Contributing page.
