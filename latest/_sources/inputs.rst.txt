-----------------------
Input data requirements
-----------------------


1) Hydrography data
------------------------

`NHDPlus v2 hydrography datasets`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Archives needed/relevant files:

**NHDPlusV21\_XX\_YY\_NHDSnapshot_\*\.7z**

* NHDFcode.dbf
* NHDFlowline.dbf, .prj, .shp, .shx

**NHDPlusV21\_XX\_YY\_NHDPlusAttributes\_\*\.7z**

* elevslope.dbf
* PlusFlow.dbf
* PlusFlowlineVAA.dbf

Notes:
######
XX is Drainage Area ID (e.g., GL for Great Lakes) and YY is the Vector Processing Unit (VPU; e.g. 04) in the  above (see NHDPlus website for details).


Custom hydrography
^^^^^^^^^^^^^^^^^^^^^^^^^
Any Polyline shapefile can be supplied in lieu of NHDPlus, but it must have the following columns, as shown in the second example:

* **flowlines\_file**: path to shapefile
* **id\_column**: unique identifier for each polyline
* **routing\_column**: downstream connection (ID), 0 if none
* **width1\_column**: channel width at start of line, in ``attr\_length\_units`` (optional)
* **width2\_column**: channel width at end of line, in ``attr_length_units`` (optional)
* **up\_elevation\_column**: streambed elevation at start of line, in ``attr_height_units``
* **dn\_elevation\_column**: streambed elevation at end of line, in ``attr_height_units``
* **name\_column**: stream name (optional)
* **attr\_length\_units**: channel width units
* **attr\_height\_units**: streambed elevation units



2) Model grid information
--------------------------
is supplied by creating a 	```flopy.discretization.Grid`` instance`_, instance or `via a shapefile`_, as shown in the :doc:`examples <examples/index>`.

.. _flopy.discretization.Grid`` instance`: https://aleaf.github.io/sfrmaker/usage.html#create-a-flopy-structuredgrid-instance-defining-the-model-grid

.. _NHDPlus v2 hydrography datasets: http://www.horizon-systems.com/NHDPlus/NHDPlusV2_data.php
.. _via a shapefile: https://aleaf.github.io/sfrmaker/usage.html#alternatively-the-model-grid-can-be-defined-with-a-sfrmaker-structuredgrid-instance