=====
Usage
=====


.. code-block:: python

   import flopy
   import sfrmaker

Create an instance of the Lines class from NHDPlus data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* alternatively, **``Lines``** can also be created from a shapefile or dataframe containing LineString features representing streams
* see input data requirements below for more details

.. code-block:: python

    lns = sfrmaker.Lines.from_nhdplus_v2(NHDFlowlines='NHDFlowlines.shp',
                                            PlusFlowlineVAA='PlusFlowlineVAA.dbf',
                                            PlusFlow='PlusFlow.dbf',
                                            elevslope='elevslope.dbf',
                                            filter='data/grid.shp')

.. note:: If your model domain encompasses multiple drainage areas, each type of NHDPlus file (e.g. NHDFlowline.shp, PlusFlow.dbf, etc. can be supplied as a list. e.g.

.. code-block:: python

    NHDFlowlines=['<path to drainage area 1>/NHDFlowlines.shp',
                  '<path to drainage area 2>/NHDFlowlines.shp'...
                  ]

create an instance of ``Lines`` from a hydrography shapefile
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
When creating ``Lines`` from a shapefile or dataframe, attribute field or column names can be supplied in lieu of the NHDPlus attribute tables (.dbf files).

.. code-block:: python

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


create a flopy ``StructuredGrid`` instance defining the model grid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
delr and delc have to specified in meters (consistent with projected CRS)

.. code-block:: python

    grid = flopy.discretization.StructuredGrid(delr=m.dis.delr.array * .3048,  # cell spacing along a row
                                               delc=m.dis.delc.array * .3048,  # cell spacing along a column
                                               xoff=682688, yoff=5139052,  # lower left corner of model grid
                                               angrot=0,  # grid is unrotated
                                               proj4='epsg:26715'
                                               # projected coordinate system of model (UTM NAD27 zone 15 North)
                                               )


Alternatively, the model grid can be defined with a ``sfrmaker.StructuredGrid`` instance
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* can be created from a shapefile or ``flopy.discretization.StructuredGrid``
* an ``active_area`` polygon defines the area within the grid where SFR will be populated
* See example scripts in ``Examples/`` for more details.

.. code-block:: python

    grid = StructuredGrid.from_shapefile(shapefile='Examples/data/badriver/grid.shp',
                                        icol='i',
                                        jcol='j',
                                        active_area='Examples/data/badriver/active_area.shp'.format(data_dir)
                                        )


Intersect the lines with the model grid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* results in an **``sfrdata``** class instance

.. code-block:: python

    sfr = lns.to_sfr(grid=grid)



write a sfr package file
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    sfr.write_package('model.sfr')


Write a MODFLOW 6 SFR package file:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    sfr.write_package('model.sfr6', version='mf6')


Write shapefiles for visualizing the SFR package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    sfr.export_cells('sfr_cells.shp')
    sfr.export_outlets('sfr_outlets.shp')
    sfr.export_transient_variable('flow', 'sfr_inlets.shp') # inflows to SFR network
    sfr.export_lines('sfr_lines.shp')
    sfr.export_routing('sfr_routing.shp')

