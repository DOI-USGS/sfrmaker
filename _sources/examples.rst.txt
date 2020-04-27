========
Examples
========


Make an SFR package from NHDPlus files and a flopy model grid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: ../../Examples/make_sfr.py
    :language: python
    :lines: 5-

Make an SFR package without a prexisting Flopy model, using a shapefile to specify the grid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Reaches will be assigned to layer 1 by default, and Flopy diagnostics that only involved the SFR package can still be run.

.. literalinclude:: ../../Examples/make_sfr_without_model.py
    :language: python
    :lines: 7-