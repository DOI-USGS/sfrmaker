The Tyler Forks watershed: Creating an SFR package from a configuration file using NHDPlus
--------------------------------------------------------------------------------------------------
This example shows how a configuration file can be used to build an SFR package from NHDPlus data, with a model grid specified using a MODFLOW model or a shapefile. The example is set in the Tyler Forks watershed in northern Wisconsin, where groundwater/surface interactions are the primary modeling interest.

In the first configuration file, the Name file and workspace for a MODFLOW-NWT model are specified. An active area denoting where the SFR package will be built is provided with a shapefile. NHDPlus data are provided as a file path to the root folder level for a drainage basin (for example, 04, The Great Lakes), assuming the files within that path are in the same structure as the `download from the NHDPlus website <https://nhdplus.com/NHDPlus/NHDPlusV2_data.php>`_. Finally, a dem is provided as a `a more accurate source of streambed elevations <notebooks/Streambed_elevation_demo.html>`_.

tf_sfrmaker_config.yml:
##############################
.. literalinclude:: ../../../examples/tylerforks/tf_sfrmaker_config.yml
    :language: yaml
    :linenos:

In the second configuration file, no model is specified, so a package version, name, output path and length units are specified. The model grid is specified from a shapefile that has attribute fields indicating the row, column location of each cell. NHPlus data a specified as individual files.

.. literalinclude:: ../../../examples/tylerforks/tf_sfrmaker_config2.yml
    :language: yaml
    :linenos:

Either of these configuration files can then be used with a python script similar to the following:

tylerforks/make_sfr.py:
##############################
.. literalinclude:: ../../../examples/tylerforks/make_sfr.py
    :language: python
    :linenos:

This will produce an sfr package for MODFLOW-NWT, csv table representations of the SFR input, and shapefiles for visualizing the SFR package.

Running the tylerforks model
##############################
The above script can be found in the `examples/tylerforks folder <https://github.com/aleaf/sfrmaker/tree/develop/examples/tylerforks>`_ of the SFRmaker repository. Assuming the script was run from that location, the resulting MODFLOW model can then be run using the MODFLOW executable packaged with SFRmaker (on Windows):

.. code-block:: doscon

    cd tylerforks
    ../../../bin/win/mfnwt.exe tf.nam

(or on OSX)

.. code-block:: bash

    cd tylerforks
    ../../../bin/mac/mfnwt tf.nam