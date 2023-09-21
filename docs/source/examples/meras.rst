MERAS 3: Creating an SFR package from a configuration file with custom hydrography
-------------------------------------------------------------------------------------
This example illustrates the use of custom hydrography with the configuration file, as well the scalability of SFRmaker. An SFR package is generated on a regular 1 km grid (Clark and others, 2018) that spans the Mississippi Embayment, a former bay of the Gulf of Mexico that includes portions of Missouri, Tennessee, Arkansas, Mississippi and Louisiana. The goal is to update the stream network for the Mississippi Embayment Regional Aquifer System (MERAS 2) model (Haugh and others, 2020; Clark and others, 2013) to include a realistic representation of the thousands of mapped streams (Figure 1). NHDPlus hydrography were preprocessed to only include the lines intersecting the MERAS footprint that had a total upstream drainage (arbolate sum) of 20 kilometers or greater (figure 1). Attribute information from NHDPlus, including routing connections, elevations at the upstream and downstream ends of line arcs and channel widths estimated from arbolate sum values (e.g. ref: Feinstein and others (2010, p 266)) were joined to the culled flowlines, which were then saved to a shapefile.

.. figure:: meras_comparison.png
    :align: left
    :scale: 70 %
    :alt: alternate text

    Figure 1: Location of the MERAS model domain with streams represented in the MERAS 2 model (A); streams mapped in NHDPlus version 2 (B).
    

Many streams enter the Mississippi Embayment with appreciable flow, which must be accounted for in the SFR package to achieve a realistic mass balance. Inflows to the SFR network for this study were derived from a combination of streamgage records and estimates from a statistical model (ref: MS water res. conference abstract). Inflow values for each site and model stress period were saved in a comma-separated-variable (CSV) format.

The configuration file for the Mississippi Embayment example is shown below. The configuration file, associate python script, and input files needed to reproduce the example are available on the `SFRmaker GitHub page <https://github.com/aleaf/sfrmaker/examples/meras/>`_. The configuration file is specified in the `YAML format <yaml.org>`_, which maps ``key: value`` pairs similar to a Python dictionary. In the SFRmaker configuration file, keys (to the left of the colons) indicate variables or groups of variables; values to apply to those variables are listed to the right of the colons. 

meras_sfrmaker_config.yml:
##############################
.. literalinclude:: ../../../examples/meras/meras_sfrmaker_config.yml
    :language: yaml
    :linenos:

Information on the model grid is specified in the ``modelgrid:`` block of the configuration file. The ``xoffset`` and ``yoffset`` arguments are the location of the lower left corner of the model grid, in units of the CRS defined by the EPSG code (typically meters). Since this is a uniform grid, only scalar values are provided for the row and column spacing. Variable row and column spacings can be specified using the list notation in YAML, or more simply by loading a flopy model instance into a Python script. SFRmaker does not require a modflow model as input, but a model can be optionally be specified in the configuration file with the ``model:`` and ``simulation:`` blocks (the latter is only needed for MODFLOW-6 models).

An advantage of specifying a model is that sfmaker will automatically assign stream reaches to the correct layer, based on the layer top and bottom elevations, and account for inactive cells. If the model has correct georeference information in the namefile header as assigned by `FloPy <https://github.com/modflowpy/flopy>`_, the ``modelgrid:`` block is not needed. An optional active_area: allows geographic area for generation of the SFR package input (that may or may not coincide with the active area of the model) to be defined. Otherwise, SFR input will be generated for the extent of the model grid.

While the MERAS hydrography includes elevations, this may not always be the case, or more accurate elevations may be available from a digital elevation model (DEM). If a ``dem:`` block is specified and ``set_streambed_top_elevations_from_dem: True``, SFRmaker will sample the DEM to the stream reaches as described in the methods section. The default buffer distance is 100 in the units of the projected CRS, or another buffer distance can be specified with an optional ``buffer_distance:`` argument.

The ``inflows:`` block allows for specification inflows to the SFR network. Inflow time series must include model stress period information (SFRmaker does not do any resampling) and column names in the input csv file must be specified. Similarly, an observations: block allows observation site locations to be specified in the coordinates of the projected CRS specified for the model grid, or using identifiers corresponding to flowlines in the input hydrography. With this information information, SFRmaker will generate input to the Gage Package (e.g. Prudic, 2004) or Modflow-6 observation utility (Langevin and others, 2017). 

The ``options:`` block can include keyword arguments to :meth:`sfrmaker.lines.Lines.to_sfr`, or other options such as ``set_streambed_top_elevations_from_dem:``.

The above configuration file can be used to generate an SFR package with the following Python code:

make_sfr.py:
##############################
.. literalinclude:: ../../../examples/meras/make_sfr.py
    :language: python
    :linenos:

This will produce an sfr package for MODFLOW-6, csv table representations of the SFR input, and shapefiles for visualizing the SFR package. The resulting SFR package is shown in Figure 2.

.. figure:: meras_results.png
    :align: left
    :scale: 100 %
    :alt: alternate text

    Figure 2: MERAS 3 SFR package, as shown by the shapefiles output from SFRmaker. Visualization of routing connections is illustrated in the map inset, where a connection crossing several cells indicates an error in the input hydrography routing attributes.