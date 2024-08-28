from collections.abc import Mapping
import os
import inspect
from pathlib import Path
import pprint
import warnings

import numpy as np
import flopy
import sfrmaker
from sfrmaker.routing import get_upsegs, make_graph
from sfrmaker.units import convert_length_units

unit_conversion = {'feetmeters': 0.3048,
                   'metersfeet': 1 / .3048}


def assign_layers(reach_data, botm_array, 
                  strtop_col='strtop', strthick_col='strthick',
                  idomain=None, pad=1., 
                  inplace=False):
    """Assigns the appropriate layer for each SFR reach,
            based on cell bottoms at location of reach.

    Parameters
    ----------
    reach_data : DataFrame
        Table of reach information, similar to SFRData.reach_data
    botm_array : ndarary
        3D numpy array of layer bottom elevations
    strtop_col : str
        Column name of streambed top elevations in reach_data.
        by default, 'strtop'.
    strthick_col : str
        Column name of streambed bottom thickness in reach_data.
        by default, 'strthick'
    idomain : ndarray  (optional)
        3D integer array of MODFLOW ibound or idomain values. Values > 0
        are considered active. Reaches in cells with values < 1 will be moved
        to the highest active cell if possible.
        by default, None (all cells are assumed to be active)
    pad : scalar
        Minimum distance that streambed bottom must be above layer bottom.
        When determining the layer or whether the streambed bottom is below
        the model bottom, streambed bottom - pad is used. Similarly, when
        lowering the model bottom to accomodate the streambed bottom,
        a value of streambed bottom - pad is used.
    inplace : bool
        If True, operate on reach_data and botm_array, otherwise,
        return 1D array of layer numbers and 2D array of new model botm elevations

    Returns
    -------
    (if inplace=True)
    layers : 1D array of layer numbers
    new_model_botms : ndarray or None
        New model bottom elevations. If no cell bottoms require adjustment, None
        is returned. If one or more cells require adjustment and ``idomain=None`` (not supplied), 
        a 2D array is returned, for backward compatibility with previous versions of SFRmaker.
        Otherwise, a 3D array of new bottom elevations is returned. In the future,
        a 3D array will always be returned.

    Notes
    -----
    Streambed bottom = strtop - strthick
    When multiple reaches occur in a cell, the lowest streambed bottom is used
    in determining the layer and any corrections to the model bottom.

    """
    i, j = reach_data.i.values, reach_data.j.values
    streambed_botm = reach_data[strtop_col].values - reach_data[strthick_col].values
    layers = get_layer(botm_array, i, j, streambed_botm - pad)
    if idomain is not None:
        active_frac = np.sum(idomain[layers, i, j] > 0)/len(i)
        if active_frac < 0.99:
            warnings.warn(f'{1-active_frac:.2%}of SFR reaches are in inactive cells. '
                          'Check the model idomain array and possibly the SFR '
                          'elevation and model length units for consistency.')
        

    # check against model bottom
    model_botm = botm_array[-1, i, j]
    # only adjust bottoms in i, j locations with at least one active cell
    in_active_area = np.array([True] * len(model_botm), dtype=bool)
    return_3d = True
    if idomain is None:
        return_3d = False
        idomain = np.ones_like(botm_array, dtype=int)
    else:
        idomain = (idomain > 0).astype(int)
    in_active_area = (idomain.sum(axis=0) > 0)[i, j]
    below = (streambed_botm - pad <= model_botm) & in_active_area
    below_inds = np.where(below)[0]
    nlay = idomain.shape[0]
    new_botm_array = None
    if np.any(below):
        new_botm_array = botm_array.copy()
        for pos in below_inds:
        #for n, (ib, jb) in enumerate(zip(below_i, below_j)):
            inds = (reach_data.i == i[pos]) & (
                    reach_data.j == j[pos])
            # reset the bottom of the lowest active layer
            # to accomodate the streambed bottom
            lowest_active_layer = nlay - np.argmax(idomain[:, i[pos], j[pos]][::-1]) - 1
            new_botm = streambed_botm[inds].min() - pad
            new_botm_array[lowest_active_layer:, i[pos], j[pos]] = new_botm
            layers[pos] = lowest_active_layer
        # check only reaches in active i, j locations
        assert not np.any(streambed_botm[in_active_area] <= \
            new_botm_array[-1, i, j][in_active_area])
    if inplace:
        if new_botm_array is not None:
            botm_array = new_botm_array
        reach_data['k'] = layers
    elif return_3d or new_botm_array is None:
        return layers, new_botm_array
    else:
        warnings.warn('Previously, assign_layers has returned a 2D (model bottom) array. '
                      'Now, if idomain is specified, a 3D array of model cell bottoms is returned. '
                      'In the future, a 3D array of model cell bottoms will always be returned.',
                      FutureWarning)
        return layers, new_botm_array[-1]


def get_layer(botm_array, i, j, elev):
    """Return the layers for elevations at i, j locations.

    Parameters
    ----------
    botm_array : 3D numpy array of layer bottom elevations
    i : scaler or sequence
        row index (zero-based)
    j : scaler or sequence
        column index
    elev : scaler or sequence
        elevation (in same units as model)

    Returns
    -------
    k : np.ndarray (1-D) or scalar
        zero-based layer index
    """

    def to_array(arg):
        if not isinstance(arg, np.ndarray):
            return np.array([arg])
        else:
            return arg

    i = to_array(i)
    j = to_array(j)
    nlay = botm_array.shape[0]
    elev = to_array(elev)
    botms = botm_array[:, i, j].tolist()
    layers = np.sum(((botms - elev) > 0), axis=0)
    # force elevations below model bottom into bottom layer
    layers[layers > nlay - 1] = nlay - 1
    layers = np.atleast_1d(np.squeeze(layers))
    if len(layers) == 1:
        layers = layers[0]
    return layers


def get_sfr_package_format(sfr_package_file):
    format = 'mf2005'
    with open(sfr_package_file) as src:
        for line in src:
            if 'being options' in line.lower():
                format = 'mf6'
                break
            elif 'begin packagedata' in line.lower():
                format = 'mf6'
                break
    return format


def arbolate_sum(segment, lengths, routing, starting_asums=None):
    """Compute the total length of all tributaries upstream from
    segment, including that segment, using the supplied lengths and
    routing connections.

    Parameters
    ----------
    segment : int or list of ints
        Segment or node number that is also a key in the lengths
        and routing dictionaries.
    lengths : dict
        Dictionary of lengths keyed by segment (node) numbers,
        including those in segment.
    routing : dict
        Dictionary describing routing connections between
        segments (nodes); values represent downstream connections.
    starting_asums : dict
        Option to supply starting arbolate sum values for any of
        the segments. By default, None.

    Returns
    -------
    asum : float or dict
        Arbolate sums for each segment.
    """
    if np.isscalar(segment):
        segment = [segment]
    graph_r = make_graph(list(routing.values()), list(routing.keys()))

    asum = {}
    for s in segment:
        upsegs = get_upsegs(graph_r, s)
        lnupsegs = [lengths[us] for us in upsegs]
        upstream_starting_asums = [0.]
        segment_starting_asum = 0.
        if starting_asums is not None:
            upstream_starting_asums = [starting_asums.get(us, 0.) for us in upsegs]
            segment_starting_asum = starting_asums.get(s, 0.)
        asum[s] = np.sum(lnupsegs) + lengths[s] + np.sum(upstream_starting_asums) + segment_starting_asum
    return asum


def get_input_arguments(kwargs, function, warn=False):
    """Return subset of keyword arguments in kwargs dict
    that are valid parameters to a function or method.

    Parameters
    ----------
    kwargs : dict (parameter names, values)
    function : function of class method

    Returns
    -------
    input_kwargs : dict
    """
    np.set_printoptions(threshold=20)
    #print('\narguments to {}:'.format(function.__qualname__))
    params = inspect.signature(function)
    input_kwargs = {}
    not_arguments = {}
    for k, v in kwargs.items():
        if k in params.parameters:
            input_kwargs[k] = v
            #print_item(k, v)
        else:
            not_arguments[k] = v
    if warn:
        #print('\nother arguments:')
        for k, v in not_arguments.items():
            # print('{}: {}'.format(k, v))
            print_item(k, v)
    #print('\n')
    return input_kwargs


def print_item(k, v):
    print('{}: '.format(k), end='')
    if isinstance(v, dict):
        # print(json.dumps(v, indent=4))
        pprint.pprint(v)
    elif isinstance(v, list):
        pprint.pprint(v)
    else:
        print(v)


def which(program):
    """Check for existance of executable.
    https://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file


def exe_exists(exe_name):
    exe_path = which(exe_name)
    if exe_path is not None:
        return os.path.exists(exe_path) and \
               os.access(which(exe_path), os.X_OK)


def update(d, u):
    """Recursively update a dictionary of varying depth
    d with items from u.
    from: https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
    """
    for k, v in u.items():
        if isinstance(d, Mapping):
            if isinstance(v, Mapping):
                d[k] = update(d.get(k, {}), v)
            else:
                d[k] = v
        else:
            d = {k: v}
    return d


def width_from_arbolate_sum(asum, a=0.1193, b=0.5032, minimum_width=1., input_units='meters',
                            output_units='meters'):
    """Estimate stream width from arbolate sum, using a power law regression
    comparing measured channel widths to arbolate sum.
    (after Leaf, 2020 and Feinstein et al. 2010, Appendix 2, p 266.)

    .. math::
        width = unit\_conversion * a * {asum_{(meters)}}^b

    Parameters
    ----------
    asum: float or 1D array
        Arbolate sum in the input units.
    a : float
        Multiplier parameter. Literature values:
        Feinstein et al (2010; Lake MI Basin): 0.1193
        Leaf (2020; Mississippi Embayment): 0.0592
    b : float
        Exponent in power law relationship. Literature values:
        Feinstein et al (2010; Lake MI Basin): 0.5032
        Leaf (2020; Mississippi Embayment): 0.5127
    minimum_width : float
        Minimum width to be returned. By default, 1.
    input_units : str, any length unit; e.g. {'m', 'meters', 'km', etc.}
        Length unit of asum
    output_units : str, any length unit; e.g. {'m', 'meters', 'ft', etc.}
        Length unit of output width

    Returns
    -------
    width: float
        Estimated width in feet

    Notes
    -----
    The original relationship described by Feinstein et al (2010) was for arbolate sum in meters
    and output widths in feet. The :math:`u` values above reflect this unit conversion. Therefore, the
    :math:`unit\_conversion` parameter above includes conversion of the input to meters, and output
    from feet to the specified output units.

    NaN arbolate sums are filled with the specified ``minimum_width``.

    References
    ----------
    see :doc:`References Cited <../references>`

    Examples
    --------
    Original equation from Feinstein et al (2010), for arbolate sum of 1,000 km:
    >>> width = width_from_arbolate_sum(1000, 0.1193, 0.5032, input_units='kilometers', output_units='feet')
    >>> float(round(width, 2))
    124.69
    """
    if not np.isscalar(asum):
        asum = np.atleast_1d(np.squeeze(asum))
    input_unit_conversion = convert_length_units(input_units, 'meters')
    output_unit_conversion = convert_length_units('feet', output_units)
    w = output_unit_conversion * a * (asum * input_unit_conversion) ** b
    if not np.isscalar(asum):
        w[w < minimum_width] = float(minimum_width)
        w[np.isnan(w)] = float(minimum_width)
    elif np.isnan(w) or w < minimum_width:
        w = minimum_width
    else:
        pass
    return w


def get_method_args_values(method, exclude=None):
    """Get arguments for a function or method and their default values.

    Parameters
    ----------
    method : function or method

    Returns
    -------
    args : dict

    """
    if exclude is None:
        exclude = set()
    args = {}
    inspected = inspect.signature(method)
    for parname, info in inspected.parameters.items():
        if parname not in {'self', 'kwargs'} and parname not in exclude:
            default_value = info.default
            if default_value == inspect._empty:
                default_value = None
            args[parname] = default_value
    return args


def make_config_summary():

    # methods called by SFRData.from_yaml
    # <block name in yaml>: [methods called by block]
    methods = {
        'modelgrid': {
            'Shapefile option (arguments to sfrmaker.grid.StructuredGrid.from_shapefile)':
            sfrmaker.grid.StructuredGrid.from_shapefile,
            'Flopy modelgrid option (arguments to flopy.discretization.StructuredGrid)':
            flopy.discretization.StructuredGrid
        },
        'flowlines': {
            'Custom hydrography option (arguments to sfrmaker.Lines.from_shapefile)':
            sfrmaker.Lines.from_shapefile,
            'NHDPlus option (arguments to sfrmaker.Lines.from_nhdplus_v2)':
            sfrmaker.Lines.from_nhdplus_v2
        },
        'simulation': {
            'Simulation information\n(MODFLOW-6 only; arguments to flopy.mf6.MFSimulation.load)':
                flopy.mf6.MFSimulation.load
        },
        'model': {
            'MODFLOW-NWT or MODFLOW-2005 information\n(arguments to flopy.modflow.Modflow.load)':
                flopy.modflow.Modflow.load,
            'MODFLOW-6 information\n(arguments to flopy.mf6.MFSimulation.get_model)':
                flopy.mf6.MFSimulation.get_model
        },
        'dem': {
            ('Option to sample streambed top elevations from a DEM\n'
             '(arguments to sfrmaker.SFRData.set_streambed_top_elevations_from_dem)'):
                sfrmaker.SFRData.set_streambed_top_elevations_from_dem,
             'keyword arguments to sfrmaker.SFRData.sample_reach_elevations':
                 sfrmaker.SFRData.sample_reach_elevations
        },
        'inflows': {
            ('Option to add specified inflows\n'
             '(arguments to sfrmaker.SFRData.add_to_perioddata for MODFLOW-6)'):
                sfrmaker.SFRData.add_to_perioddata,
            '(arguments to sfrmaker.SFRData.add_to_segment_data for MODFLOW-2005)':
                sfrmaker.SFRData.add_to_segment_data
        },
        'observations': {
            ('Option to add SFR observation input\n'
             '(arguments to sfrmaker.SFRData.add_observations)'):
                sfrmaker.SFRData.add_observations
        },
        'to_riv': {
            ('Option to convert parts of the SFR network\n'
             'downstream of specified line_ids to the RIV package'):
                sfrmaker.SFRData.to_riv
        },
        'options': {
            ('SFR options (arguments to sfrmaker.Lines.to_sfr)'):
                sfrmaker.Lines.to_sfr,
            ('keyword arguments to sfrmaker.SFRData'):
                sfrmaker.SFRData,
        }
    }

    # deprecated arguments
    # or arguments that are redundant/invalid in the context of the from_yaml
    deprecated = {
        'dem': {'dem', 'dem_z_units'},
        'options': {'model', 'grid',
                    'reach_data', 'segment_data',
                    'package_name'},
        'to_riv': {'segments', 'rno'}
    }

    config_summary = {
        'package_version': 'mf6  # mf6, mfnwt or mf2005',
        'package_name': ('model  # (defaults to model or '
                         'simulation name if not specified)'),
        'output_path': '"."  # where the output will be written'
    }
    for blockname, methods_called in methods.items():
        config_summary[blockname] = {}
        written = set()
        for comment, method in methods_called.items():
            exclude = deprecated.get(blockname, set())
            exclude.update(written)
            args = get_method_args_values(method, exclude=exclude)
            config_summary[blockname][f'{comment}'] = args
            written.update(args.keys())

    # edits to config summary generated by introspection
    # use 'namefile' instead of 'f' for namefile
    k = 'MODFLOW-NWT or MODFLOW-2005 information\n(arguments to flopy.modflow.Modflow.load)'
    del config_summary['model'][k]['f']
    # add namefile arg to first pos in dictionary
    config_summary['model'][k] = {'namefile': 'None  # "f" argument in flopy',
                                  **config_summary['model'][k]}

    project_root_path = Path(__file__).parent.parent
    config_summary_filename = project_root_path / 'docs/source/config-summary.yaml'
    with open(config_summary_filename, 'w') as dest:
        dest.write("# Summary of SFRmaker configuration file options\n\n")
        for blockname, subblock in config_summary.items():

            if subblock is None:
                subblock = 'None'
            if isinstance(subblock, str):
                dest.write(f'{blockname}: {subblock}\n')
                continue
            dest.write(f'\n{blockname}:\n')
            for comment, args in subblock.items():
                comment_text = '\n  # '.join(comment.split('\n'))
                dest.write(f"  # {comment_text}:\n")
                for argname, default_value in args.items():
                    dest.write(f'  {argname}: {default_value}\n')
    j=2