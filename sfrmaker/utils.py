import collections
import os
import inspect
import pprint

import numpy as np

from sfrmaker.routing import get_upsegs, make_graph

unit_conversion = {'feetmeters': 0.3048,
                   'metersfeet': 1 / .3048}


def assign_layers(reach_data, botm_array, idomain=None, pad=1., inplace=False):
    """Assigns the appropriate layer for each SFR reach,
            based on cell bottoms at location of reach.

    Parameters
    ----------
    reach_data : DataFrame
        Table of reach information, similar to SFRData.reach_data
    botm : ndarary
        3D numpy array of layer bottom elevations
    idomain : ndarray
        3D integer array of MODFLOW ibound or idomain values. Values >=1
        are considered active. Reaches in cells with values < 1 will be moved
        to the highest active cell if possible.
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
    new_model_botms : 2D array of new model bottom elevations

    Notes
    -----
    Streambed bottom = strtop - strthick
    When multiple reaches occur in a cell, the lowest streambed bottom is used
    in determining the layer and any corrections to the model bottom.

    """
    i, j = reach_data.i.values, reach_data.j.values
    streambotms = reach_data.strtop.values - reach_data.strthick.values
    layers = get_layer(botm_array, i, j, streambotms - pad)

    # check against model bottom
    model_botm = botm_array[-1, i, j]
    below = streambotms - pad <= model_botm
    below_i = i[below]
    below_j = j[below]
    new_model_botm = None
    if np.any(below):
        new_model_botm = botm_array[-1].copy()
        for ib, jb in zip(below_i, below_j):
            inds = (reach_data.i == ib) & (
                    reach_data.j == jb)
            new_model_botm[ib, jb] = streambotms[inds].min() - pad
        assert not np.any(streambotms <= new_model_botm[i, j])
    if inplace:
        if new_model_botm is not None:
            botm_array[-1] = new_model_botm
        reach_data['k'] = layers
    else:
        return layers, new_model_botm


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


def arbolate_sum(segment, lengths, routing):
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

    Returns
    -------
    asum : float or dict
        Arbolate sums for each segment.
    """
    scalar = False
    if np.isscalar(segment):
        scalar = True
        segment = [segment]
    graph_r = make_graph(list(routing.values()), list(routing.keys()))

    asum = {}
    for s in segment:
        upsegs = get_upsegs(graph_r, s)
        lnupsegs = [lengths[s] for s in upsegs]
        asum[s] = np.sum(lnupsegs) + lengths[s]
    return asum


def width_from_arbolate_sum(arbolate_sum, minimum_width=1):
    """Estimate stream width in feet from arbolate sum in meters, using relationship
    described by Feinstein et al (2010), Appendix 2, p 266.

    Parameters
    ----------
    arbolate: float
        Arbolate sum in meters.

    Returns
    -------
    width: float
        Estimated width in meters (original formula returned width in ft)
    """
    scalar = np.isscalar(arbolate_sum)
    if scalar:
        arbolate_sum = np.array([arbolate_sum])
    w = 0.3048 * 0.1193 * arbolate_sum ** 0.5032
    fill = np.isnan(w) | (w == 0)
    w[fill] = minimum_width
    if scalar:
        return w[0]
    return w


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
        if isinstance(d, collections.Mapping):
            if isinstance(v, collections.Mapping):
                d[k] = update(d.get(k, {}), v)
            else:
                d[k] = v
        else:
            d = {k: v}
    return d
