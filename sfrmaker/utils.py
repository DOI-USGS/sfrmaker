import os
import inspect
import json
import operator
import pprint
import time

import numpy as np
import pandas as pd
from shapely.geometry import Point

from sfrmaker.routing import get_nextupsegs, get_upsegs, make_graph

unit_conversion = {'feetmeters': 0.3048,
                   'metersfeet': 1 / .3048}


def assign_layers(reach_data, botm_array, pad=1., inplace=False):
    """Assigns the appropriate layer for each SFR reach,
            based on cell bottoms at location of reach.

    Parameters
    ----------
    reach_data : DataFrame with reach information
    botm : 3D numpy array of bottom elevations
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


def consolidate_reach_conductances(rd, keep_only_dominant=False):
    """For model cells with multiple SFR reaches, shift all conductance to widest reach,
    by adjusting the length, and setting the lengths in all smaller collocated reaches to 1,
    and the K-values in these to bedKmin

    Parameters
    ----------
    rd : DataFrame
        Reach data table.
    keep_only_dominant : bool
        Option to only retain the most downstream reach in each cell.
        (default False)

    Returns
    -------
    rd : DataFrame
        Modified version of rd, either culled to just one reach per cell,
        or with non-dominant reaches assigned streambed K values of zero.
    """
    print('\nAssigning total SFR conductance to dominant reach in cells with multiple reaches...')
    # use value of 1 for streambed k, since k would be constant within a cell
    rd['cond'] = rd.width * rd.rchlen * rd.strhc1  # assume value of 1 for strthick

    # make a new column that designates whether a reach is dominant in each cell
    # dominant reaches include those not collocated with other reaches, and the longest collocated reach
    rd['Dominant'] = [True] * len(rd)

    for c, nreaches in rd.node.value_counts().iteritems():
        # this is apparently nearly twice as fast as
        # a vectorized approach on all of m1.
        # apparently because it only operates on collocated reaches
        if nreaches > 1:
            # select the collocated reaches for this cell
            df = rd[rd.node == c].sort_values(by='width', ascending=False)

            # set all of these reaches except the largest to not Dominant
            rd.loc[df.index[1:], 'Dominant'] = False

    # Sum up the conductances for all of the collocated reaches
    # returns a series of conductance sums by model cell, put these into a new column in Mat1
    Cond_sums = rd[['node', 'cond']].groupby('node').agg('sum').cond
    rd['Cond_sum'] = [Cond_sums[c] for c in rd.node]

    # Calculate a new streambed Kv for widest reaches, set streambed Kv in secondary collocated reaches to 0
    m1d = rd.loc[rd.Dominant]
    rd.loc[rd.Dominant, 'strhc1'] = m1d.Cond_sum / (m1d.rchlen * m1d.width)
    rd.loc[~rd.Dominant, 'strhc1'] = 0.
    if keep_only_dominant:
        print('Dropping {} non-dominant reaches...'.format(np.sum(~rd.Dominant)))
        return rd.loc[rd.Dominant].copy()
    return rd


def interpolate_to_reaches(reach_data, segment_data,
                           segvar1, segvar2,
                           reach_data_group_col='iseg', segment_data_group_col='nseg'
                           ):
    """Interpolate values in datasets 6b and 6c to each reach in stream segment

    Parameters
    ----------
    segvar1 : str
        Column/variable name in segment_data array for representing start of segment
        (e.g. hcond1 for hydraulic conductivity)
        For segments with icalc=2 (specified channel geometry); if width1 is given,
        the eigth distance point (XCPT8) from dataset 6d will be used as the stream width.
        For icalc=3, an abitrary width of 5 is assigned.
        For icalc=4, the mean value for width given in item 6e is used.
    segvar2 : str
        Column/variable name in segment_data array for representing start of segment
        (e.g. hcond2 for hydraulic conductivity)

    Returns
    -------
    reach_values : 1D array
        One dimmensional array of interpolated values of same length as reach_data array.
        For example, hcond1 and hcond2 could be entered as inputs to get values for the
        strhc1 (hydraulic conductivity) column in reach_data.

    """
    if 'per' in segment_data.columns:
        segment_data = segment_data.loc[segment_data.per == 0]
    assert len(segment_data[segment_data_group_col].unique()) == len(segment_data), \
        "Segment ID column: {} has not non-unique values."
    rd_groups = reach_data.groupby(reach_data_group_col)
    sd_groups = segment_data.groupby(segment_data_group_col)

    reach_values = []
    for seg in segment_data[segment_data_group_col]:
        reaches = rd_groups.get_group(seg)
        segdata = sd_groups.get_group(seg)
        segment_length = reaches.rchlen.sum()
        # reach midpoint locations (to interpolate to)
        dist = (np.cumsum(reaches.rchlen) - 0.5 * reaches.rchlen).values
        fp = [segdata[segvar1].values[0],  # values at segment ends
              segdata[segvar2].values[0]]
        xp = [0, segment_length]  # segment start/end distances
        reach_values += np.interp(dist, xp, fp).tolist()
    assert len(reach_values) == len(reach_data)
    return np.array(reach_values)


def smooth_elevations(fromcomids, tocomids, elevup, elevdn):
    # make forward and reverse dictionaries with routing info
    graph = dict(zip(fromcomids, tocomids))
    assert 0 in set(graph.values()), 'No outlets in routing network!'
    graph_r = make_graph(tocomids, fromcomids)

    # make dictionaries of segment end elevations
    elevmin = dict(zip(fromcomids, elevdn))
    elevmax = dict(zip(fromcomids, elevup))

    def get_upseg_levels(seg):
        """Traverse routing network, returning a list of segments
        at each level upstream from the outlets. (level 0 route to seg;
        segments in level 1 route to a segment in level 0, etc.)

        Parameters:
        -----------
        seg : int
            Starting segment number

        Returns
        -------
        all_upsegs : list
            List with list of segments at each level
        """
        upsegs = graph_r[seg].copy()
        all_upsegs = [upsegs]
        for i in range(len(fromcomids)):
            upsegs = get_nextupsegs(graph_r, upsegs)
            if len(upsegs) > 0:
                all_upsegs.append(upsegs)
            else:
                break
        return all_upsegs

    def reset_elevations(seg):
        # reset segment elevations above (upsegs) and below (outseg) a node
        oseg = graph[seg]
        all_upsegs = np.array(list(get_upsegs(graph_r, seg)) + [seg])  # all segments upstream of node
        elevmin_s = np.min([elevmin[s] for s in all_upsegs])  # minimum current elevation upstream of node
        oldmin_s = elevmin[seg]
        elevs = [elevmin_s, oldmin_s]
        if oseg > 0:  # if segment is not an outlet,
            elevs.append(elevmax[oseg])  # outseg start elevation (already updated)
        # set segment end elevation as min of
        # upstream elevations, current elevation, outseg start elevation
        elevmin[seg] = np.min(elevs)
        # if the node is not an outlet, reset the outseg max if the current min is lower
        if oseg > 0:
            outseg_max = elevmax[oseg]
            elevmax[graph[seg]] = np.min([elevmin_s, outseg_max])

    print('\nSmoothing elevations...')
    ta = time.time()
    # get list of segments at each level, starting with 0 (outlet)
    segment_levels = get_upseg_levels(0)
    # at each level, reset all of the segment elevations as necessary
    for level in segment_levels:
        [reset_elevations(s) for s in level]
    print("finished in {:.2f}s".format(time.time() - ta))
    return elevmin, elevmax


def setup_reach_data(flowline_geoms, fl_comids, grid_intersections,
                     grid_geoms, tol=0.01):
    """Create prelimnary stream reaches from lists of grid cell intersections
    for each flowline.

    Parameters
    ----------
    flowline_geoms : list of shapely LineStrings
        LineStrings representing streams (e.g. read from a shapefile).
    fl_comids : list of integers
        List of unique ID numbers, one for each feature in
        flowline_geoms.
    grid_intersections : list of lists
        A list of intersecting grid cell node numbers for
        each feature in flowline_geoms and fl_comids.
    grid_geoms : list of shapely Polygons
        List of Polygons for each grid cell; each node number
        in grid_intersections refers to a polygon in this list.
    tol : float
        Maximum distance for identifying a connecting reach. Routing will
        not occur to reaches beyond this distance. This number should be small,
        because the ends of consecutive reaches should be touching if they were
        created via intersection with the model grid. (default 0.01)

    Returns
    -------
    m1 : DataFrame of reaches (one row per reach), with the following columns:
        ireach : int
            Reach number within a segment.
        iseg : int
            Segment number.
        node : int
            Grid cell number.
        line_id : int
            Original ID number for the intersected line.
        geometry : LineString
            LineString representing the intersected reach.
    """
    print("\nSetting up reach data... (may take a few minutes for large grids)")
    ta = time.time()
    fl_segments = np.arange(1, len(flowline_geoms) + 1)
    reach = []
    segment = []
    node = []
    geometry = []
    comids = []

    for i in range(len(flowline_geoms)):
        segment_geom = flowline_geoms[i]
        segment_nodes = grid_intersections[i]
        if segment_geom.type != 'MultiLineString' and segment_geom.type != 'GeometryCollection':
            ordered_reach_geoms, ordered_node_numbers = create_reaches(segment_geom, segment_nodes, grid_geoms, tol=tol)
            reach += list(np.arange(len(ordered_reach_geoms)) + 1)
            geometry += ordered_reach_geoms
            node += ordered_node_numbers
            segment += [fl_segments[i]] * len(ordered_reach_geoms)
            comids += [fl_comids[i]] * len(ordered_reach_geoms)
        else:
            start_reach = 0
            for j, part in enumerate(list(segment_geom.geoms)):
                geoms, node_numbers = create_reaches(part, segment_nodes, grid_geoms)
                if j > 0:
                    start_reach = reach[-1]
                reach += list(np.arange(start_reach, start_reach + len(geoms)) + 1)
                geometry += geoms
                node += node_numbers
                segment += [fl_segments[i]] * len(geoms)
                comids += [fl_comids[i]] * len(geoms)
        if len(reach) != len(segment):
            print('bad reach assignment!')
            break

    m1 = pd.DataFrame({'ireach': reach, 'iseg': segment, 'node': node,
                       'geometry': geometry, 'line_id': comids})
    m1.sort_values(by=['iseg', 'ireach'], inplace=True)
    m1['rno'] = np.arange(len(m1)) + 1
    print("finished in {:.2f}s\n".format(time.time() - ta))
    return m1


def create_reaches(part, segment_nodes, grid_geoms, tol=0.01):
    """Creates SFR reaches for a segment by ordering model cells
    intersected by a LineString part. Reaches within a part are
    ordered by proximity.

    Parameters
    ----------
    part: LineString
        Shapely LineString object (or a part of a MultiLineString)
    segment_nodes: list of integers
        Model node numbers intersected by *part*
    grid_geoms: list of Polygons
        List of shapely Polygon objects for the model grid cells, sorted by node number
    tol : float
        Maximum distance for identifying a connecting reach. Routing will
        not occur to reaches beyond this distance. This number should be small,
        because the ends of consecutive reaches should be touching if they were
        created via intersection with the model grid. (default 0.01)

    Returns
    -------
    ordered_reach_geoms: list of LineStrings
        List of LineString objects representing the SFR reaches for the segment.
    ordered_node_numbers: list of ints
        List of model cells containing the SFR reaches for the segment
    """
    reach_nodes = {}
    reach_geoms = {}
    # interesct flowline part with grid nodes
    reach_intersections = {c: part.intersection(grid_geoms[c]) for c in segment_nodes}
    reach_intersections = {c: g for c, g in reach_intersections.items() if
                           g.length > 0}  # drops points and empty geometries
    # empty geometries are created when segment_nodes variable includes nodes intersected by
    # other parts of a multipart line. Not sure what causes points besides duplicate vertices.

    # "flatten" all grid cell intersections to single part geometries
    n = 1
    for node, g in reach_intersections.items():
        if g.type == 'LineString':
            reach_nodes[n] = node
            reach_geoms[n] = g
            n += 1
        else:
            geoms = [gg for gg in g.geoms if gg.type == 'LineString']
            reach_nodes.update({n + nn: node for nn, gg in enumerate(geoms)})
            reach_geoms.update({n + nn: gg for nn, gg in enumerate(geoms)})
            n += len(geoms)

    # make point features for start and end of flowline part
    start = Point(part.coords[0])
    end = Point(part.coords[-1])

    ordered_reach_geoms = []
    ordered_node_numbers = []
    current_reach = start
    nreaches = len(reach_geoms)  # length before entries are popped

    # for each flowline part (reach)
    for i in range(nreaches):
        # find the next flowline part (reach) that touches the current reach
        dist = {j: g.distance(current_reach) for j, g in reach_geoms.items()}
        dist_sorted = sorted(dist.items(), key=operator.itemgetter(1))
        r = dist_sorted[0][0]
        next_reach = reach_geoms.pop(r)
        ordered_reach_geoms.append(next_reach)
        ordered_node_numbers.append(reach_nodes[r])
        current_reach = next_reach

        if current_reach.touches(end.buffer(tol)) and len(ordered_node_numbers) == nreaches:
            break
    assert len(ordered_node_numbers) == nreaches  # new list of ordered node numbers must include all flowline parts
    return ordered_reach_geoms, ordered_node_numbers


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


def load_json(jsonfile):
    """Convenience function to load a json file; replacing
    some escaped characters."""
    with open(jsonfile) as f:
        return json.load(f)


def load_sr(jsonfile):
    """Create a SpatialReference instance from model config json file."""
    from flopy.utils import SpatialReference
    with open(jsonfile) as input:
        cfg = json.load(input)

    rename = {'xoff': 'xll',
              'yoff': 'yll',
              }
    for k, v in rename.items():
        if k in cfg:
            cfg[v] = cfg.pop(k)
    cfg['delr'] = np.ones(cfg['ncol']) * cfg['delr']
    cfg['delc'] = np.ones(cfg['nrow']) * cfg['delc']
    kwargs = get_input_arguments(cfg, SpatialReference)
    return SpatialReference(**kwargs)


def get_input_arguments(kwargs, function, warn=True):
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
    print('\narguments to {}:'.format(function.__qualname__))
    params = inspect.signature(function)
    input_kwargs = {}
    not_arguments = {}
    for k, v in kwargs.items():
        if k in params.parameters:
            input_kwargs[k] = v
            print_item(k, v)
        else:
            not_arguments[k] = v
    if warn:
        print('\nother arguments:')
        for k, v in not_arguments.items():
            # print('{}: {}'.format(k, v))
            print_item(k, v)
    print('\n')
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