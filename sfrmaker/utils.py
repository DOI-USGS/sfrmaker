import time
import operator
import numpy as np
import pandas as pd
from shapely.geometry import Point

unit_conversion = {'ftm': 0.3048,
                   'mft': 1/.3048}

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
    rd['cond'] = rd.width * rd.length * rd.strhc1 / rd.strthick

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

    # Calculate a new length for widest reaches, set length in secondary collocated reaches to 1
    # also set the K values in the secondary cells to 0
    m1d = rd.loc[rd.Dominant]
    rd.loc[rd.Dominant, 'strhc1'] = m1d.Cond_sum * m1d.strthick / (m1d.length * m1d.width)
    rd.loc[~rd.Dominant, 'strhc1'] = 0.
    if keep_only_dominant:
        print('Dropping {} non-dominant reaches...'.format(np.sum(~rd.Dominant)))
        return rd.loc[rd.Dominant].copy()
    return rd

def pick_toids(routing, elevations):
    """Reduce routing connections to one per ID (no divergences).
    Select the downstream ID based on elevation, or first position in
    downstream ID list.

    Parameters
    ----------
    routing : dict
        Dictionary of id ints (keys) and to_id lists (values).
    elevations : dict
        Dictionary of starting elevations (values) for each id (key)

    Returns
    -------
    routing2 : dict
        Same is input routing dictionary, except values have been
        reduced from lists to integers identifying the downstream
        connection.
    """
    print('\nPicking routing connections at divergences...')
    ta = time.time()
    routing2 = {}
    for k, v in routing.items():
        if isinstance(v, list):
            elevs = [elevations.get(vv, 1e5) for vv in v]
            routing2[k] = v[np.argmin(elevs)]
        elif np.isscalar(v):
            routing2[k] = v
    print("finished in {:.2f}s\n".format(time.time() - ta))
    return routing2

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


def get_nextupsegs(graph_r, upsegs):
    """Get adjacent upsegs for a list of segments
    as a single flat list.

    Parameters
    ----------
    graph_r : dict
        Dictionary of upstream routing connections.
        (keys=segments, values=adjacent upstream segments)
    upsegs : list
        List of segments

    Returns
    -------
    nextupsegs : list
        Flat list of next segments upstream from upsegs
    """
    nextupsegs = []
    for s in upsegs:
        nextupsegs += graph_r[s]
    return nextupsegs


def get_upsegs(graph_r, seg):
    """Get all segments upstream of seg as a single flat set.

    Parameters
    ----------
    graph_r : dict
        Dictionary of upstream routing connections.
        (keys=segments, values=adjacent upstream segments)
    seg : int
        Segment number (must be in graph.keys())

    Returns
    -------
    all_upsegs : set
        Flat set of all segments upstream from seg.
    """
    upsegs = graph_r[seg].copy()
    all_upsegs = upsegs
    for i in range(len(graph_r)):
        upsegs = get_nextupsegs(graph_r, upsegs)
        if len(upsegs) > 0:
            all_upsegs.update(upsegs)
        else:
            break
    return all_upsegs


def find_path(graph, start, end=0, path=[]):
    """Get a path through the routing network,
    from a segment to an outlet.

    Parameters
    ----------
    graph : dict
        Dictionary of seg : outseg numbers
    start : int
        Starting segment
    end : int
        Ending segment (default 0)

    Returns
    -------
    path : list
        List of segment numbers along routing path.
    """
    path = path + [start]
    if start == end:
        return path
    if start not in graph:
        return None
    if not isinstance(graph[start], list):
        graph[start] = [graph[start]]
    for node in graph[start]:
        if node not in path:
            newpath = find_path(graph, node, end, path)
            if newpath: return newpath
    return None


def make_graph(fromcomids, tocomids, one_to_many=True):
    """Make a dictionary of routing connections
    from fromcomids to tocomids.

    Parameters
    ----------
    fromcomids : list or 1D array
        Sequence of from nodes. The same value can
        appear more than once to describe convergence.
    tocomids : list or 1D array
        Sequence of to nodes. The same value can
        appear more than once to describe divergence.
    one_to_many : bool
        If True, values returned in graph are sets containing
        the tocomids associated with a fromcomid. If False,
        values are ints; each fromcomid only has one to comid.
    Returns
    -------
    graph : defaultdict
        Dictionary of lists (tocomids) keyed by values
        in fromcomids.
    """
    from collections import defaultdict
    fromcomids = np.array(fromcomids).astype(int)
    # convert tocomids to ints regardless of (enclosing) dtypes
    #tocomids = [a.astype(int).tolist() for a in map(np.array, tocomids)]
    tocomids = np.array(tocomids).astype(int)
    tuples = zip(fromcomids, tocomids)
    graph = defaultdict(set)
    for fromcomid, tocomid in tuples:
        graph[fromcomid].update({tocomid})
    if not one_to_many:
        graph121 = {}
        for k, v in graph.items():
            assert len(v) == 1, "one_to_many=False but node {} connects to {}".format(k, v)
            graph121[k] = v.pop()
        return graph121
    return graph

def renumber_segments(nseg, outseg):
    """Renumber segments so that segment numbering is continuous, starts at 1, and always increases
        in the downstream direction. Experience suggests that this can substantially speed
        convergence for some models using the NWT solver.

    Parameters
    ----------
    nseg : 1-D array
        Array of segment numbers
    outseg : 1-D array
        Array of outsegs for segments in nseg.

    Returns
    -------
    r : dict
        Dictionary mapping old segment numbers (keys) to new segment numbers (values). r only
        contains entries for number that were remapped.
    """
    if not isinstance(nseg, np.ndarray):
        nseg = np.array(nseg)
    if not isinstance(outseg, np.ndarray):
        outseg = np.array(outseg)

    def reassign_upsegs(r, nexts, upsegs):
        nextupsegs = []
        for u in upsegs:
            r[u] = nexts if u > 0 else u # handle lakes
            nexts -= 1
            nextupsegs += list(nseg[outseg == u])
        return r, nexts, nextupsegs

    print('enforcing best segment numbering...')
    # enforce that all outsegs not listed in nseg are converted to 0
    # but leave lakes alone
    r = {0: 0}
    r.update({o: 0 for o in outseg if o > 0 and o not in nseg})
    outseg = np.array([o if o in nseg or o < 0 else 0 for o in outseg])

    # if reach data are supplied, segment/outseg pairs may be listed more than once
    if len(nseg) != len(np.unique(nseg)):
        d = dict(zip(nseg, outseg))
        nseg, outseg = np.array(d.keys()), np.array(d.values())
    ns = len(nseg)

    nexts = ns
    nextupsegs = nseg[outseg == 0]
    for i in range(ns):
        r, nexts, nextupsegs = reassign_upsegs(r, nexts, nextupsegs)
        if len(nextupsegs) == 0:
            break
    return r

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
    fl_segments = np.arange(1, len(flowline_geoms)+1)
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
                reach += list(np.arange(start_reach, start_reach+len(geoms)) + 1)
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
    reach_intersections = {c:g for c, g in reach_intersections.items() if g.length > 0} #drops points and empty geometries
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
            reach_nodes.update({n+nn: node for nn, gg in enumerate(geoms)})
            reach_geoms.update({n+nn: gg for nn, gg in enumerate(geoms)})
            n += len(geoms)

    # make point features for start and end of flowline part
    start = Point(part.coords[0])
    end = Point(part.coords[-1])

    ordered_reach_geoms = []
    ordered_node_numbers = []
    current_reach = start
    nreaches = len(reach_geoms) # length before entries are popped

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
    w = 0.3048 * 0.1193 * (arbolate_sum) ** 0.5032
    if isinstance(arbolate_sum, np.ndarray):
        w[np.isnan(w)] = float(minimum_width)
    elif np.isnan(w):
        w = minimum_width
    else:
        pass
    return w


