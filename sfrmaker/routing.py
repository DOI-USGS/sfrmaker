import time

import numpy as np


def pick_toids(routing, elevations):
    """Reduce routing connections to one per ID (no divergences).
    Select the downstream ID based on elevation, or first position in
    downstream ID list.

    Parameters
    ----------
    routing : dict
        Dictionary of id ints (keys) and to_id lists or sets (values).
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
        if isinstance(v, set):
            v = list(v)
        if isinstance(v, list):
            elevs = [elevations.get(vv, 1e5) for vv in v]
            # check for empty list
            if len(elevs) == 0:
                routing2[k] = 0
            else:
                routing2[k] = v[np.argmin(elevs)]
        elif np.isscalar(v):
            routing2[k] = v
    print("finished in {:.2f}s\n".format(time.time() - ta))
    return routing2


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
    nextupsegs = set()
    for s in upsegs:
        nextupsegs.update(graph_r.get(s, {}))
        #nextupsegs.update(graph_r[s])
    return nextupsegs


def get_upsegs(graph_r, seg):
    """Get all segments upstream of seg as a single flat set,
    by performing a breadth-first search of the routing graph,
    going in the upstream direction.

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
        if set(upsegs) == {0}:
            break
        elif len(upsegs) > 0:
            all_upsegs.update(upsegs)
        else:
            break
    if 0 in all_upsegs:
        all_upsegs.remove(0)
    return all_upsegs


def find_path(graph, start, end='0', limit=None):
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
    limit : int
        Option to limit the length of the path returned.
        By default, None (path is traced to the end routing number).

    Returns
    -------
    path : list
        List of segment numbers along routing path.
    """
    if limit is None:
        limit = len(graph)
    path = [start]
    if str(start) == str(end):
        return path
    next_id = start
    for i in range(limit):
        next_id = graph[next_id]
        path.append(next_id)
        if str(next_id) == str(end):
            break
    return path


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
        Dictionary of lists or ints (tocomids) keyed by values
        in fromcomids.
    """
    # from collections import defaultdict
    # fromcomids = np.array(fromcomids).astype(int)
    ## convert tocomid values to ints regardless of (enclosing) dtypes
    # tocomids = [a.astype(int).tolist() for a in map(np.array, tocomids)]
    # tuples = zip(fromcomids, tocomids)
    # graph = defaultdict(list)
    # if one_to_many: # tocomids should all be lists (not ints)
    #    for fromcomid, tocomid in tuples:
    #        v = graph[fromcomid] + tocomid
    #        graph[fromcomid] = list(set(v))
    # else: # tocomids should all be ints
    #    for fromcomid, tocomid in tuples:
    #        graph[fromcomid] = [tocomid]
    #    graph121 = {}
    #    for k, v in graph.items():
    #        assert len(v) == 1, "one_to_many=False but node {} connects to {}".format(k, v)
    #        graph121[k] = v.pop()
    #    return graph121
    # return graph
    from collections import defaultdict
    fromcomids = np.array(list(fromcomids))
    scalar_tocomids = np.all([np.isscalar(v) for v in tocomids])
    if scalar_tocomids:
        tocomid_sets = [{v} for v in tocomids]
    else:
        tocomid_sets = [set(a) for a in tocomids]
    tuples = zip(fromcomids, tocomid_sets)
    graph = defaultdict(set)
    for fromcomid, tocomid in tuples:
        graph[fromcomid].update(set(tocomid))
    if not one_to_many:
        graph121 = {}
        for k, v in graph.items():
            assert len(v) == 1, "one_to_many=False but node {} connects to {}".format(k, v)
            graph121[k] = v.pop()
        return graph121
    return graph


def make_reverse_graph(graph):
    """Make a reverse routing graph from a forward routing
    graph of {fromcomid: tocomid} connections.

    Parameters
    ----------
    graph : dict
        {fromcomid: tocomid} connections

    Returns
    -------
    graph_r : dict
        {tocomid: {fromcomid1, fromcomid2,...}} connections. Values
        are sets because tocomids will often have multiple fromcomids
        (tributaries).

    Examples
    --------
    >>> make_reverse_graph({1:2, 2:4, 3:4})
    {2: {1}, 4: {2, 3}}

    """
    graph_r = {}
    for fromid, toids in graph.items():
        if np.isscalar(toids):
            toids = {toids}
        for toid in toids:
            if toid not in graph_r:
                graph_r[toid] = {fromid}
            else:
                graph_r[toid].add(fromid)
    return graph_r


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
            r[u] = nexts if u > 0 else u  # handle lakes
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
        nseg, outseg = np.array(list(d.keys())), np.array(list(d.values()))
    ns = len(nseg)

    nexts = ns
    nextupsegs = nseg[outseg == 0]
    for i in range(ns):
        r, nexts, nextupsegs = reassign_upsegs(r, nexts, nextupsegs)
        if len(nextupsegs) == 0:
            break
    return r


def get_next_id_in_subset(subset, routing, ids):
    """If source linework are consolidated in the creation of
    SFR reaches (e.g. with lines.to_sfr(one_reach_per_cell=True)),
    not all line_ids in the source hydrography will be associated
    with a reach in the SFR dataset. This method finds the next downstream
    source line that is referenced in the reach data table (line_id column).

    Parameters
    ----------
    subset : list of ids that is a subset of the ids in routing
    routing : dict
        of id: to_id connections
    ids : iterable
        List of ids that are in routing but may not be in subset

    Returns
    -------
    ids : revised list of first values downstream of the values in ids (determined by routing)
        that are also in subset.
    """
    subset = set(subset).union({'0'})
    routing = routing.copy()
    if np.isscalar(ids):
        ids = [ids]
    paths = [find_path(routing, i) for i in ids]
    new_ids = []
    for p in paths:
        for id in p:
            if id in subset:
                new_ids.append(id)
                break
    assert len(new_ids) == len(ids)
    return new_ids


def get_previous_ids_in_subset(subset, routing, ids):
    """If source linework are consolidated in the creation of
    SFR reaches (e.g. with lines.to_sfr(one_reach_per_cell=True)),
    not all line_ids in the source hydrography will be associated
    with a reach in the SFR dataset. This method finds the previous (upstream)
    source line(s) that are referenced in the reach data table (line_id column).

    Parameters
    ----------
    subset : list of ids that is a subset of the ids in routing
    routing : dict
        of id: to_id connections
    ids : iterable
        List of ids that are in routing but may not be in subset

    Returns
    -------
    ids : revised list of first values upstream of the values in ids (determined by routing)
        that are also in subset.
    """
    subset = set(subset).union({0})
    routing = routing.copy()
    if np.isscalar(ids):
        ids = [ids]
    else:
        ids = ids.copy()
    graph_r = make_reverse_graph(routing)

    new_ids = set()
    nextupsegs = ids
    for i in range(len(graph_r)):
        for upseg in nextupsegs:
            if upseg in subset:
                new_ids.add(upseg)
        if len(new_ids) == len(nextupsegs):
            break
        nextupsegs = set(nextupsegs).difference(new_ids)
        if len(nextupsegs) > 0:
            nextupsegs = get_nextupsegs(graph_r, nextupsegs)
    return new_ids


def route_lines_by_proximity(flowline_geometries, line_ids=None, distance_tol=100):
    """Route flowlines based on the proximity of their start and end-points.
    
    Parameters
    ----------
    flowline_geometries : sequence of shapely linestring geometries
    line_ids : (optional) sequence of line ID numbers
        If None, the one-based index position of the lines is used, 
        with zero values indicating an outlet.
    distance_tol : float
        Maximum distance from a line end to the start of the next line. Lines
        outside of this distance will not be routed to. 
        
    Returns
    -------
    flowline_routing : sequence of downstream line ID numbers
        The next downstream line for each line in line_ids.
    """
    line_id_dtype = str
    if line_ids is None:
        line_ids = np.arange(1, len(flowline_geometries) + 1, dtype=int)
        line_id_dtype = int
    else:
        line_ids = np.array(list(line_ids))
        if len(line_ids) != len(flowline_geometries):
            raise ValueError(f"{len(line_ids)} line_ids for {len(flowline_geometries)} flowlines!")
        if len(set(line_ids)) != len(line_ids):
            raise ValueError(f"{len(line_ids) - len(set(line_ids))} duplicate line_ids!")
        line_id_dtype = type(line_ids[0])
        outlet_id = line_id_dtype(0)
        if outlet_id in line_ids:
            raise ValueError(f"'0' is reserved for an outlet condition, not allowed as a line_id.")
    start_xy = np.array([line.coords[0] for line in flowline_geometries])
    end_xy = np.array([line.coords[-1] for line in flowline_geometries])
    
    flowline_routing = []
    for line_id, (end_x, end_y) in zip(line_ids, end_xy):
        other_line_start_xy = start_xy[line_ids != line_id]
        other_line_ids = line_ids[line_ids != line_id]
        dist_to_start_xys = np.sqrt((end_x - other_line_start_xy[:, 0])**2 +\
            (end_y - other_line_start_xy[:, 1])**2)
        if np.min(dist_to_start_xys) <= distance_tol:
            next_start_id = other_line_ids[np.argmin(dist_to_start_xys)]
        else:
            next_start_id = line_id_dtype(0)
        flowline_routing.append(next_start_id)
        
    # go thru the paths and fix any instances of circular routing
    graph = make_graph(line_ids, flowline_routing, one_to_many=False)
    paths = {fid: find_path(graph, fid) for fid in graph.keys()}
    for line_id, routing_path in paths.items():
        if routing_path.count(line_id) > 1:
            upstream_ids = routing_path[:routing_path[1:].index(line_id)+1]
            fix_from_id = upstream_ids[-1]
            lines_not_upstream = np.array([True if lid in upstream_ids else False for lid in line_ids])#~np.isin(line_ids, upstream_ids)
            other_line_start_xy = start_xy[lines_not_upstream]
            other_line_ids = line_ids[lines_not_upstream]
            dist_to_start_xys = np.sqrt((end_x - other_line_start_xy[:, 0])**2 +\
                (end_y - other_line_start_xy[:, 1])**2)
            if np.min(dist_to_start_xys) <= distance_tol:
                next_start_id = other_line_ids[np.argmin(dist_to_start_xys)]
            else:
                next_start_id = outlet_id
            graph[fix_from_id] = next_start_id
        new_path = find_path(graph, line_id)
        if new_path.count(line_id) > 1:
            raise ValueError("Circular routing")
    for from_line, to_line in graph.items():
        if to_line in graph.keys() and to_line not in {0, '0'}:
            start_x, start_y = start_xy[line_ids == to_line][0]
            end_x, end_y = end_xy[line_ids == from_line][0]
            routing_dist = np.sqrt((start_x - end_x)**2 + (start_y - end_y)**2)
            assert routing_dist <= distance_tol
    
    flowline_routing = [graph[line_id] for line_id in line_ids]
    return flowline_routing
