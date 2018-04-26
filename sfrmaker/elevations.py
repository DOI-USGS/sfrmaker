"""Methods related to sampling and smoothing elevations."""
import time
import numpy as np
from .utils import make_graph, get_nextupsegs, get_upsegs

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