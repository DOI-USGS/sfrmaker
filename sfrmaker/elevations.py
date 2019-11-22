"""Methods related to sampling and smoothing elevations."""
import time

import numpy as np

from sfrmaker.routing import get_nextupsegs, get_upsegs, make_graph


def smooth_elevations(fromids, toids, elevations):  # elevup, elevdn):
    # make forward and reverse dictionaries with routing info
    graph = dict(zip(fromids, toids))
    assert 0 in set(graph.values()), 'No outlets in routing network!'
    graph_r = make_graph(toids, fromids)

    # make dictionaries of segment end elevations
    elevations = dict(zip(fromids, elevations))

    # elevmax = dict(zip(fromids, elevup))

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
        for i in range(len(fromids)):
            upsegs = get_nextupsegs(graph_r, upsegs)
            if len(upsegs) > 0:
                all_upsegs.append(upsegs)
            else:
                break
        return all_upsegs

    def reset_elevations(seg):
        """Reset segment elevations above (upsegs) and below (outseg) a node.
        """
        oseg = graph[seg]
        all_upsegs = np.array(list(get_upsegs(graph_r, seg)) + [seg])  # all segments upstream of node
        elevmin_s = np.min([elevations[s] for s in all_upsegs])  # minimum current elevation upstream of node
        oldmin_s = elevations[seg]
        elevs = [elevmin_s, oldmin_s]
        if oseg > 0:  # if segment is not an outlet,
            pass  # elevs.append(elevmax[oseg])  # outseg start elevation (already updated)
        # set segment end elevation as min of
        # upstream elevations, current elevation, outseg start elevation
        elevations[seg] = np.min(elevs)
        # if the node is not an outlet, reset the outseg max if the current min is lower
        if oseg > 0:
            # outseg_max = elevmax[oseg]
            next_reach_elev = elevations[oseg]
            elevations[graph[seg]] = np.min([elevmin_s, next_reach_elev])

    print('\nSmoothing elevations...')
    ta = time.time()
    # get list of segments at each level, starting with 0 (outlet)
    segment_levels = get_upseg_levels(0)
    # at each level, reset all of the segment elevations as necessary
    for level in segment_levels:
        [reset_elevations(s) for s in level]
    print("finished in {:.2f}s".format(time.time() - ta))
    return elevations
