"""Methods related to sampling and smoothing elevations."""
import time

import numpy as np

from sfrmaker.routing import get_nextupsegs, get_upsegs, make_graph


def get_slopes(streambed_tops, reach_lengths, reach_numbers, outreach_numbers, 
               default_slope=0.001, minimum_slope=0.0001,
                maximum_slope=1.):
    """Compute slopes by reach using values in strtop (streambed top) and rchlen (reach length)
    columns of reach_data. The slope for a reach n is computed as strtop(n) - strtop(n+1) / rchlen(n).
    Slopes for outlet reaches are set equal to a default value (default_slope).
    Populates the slope column in reach_data.

    Parameters
    ----------
    streambed_top : sequence
        Streambed top elevations
    reach_lengths : sequence
        Reach lengths
    reach_numbers : sequence
        Unique identifiers for each reach.
    outreach_numbers : sequence
        Unique identifier of next downtream reach.
    default_slope : float
        Slope value applied to outlet reaches (where water leaves the model).
        Default value is 0.001
    minimum_slope : float
        Assigned to reaches with computed slopes less than this value.
        This ensures that the Manning's equation won't produce unreasonable values of stage
        (in other words, that stage is consistent with assumption that
        streamflow is primarily drive by the streambed gradient).
        Default value is 0.0001.
    maximum_slope : float
        Assigned to reaches with computed slopes more than this value.
        Default value is 1.
    """
    # cast everything to lists to avoid confusion with numpy vs. pandas indexers
    streambed_tops = list(streambed_tops)
    reach_lengths = list(reach_lengths)
    reach_numbers = list(reach_numbers)
    outreach_numbers = list(outreach_numbers)
    assert np.sum(outreach_numbers) > 0, \
        ("outreach_numbers appear to be invalid; make sure outreaches are popluated, "
         "for example by running SFRData.set_outreaches()")
    elev = dict(zip(reach_numbers, streambed_tops))
    dist = dict(zip(reach_numbers, reach_lengths))
    dnelev = {rno: elev[outreach_numbers[i]] if outreach_numbers[i] != 0
              else -9999 for i, rno in enumerate(reach_numbers)}
    slopes = np.array(
        [(elev[rno] - dnelev[rno]) / dist[rno] if dnelev[rno] != -9999 and dist[rno] > 0
            else default_slope for rno in reach_numbers])
    slopes[slopes < minimum_slope] = minimum_slope
    slopes[slopes > maximum_slope] = maximum_slope
    return slopes
        
        
def smooth_elevations(fromids, toids, elevations, start_elevations=None):  # elevup, elevdn):
    """

    Parameters
    ----------
    fromids : sequence of hashables
    toids : sequence of hashables
        Downstream connections of fromids
    elevations : sequence of floats
        Elevation for each edge (line) in a stream network, or if start_elevations
        are specified, the end elevation for each edge.
    start_elevations : sequence of floats, optional
        Start elevation for edge (line) in a stream network.
        By default, None.

    Returns
    -------
    Elevations : dict or tuple
        Dictionary of smoothed edge elevations,
        or smoothed end elevations, start elevations
    """
    # make forward and reverse dictionaries with routing info
    graph = dict(zip(fromids, toids))
    assert 0 in set(graph.values()), 'No outlets in routing network!'
    graph_r = make_graph(toids, fromids)

    # make dictionaries of segment end elevations
    elevations = dict(zip(fromids, elevations))
    if start_elevations is not None:
        elevmax = dict(zip(fromids, start_elevations))

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
            if start_elevations is not None:
                elevs.append(elevmax[oseg]) # outseg start elevation (already updated)
        # set segment end elevation as min of
        # upstream elevations, current elevation, outseg start elevation
        elevations[seg] = np.min(elevs)
        # if the node is not an outlet, reset the outseg max if the current min is lower
        if oseg > 0:
            if start_elevations is not None:
                next_reach_elev = elevmax[oseg]
                elevmax[graph[seg]] = np.min([elevmin_s, next_reach_elev])
            else:
                next_reach_elev = elevations[oseg]
                elevations[graph[seg]] = np.min([elevmin_s, next_reach_elev])

    print('\nSmoothing elevations...')
    ta = time.time()
    # get list of segments at each level, starting with 0 (outlet)
    segment_levels = get_upseg_levels(0)
    # at each level, reset all of the segment elevations as necessary
    for level in segment_levels:
        for s in level:
            if 0 in level:
                j=2
            reset_elevations(s)
    print("finished in {:.2f}s".format(time.time() - ta))
    if start_elevations is not None:
        return elevations, elevmax
    return elevations
