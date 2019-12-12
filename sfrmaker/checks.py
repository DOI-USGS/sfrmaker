import numpy as np
import pandas as pd

from sfrmaker.routing import find_path, make_graph


def valid_rnos(rnos):
    """Check that unique reach numbers (rno in MODFLOW 6)
    are consecutive and start at 1.
    """
    sorted_reaches = sorted(rnos)
    consecutive = np.diff(sorted_reaches).sum() \
                  == len(rnos) - 1
    onebased = np.min(sorted_reaches) == 1
    return consecutive & onebased


def valid_nsegs(nsegs, outsegs=None, increasing=True):
    """Check that segment numbers are valid.

    Parameters
    ----------
    nsegs : list of segment numbers
    outsegs : list of corresponding routing connections
        Required if increasing=True.
    increasing : bool
        If True, segment numbers must also only increase downstream.
    """
    # cast to array if list or series
    nsegs = np.atleast_1d(nsegs)
    outsegs = np.atleast_1d(outsegs)
    consecutive_and_onebased = valid_rnos(nsegs)
    if increasing:
        assert outsegs is not None
        graph = make_graph(nsegs, outsegs, one_to_many=False)
        monotonic = []
        for s in nsegs:
            seg_sequence = find_path(graph.copy(), s)[:-1]  # last number is 0 for outlet
            monotonic.append(np.all(np.diff(np.array(seg_sequence)) > 0))
        monotonic = np.all(monotonic)
        return consecutive_and_onebased & monotonic
    else:
        return consecutive_and_onebased


def is_to_one(toid_sequence):
    squeezed = np.squeeze(list(toid_sequence))
    if squeezed.ndim == 0:
        return True
    else:
        return np.isscalar(np.squeeze(list(toid_sequence))[0])


def rno_nseg_routing_consistent(nseg, outseg, iseg, ireach, rno, outreach):
    """Check that routing of segments (MODFLOW-2005 style) is consistent
    with routing between unique reach numbers (rno; MODFLOW 6 style)

    Parameters
    ----------
    nseg : list or 1D numpy array
    outseg : list or 1D numpy array
    iseg : list or 1D numpy array
    ireach : list or 1D numpy array
    rno : list or 1D numpy array
    outreach : list or 1D numpy array

    Returns
    -------
    consistent : bool
    """
    df = pd.DataFrame({'nseg': nseg,
                       'outseg': outseg,
                       'iseg': iseg,
                       'ireach': ireach,
                       'rno': rno,
                       'outreach': outreach})
    df.sort_values(by=['iseg', 'ireach'], inplace=True)
    segment_routing = dict(zip(nseg, outseg))
    seg_groups = df.groupby('iseg')

    # segments associated with reach numbers that are first reaches
    rno1_segments = {g.rno.first(): s for s, g in seg_groups}

    segments_consistent = []
    for s, g in seg_groups:
        # since the segment is sorted,
        # rno[i+1] should == outreach[i]
        preceding_consistent = np.array_equal(g.rno.values[1:],
                                              g.outreach.values[:1])

        # check that last reach goes to same segment
        # as outseg in segment_routing
        last_outreach = g.outreach.last()
        next_segment = rno1_segments[last_outreach]
        last_consistent = next_segment == segment_routing[s]
        segments_consistent.append(preceding_consistent &
                                   last_consistent)
    return np.all(segments_consistent)


def routing_numbering_is_valid(nseg, outseg, iseg, ireach,
                               rno, outreach, increasing_nseg=True):
    """Check that routing numbering for an SFR dataset is valid.

    * verify that segment numbering is consecutive and starts at 1
        * optionally verify that segment number only increase downstream
    * verify that unique reach numbering (e.g. rno in MODFLOW 6)
        is consecutive and starts at 1
    * check that routing is consistent between segment connections
        (MODFLOW-2005 convention of nseg -> outseg)
        and reach connections (MODFLOW 6 convention based on rno)

    An additional check would be all non-outlet connections are
    listed in nseg and/or rno, but these can be assumed to be outlets
    (converted to 0) without modifying the nseg or rno.

    Parameters
    ----------
    nseg : list or 1D numpy array
    outseg : list or 1D numpy array
    iseg : list or 1D numpy array
    ireach : list or 1D numpy array
    rno : list or 1D numpy array
    outreach : list or 1D numpy array
    increasing_nseg : bool
        If True, segment numbers must also only increase downstream.

    Returns
    -------
    valid
    """
    return valid_rnos(rno) & \
           valid_nsegs(nseg, outseg, increasing=increasing_nseg) & \
           rno_nseg_routing_consistent(nseg, outseg, iseg, ireach,
                                       rno, outreach)


def routing_is_circular(fromid, toid):
    """Verify that segments or reaches never route to themselves.

    Parameters
    ----------
    fromid : list or 1D array
        e.g. COMIDS, segments, or rnos
    toid : list or 1D array
        routing connections
    """
    fromid = np.atleast_1d(fromid)
    toid = np.atleast_1d(toid)

    graph = make_graph(fromid, toid, one_to_many=False)
    paths = {fid: find_path(graph, fid) for fid in graph.keys()}
    # a fromid should not appear more than once in its sequence
    for k, v in paths.items():
        if v.count(k) > 1:
            return True
    return False


def same_sfr_numbering(reach_data1, reach_data2):
    """Compare two sets of reach data.

    Parameters
    ----------
    reach_data1 : DataFrame
        Must have columns:
        i : zero-based row
        j : zero-based column
        iseg : segment number
        ireach : reach number
    reach_data2 : DataFrame
        Must have same columns as reach_data1

    Returns
    -------
    issame : bool
        Whether the two datasets have the same numbering for i, j andn iseg/ireach.

    Notes
    -----
    k (layer) is not tested because k can be different for the same SFR package depending on the context.
    For example, a reach might have k=1 in the input file, and k=3 in the output file if
    the flux was placed in the highest active layer.

    """
    cols = ['i', 'j', 'iseg', 'ireach']
    rd1 = reach_data1[cols].sort_values(by=['iseg', 'ireach']).copy()
    rd2 = reach_data2[cols].sort_values(by=['iseg', 'ireach']).copy()
    col_equal = []
    for c in rd1.columns:
        col_equal.append(np.array_equal(rd1[c], rd2[c]))
    return np.all(col_equal)


def reach_elevations_decrease_downstream(reach_data):
    """Verify that reach elevations decrease monotonically in the downstream direction."""
    rd = reach_data.reset_index()
    elev = dict(zip(rd.rno, rd.strtop))
    dnelev = {rid: elev[rd.outreach[i]] if rd.outreach[i] != 0
    else -9999 for i, rid in enumerate(rd.rno)}
    diffs = np.array([(dnelev[i] - elev[i]) if dnelev[i] != -9999
                      else -.001 for i in rd.rno])
    return np.max(diffs) <= 0
