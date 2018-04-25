import numpy as np
import pandas as pd
from .utils import make_graph, find_path

def valid_rnos(rnos):
    """Check that unique reach numbers (rno in MODFLOW 6)
    are consecutive and start at 1.
    """
    sorted_reaches = sorted(rnos)
    consecutive = np.diff(sorted_reaches).sum() \
                  == len(rnos) - 1
    onebased = np.min(sorted_reaches) == 1
    return consecutive & onebased

def valid_nsegs(nsegs, outsegs=None, increasing=False):
    """Check that segment numbers are valid.

    Parameters
    ----------
    nsegs : list of segment numbers
    outsegs : list of corresponding routing connections
        Required if increasing=True.
    increasing : bool
        If True, segment numbers must also only increase downstream.
    """
    consecutive_and_onebased = valid_rnos(nsegs)
    if increasing:
        assert outsegs is not None
        graph = make_graph(nsegs, outsegs, one_to_many=False)
        monotonic = []
        for s in nsegs:
            seg_sequence = find_path(graph, s)
            monotonic.append(np.all(np.array(seg_sequence).diff() > 0))
        return consecutive_and_onebased & monotonic
    else:
        return consecutive_and_onebased

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

def routing_is_valid(nseg, outseg, iseg, ireach, rno, outreach):
    """Check that routing for an SFR dataset is valid.

    * verify that segment numbering is consecutive and starts at 1
        * optionally verify that segment number only increase downstream
    * verify that unique reach numbering (e.g. rno in MODFLOW 6)
        is consecutive and starts at 1

    :param nseg:
    :param outseg:
    :param iseg:
    :param ireach:
    :param rno:
    :param outreach:
    :return:
    """




