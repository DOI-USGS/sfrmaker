import numpy as np
from sfrmaker.routing import make_graph, get_upsegs

from ..checks import routing_is_circular, valid_nsegs
from ..routing import (get_next_id_in_subset, renumber_segments, find_path,
                       get_previous_ids_in_subset)


def add_line_sequence(routing, nlines=4, string_ids=False):
    headwater_line_ids = set(routing.keys()).difference(routing.values())
    # add in some lines upstream of a headwater
    n = nlines
    new_lines = set(range(n + 2)).difference(routing)
    new_routing = {}
    lid = new_lines.pop()
    sequence = [lid]
    for i in range(n):
        nextid = new_lines.pop()
        new_routing[lid] = nextid
        lid = nextid
        sequence.append(nextid)
    end = headwater_line_ids.pop()
    new_routing[lid] = end
    if string_ids:
        new_routing = {str(k): str(v) for k, v in new_routing.items()}
    sequence.append(end)
    routing.update(new_routing)
    return sequence


def test_get_next_id_in_subset(shellmound_sfrdata):
    rd = shellmound_sfrdata.reach_data.copy()
    line_id = dict(zip(rd.iseg, rd.line_id))
    sfr_routing = shellmound_sfrdata.segment_routing.copy()

    # routing for source hydrography
    routing = {str(line_id.get(k, 0)): line_id.get(v, '0')
               for k, v in sfr_routing.items()}
    nlines = 4
    seq = add_line_sequence(routing, nlines=nlines, string_ids=True)

    ids = {str(lid) for lid in set(range(nlines + 1))}
    result = get_next_id_in_subset(rd.line_id, routing, ids)
    assert set(result).difference(rd.line_id) == {'0'}
    assert set(result) == {'0', seq[-1]}
    assert len(result) == len(range(nlines + 1))


def test_get_previous_ids_in_subset():
    routing = {17955471: 17957799,
               17957799: 17955535,
               17955535: 17955537,
               17955537: 17955541,
               17955541: 17955663,
               17955663: 17955689,
               17955543: 17955663,
               17955689: 10001,
               10001: 0
               }
    # delete some segments
    remove_ids = {17955663, 17955689, 17955541, 17955471, 10001}
    subset = set(routing.keys()).difference(remove_ids)

    result = get_previous_ids_in_subset(subset, routing, remove_ids)
    # verify that the new ids are in the subset
    assert result == {17955537, 17955543}


def test_renumber_segments():
    non_consecutive_nseg = [1, 2, 3, 4, 5, 6, 7, 13, 14, 15, 20, 27, 31, 32, 33, 34, 35, 42, 43, 46, 48, 49, 52, 53, 54, 58, 70, 71, 84,
            99, 100, 101, 114, 115, 116, 132, 133, 167, 168, 169, 170, 180, 189, 190, 203, 204, 219, 220, 221, 233, 234,
            235, 245, 246, 247, 248, 249, 250, 259, 275, 276, 277, 278, 285, 294, 295, 304, 312, 325, 326, 327, 328,
            329, 330, 331, 332, 343, 344, 366, 367, 420, 446, 447, 448, 449, 450, 451, 452]
    non_consecutive_outseg = [2, 3, 4, 5, 6, 7, 13, 14, 15, 20, 27, 31, 32, 33, 34, 35, 42, 43, 46, 48, 49, 52, 53, 54, 58, 70, 71, 84,
              99, 100, 101, 114, 115, 116, 132, 133, 167, 168, 169, 170, 180, 189, 190, 203, 204, 219, 220, 221, 233,
              234, 235, 245, 246, 247, 248, 249, 250, 259, 275, 276, 277, 278, 285, 294, 295, 304, 312, 325, 326, 327,
              328, 329, 330, 331, 332, 343, 344, 366, 367, 420, 446, 447, 448, 449, 450, 451, 452, 0]
    assert not routing_is_circular(non_consecutive_nseg, non_consecutive_outseg)
    assert not valid_nsegs(non_consecutive_nseg, non_consecutive_outseg)
    nseg = list(range(1, len(non_consecutive_nseg) + 1))
    renumber = dict(zip(non_consecutive_nseg, nseg))
    outseg = [renumber.get(s, 0) for s in non_consecutive_outseg]
    assert not routing_is_circular(nseg, outseg)
    assert valid_nsegs(nseg, outseg)
    new_numbering1 = renumber_segments(non_consecutive_nseg, non_consecutive_outseg)
    new_numbering2 = renumber_segments(nseg, outseg)
    nseg1 = [new_numbering1[s] for s in non_consecutive_nseg]
    outseg1 = [new_numbering1[s] for s in non_consecutive_outseg]
    nseg2 = [new_numbering2[s] for s in nseg]
    outseg2 = [new_numbering2[s] for s in outseg]
    assert nseg1 == nseg2
    assert outseg1 == outseg2
    assert not routing_is_circular(nseg1, outseg1)
    assert valid_nsegs(nseg1, outseg1)


def test_get_upsegs(sfr_test_numbering):
    rd, sd = sfr_test_numbering
    graph = dict(zip(sd.nseg, sd.outseg))
    graph_r = make_graph(list(graph.values()),
                         list(graph.keys()))
    upsegs = []
    for s in sd.nseg:
        upsegs.append(get_upsegs(graph_r, s))
    assert isinstance(upsegs[1], set)
    assert upsegs[1] == {1, 3, 4, 5, 6, 7, 8, 9}


def test_find_path():
    n = 1e5
    rnos = np.arange(1, n+1)
    out_rnos = np.arange(2, n+2)
    out_rnos[-1] = 0
    routing = dict(zip(rnos, out_rnos))
    path = find_path(routing, start=1)
    assert path[0] == 1
    assert path[-1] == 0
