from ..checks import routing_is_circular, valid_nsegs
from ..routing import (get_next_id_in_subset, renumber_segments)


def add_line_sequence(routing, nlines=4):
    headwater_line_ids = set(routing.keys()).difference(routing.values())
    # add in some lines upstream of a headwater
    n = nlines
    new_lines = set(range(n + 2)).difference(routing)
    new_routing = {}
    id = new_lines.pop()
    sequence = [id]
    for i in range(n):
        nextid = new_lines.pop()
        new_routing[id] = nextid
        id = nextid
        sequence.append(nextid)
    end = headwater_line_ids.pop()
    new_routing[id] = end
    sequence.append(end)
    routing.update(new_routing)
    return sequence


def test_get_next_id_in_subset(shellmound_sfrdata):
    rd = shellmound_sfrdata.reach_data.copy()
    line_id = dict(zip(rd.iseg, rd.line_id))
    sfr_routing = shellmound_sfrdata.segment_routing.copy()

    # routing for source hydrography
    routing = {line_id.get(k, 0): line_id.get(v, 0)
               for k, v in sfr_routing.items()}
    nlines = 4
    seq = add_line_sequence(routing, nlines=nlines)

    result = get_next_id_in_subset(rd.line_id, routing, set(range(nlines + 2)))
    assert set(result).difference(rd.line_id) == {0}
    assert set(result) == {0, seq[-1]}
    assert len(result) == len(range(nlines + 2))


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
