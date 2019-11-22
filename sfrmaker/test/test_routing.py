from ..routing import get_next_id_in_subset


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
