import sys

sys.path.append('..')
import numpy as np
from sfrmaker.utils import arbolate_sum
from sfrmaker.routing import get_upsegs, make_graph


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


def test_asum(sfr_test_numbering):
    rd, sd = sfr_test_numbering
    graph = dict(zip(sd.nseg, sd.outseg))
    lengths = dict(zip(sd.nseg, np.arange(len(sd))))
    asum = arbolate_sum(sd.nseg, lengths, graph)
    assert (asum[1] == 6) & (asum[2] == np.arange(len(sd)).sum())
