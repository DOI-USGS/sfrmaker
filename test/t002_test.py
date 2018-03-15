import sys
sys.path.append('..')
import numpy as np
import pandas as pd
from sfrmaker.utils import arbolate_sum, get_upsegs, make_graph

def create_sfr_data():
    rd = pd.DataFrame()
    rd['i'] = [3, 4, 5,
              7, 8, 9,
              0, 1, 2,
              4, 4, 5,
              0, 0, 0,
              3, 4, 5,
              0, 1, 2,
              4, 5, 6,
              2, 2, 2]
    rd['j'] = [0, 1, 2,
              6, 6, 6,
              6, 6, 6,
              3, 4, 5,
              9, 8, 7,
              6, 6, 6,
              0, 0, 0,
              6, 6, 6,
              9, 8, 7]
    rd['iseg'] = sorted(list(range(1, 10)) * 3)
    rd['ireach'] = [1, 2, 3] * 9

    sd = pd.DataFrame()
    sd['nseg'] = range(1, 10)
    sd['outseg'] = [4, 0, 6, 8, 3, 8, 1, 2, 8]
    return rd, sd

def test_get_upsegs():
    rd, sd = create_sfr_data()
    graph = dict(zip(sd.nseg, sd.outseg))
    graph_r = make_graph(list(graph.values()), list(graph.keys()))
    upsegs = []
    for s in sd.nseg:
        upsegs.append(get_upsegs(graph_r, s))
    assert upsegs[1] == {1, 3, 4, 5, 6, 7, 8, 9}

def test_asum():

    rd, sd = create_sfr_data()

    graph = dict(zip(sd.nseg, sd.outseg))
    lengths = dict(zip(sd.nseg, np.arange(len(sd))))
    asum = arbolate_sum(sd.nseg, lengths, graph)
    assert (asum[1] == 6) & (asum[2] == np.arange(len(sd)).sum())

if __name__ == '__main__':
    test_get_upsegs()
    test_asum()