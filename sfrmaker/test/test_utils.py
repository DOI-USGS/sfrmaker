# TODO: add unit tests for utils.py
import numpy as np
import pandas as pd
import pytest
from sfrmaker.utils import (assign_layers, width_from_arbolate_sum,
                            arbolate_sum, make_config_summary)


def test_assign_layers(shellmound_sfrdata, shellmound_model):
    
    reach_data = shellmound_sfrdata.reach_data.copy()
    reach_data.reset_index(inplace=True, drop=True)
    botm = shellmound_model.dis.botm.array.copy()
    idomain = shellmound_model.dis.idomain.array.copy()
    
    # test cases
    # invalid bottom in stack of all inactive celle
    i, j = reach_data['i'], reach_data['j']
    original_botm_at_rno0 = botm[-1, i[0], j[0]]
    reach_data.loc[0, 'strtop'] = original_botm_at_rno0 - 10
    idomain[:, i[0], j[0]] = 0
    # active cells on top of inactive cells
    # lowest active layer should be pushed downward,
    # not lower (inactive) layer
    reach_data.loc[1, 'strtop'] = botm[-1, i[1], j[1]] - 10
    lowest_active_layer = 9
    botm[lowest_active_layer+1:, i[1], j[1]] = botm[-1, i[1], j[1]]
    idomain[lowest_active_layer+1:, i[1], j[1]] = 0
    
    k, new_botm = assign_layers(reach_data,
                                     botm_array=botm,
                                     idomain=idomain)
    reach_data['k'] = k
    active_reaches = idomain[k, i, j] > 0
    shellmound_model.dis.botm = new_botm
    assert np.all(reach_data.loc[active_reaches].strtop.values >
                  shellmound_model.dis.botm.array[reach_data.k.values,
                                                  reach_data.i.values,
                                                  reach_data.j.values][active_reaches])
    # model bottom at reach 0 should have been left alone (all inactive cells)
    assert new_botm[-1, i[0], j[0]] == original_botm_at_rno0
    # reach 1 should be in lowest active layer, not lowest layer
    assert reach_data.loc[1, 'k'] == lowest_active_layer
    # inactive layers below should all be zero-thickness
    assert np.allclose(new_botm[lowest_active_layer:, i[1], j[1]], new_botm[-1, i[1], j[1]])
    
    # no idomain supplied
    k2, new_botm2 = assign_layers(reach_data,
                                  botm_array=botm)
    # for now, return a 2D (model bottom) if no idomain
    # (for backwards compatibility)
    assert len(new_botm2.shape) == 2
    reach_data['k'] = k2
    botm[-1] = new_botm2
    assert np.all(reach_data.strtop.values > botm[k2, i, j])


@pytest.mark.parametrize('asum,expected',
                         ((1e6, 124.7*.3048),)
                         )
def test_width_from_arbolate_sum_defaults(asum, expected):
    result = width_from_arbolate_sum(asum)
    assert np.allclose(result, expected, rtol=1e-4)


@pytest.mark.parametrize('asum,a,b,input_units,output_units,expected',
                         ((1e6, 0.1193, 0.5032, 'meters', 'feet', 124.7),
                          (1e3, 0.1193, 0.5032, 'km', 'meters', 124.7*.3048),
                          (1e6, 0.0628, 0.5099, 'meters', 'feet', 72.00),
                          (1e3, 0.0628, 0.5099, 'km', 'meters', 72.00*.3048),
                          (0, 0.0628, 0.5099, 'km', 'meters', 1),
                          ([1e3], 0.0628, 0.5099, 'km', 'meters', [72.00*.3048]),
                          (np.array([1e3]), 0.0628, 0.5099, 'km', 'meters', [72.00*.3048]),
                          (np.array([1e3, 0]), 0.0628, 0.5099, 'km', 'meters', [72.00*.3048, 1]))
                         )
def test_width_from_arbolate_sum(asum, a, b, input_units, output_units, expected):
    result = width_from_arbolate_sum(asum, a, b, minimum_width=1,
                                     input_units=input_units, output_units=output_units)
    assert np.allclose(result, expected, rtol=1e-4)


def test_asum(sfr_test_numbering):
    rd, sd = sfr_test_numbering
    graph = dict(zip(sd.nseg, sd.outseg))
    lengths = dict(zip(sd.nseg, np.arange(len(sd))))
    asum = arbolate_sum(sd.nseg, lengths, graph)
    assert (asum[1] == 6) & (asum[2] == np.arange(len(sd)).sum())


def test_make_config_summary():
    results = make_config_summary()
