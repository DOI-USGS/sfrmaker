# TODO: add unit tests for utils.py
import numpy as np
import pytest
from sfrmaker.utils import (assign_layers, width_from_arbolate_sum,
                            arbolate_sum, make_config_summary)


def test_assign_layers(shellmound_sfrdata, shellmound_model):
    layers, new_botm = assign_layers(shellmound_sfrdata.reach_data,
                                     botm_array=shellmound_model.dis.botm.array)
    shellmound_sfrdata.reach_data['k'] = layers
    botm = shellmound_model.dis.botm.array.copy()
    botm[-1] = new_botm
    shellmound_model.dis.botm = botm
    rd = shellmound_sfrdata.reach_data
    assert np.all(rd.strtop.values >
                  shellmound_model.dis.botm.array[rd.k.values,
                                                  rd.i.values,
                                                  rd.j.values])


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
