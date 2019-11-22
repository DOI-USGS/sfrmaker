# TODO: add unit tests for utils.py
import numpy as np

from ..utils import assign_layers


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
