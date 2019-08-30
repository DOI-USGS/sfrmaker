import pytest
import flopy
import sfrmaker


@pytest.fixture(scope="function")
def sfrdata(sfr_test_numbering):
    rd, sd = sfr_test_numbering
    sd['width1'] = 1
    sd['width2'] = 1
    return sfrmaker.sfrdata(reach_data=rd, segment_data=sd)


def test_const(shellmound_sfrdata, sfrdata):
    assert shellmound_sfrdata._lenuni == 2  # meters
    assert shellmound_sfrdata._itmuni == 4  # days
    assert shellmound_sfrdata.const == 86400

    sfrdata._model_length_units = 'feet'
    sfrdata._model_time_units = 'days'
    assert sfrdata.const == 86400 * 1.486


def test_create_mf6sfr(shellmound_sfrdata, shellmound_model):
    mf6sfr = shellmound_sfrdata.create_mf6sfr(model=shellmound_model)
    assert isinstance(mf6sfr, flopy.mf6.ModflowGwfsfr)
    assert hasattr(shellmound_model, 'sfr')