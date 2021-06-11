"""
Tests for units.py module
"""
import numpy as np
import flopy
from ..units import (convert_flux_units, convert_length_units,
                     convert_time_units, get_model_length_units)


def test_convert_flux():
    result = convert_flux_units('inches', 'years',
                                'feet', 'days')
    assert np.allclose(result, 1/12 * 1/365.25)
    result = convert_flux_units('feet', 'second',
                                'meters', 'days')
    assert np.allclose(result, 0.3048 * 86400)


def test_convert_length_units():
    assert np.allclose(convert_length_units(2, 1), 1/.3048)
    assert np.allclose(convert_length_units(1, 2), .3048)
    assert np.allclose(convert_length_units('meters', 'feet'), 1/.3048)
    assert np.allclose(convert_length_units('feet', 'meters'), .3048)
    assert np.allclose(convert_length_units('m', 'ft'), 1/.3048)
    assert np.allclose(convert_length_units('ft', 'm'), .3048)
    assert np.allclose(convert_length_units(None, 'm'), 1.)
    assert np.allclose(convert_length_units('millimeters', 'meters'), 1/1000)
    assert np.allclose(convert_length_units('meters', 'millimeters'), 1000)
    assert np.allclose(convert_length_units('meters', 'km'), 0.001)
    assert np.allclose(convert_length_units('kilometers', 'meters'), 1000)
    assert np.allclose(convert_length_units('kilometers', 'cm'), 1000*100)


def test_convert_time_units():
    assert np.allclose(convert_time_units(4, 1), 86400)
    assert np.allclose(convert_time_units('days', 'seconds'), 86400)
    assert np.allclose(convert_time_units('d', 's'), 86400)
    assert np.allclose(convert_time_units(1, 4), 1/86400)
    assert np.allclose(convert_time_units(None, 'd'), 1.)
    assert np.allclose(convert_time_units(5, 4), 365.25)
    assert np.allclose(convert_time_units(4, 5), 1/365.25)
    assert np.allclose(convert_time_units('years', 'days'), 365.25)


def test_get_model_length_units():
    mf = flopy.modflow.Modflow()
    dis = flopy.modflow.ModflowDis(mf)
    units = get_model_length_units(mf)
    assert units == 'meters'  # flopy default
    dis.lenuni = 0
    units = get_model_length_units(mf)
    assert units == 'undefined'  # flopy default

    sim = flopy.mf6.MFSimulation()
    #tdis = flopy.mf6.modflow.mftdis.ModflowTdis(sim)
    gwf = flopy.mf6.ModflowGwf(sim)
    dis = flopy.mf6.modflow.mfgwfdis.ModflowGwfdis(gwf)
    units = get_model_length_units(gwf)
    assert units == 'undefined'  # flopy default
    dis.length_units = 'feet'
    units = get_model_length_units(gwf)
    assert units == 'feet'
