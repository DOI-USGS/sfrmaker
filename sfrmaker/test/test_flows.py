import numpy as np
import pandas as pd
import pytest

from .test_routing import add_line_sequence
from ..flows import (get_inflow_locations_from_parent_model, add_to_perioddata)
from gisutils import shp2df


@pytest.fixture()
def parent_model_sfr_flow_results():
    parent_model_sfrlines = 'sfrmaker/test/data/shellmound/merasnwt_sfrlines.shp'
    rd = shp2df(parent_model_sfrlines)
    rd['Qout'] = 0.
    rd = rd.rename(columns={'iseg': 'segment', 'ireach': 'reach'}) \
        [['rno', 'segment', 'reach', 'Qout']].copy()
    rd.loc[rd.rno == 13933, 'Qout'] = 353146.667
    rd.loc[rd.rno == 11780, 'Qout'] = 3531.46667
    rd.loc[rd.rno == 11949, 'Qout'] = 3.53146667
    rd.loc[rd.rno == 11483, 'Qout'] = 353.146667
    rd.loc[rd.rno == 13070, 'Qout'] = 7062.93334
    rd.loc[rd.rno == 15682, 'Qout'] = 7.06293334
    rd.loc[rd.rno == 15684, 'Qout'] = 35314.6667
    rd['kstpkper'] = [(0, 0)] * len(rd)
    rd2 = rd.copy()
    rd2['kstpkper'] = [(1, 1)] * len(rd)
    rd2['Qout'] *= 2
    rd = rd.append(rd2).copy()
    return rd


@pytest.mark.xfail(reason="still need to replace sr with modelgrid")
def test_get_inflow_locations_from_parent_model(outdir, shellmound_grid, shellmound_sfrdata):
    parent_reach_data = 'sfrmaker/test/data/shellmound/merasnwt_sfrlines.shp'
    inset_reach_data = shellmound_sfrdata.reach_data
    inset_sr = shellmound_grid
    df = get_inflow_locations_from_parent_model(parent_reach_data, inset_reach_data,
                                                inset_sr, active_area=None
                                                )
    shellmound_sfrdata.export_cells('{}/shellmound_sfr_cells.shp'.format(outdir))
    shellmound_sfrdata.export_lines('{}/shellmound_sfr_lines.shp'.format(outdir))
    expected_line_ids = {1000005,  # yazoo
                         18046670,  # Tallahatchie River
                         938030409,  # unnamed stream along southern boundary
                         17955337,  # unnamed stream in NW corner
                         17955371,  # downstream of "Wild Bill Bayou"
                         17955445,  # unnamed stream just east of 17955281
                         17958187,  # big sunflower river
                         17956547,  # unnamed trib to BSR on west side
                         }
    symmetric_diff = expected_line_ids ^ set(df.line_id)
    assert len(symmetric_diff) == 0


def test_add_to_perioddata(shellmound_sfrdata):
    sfrd = shellmound_sfrdata  # copy.deepcopy(shellmound_sfrdata)
    rd = shellmound_sfrdata.reach_data.copy()
    line_id = dict(zip(rd.iseg, rd.line_id))
    sfr_routing = shellmound_sfrdata.segment_routing.copy()

    # routing for source hydrography
    flowline_routing = {line_id.get(k, 0): line_id.get(v, 0)
                        for k, v in sfr_routing.items()}
    nlines = 4
    seq = add_line_sequence(flowline_routing, nlines=nlines)
    flows = pd.DataFrame({'Q_avg': [100., 10., 200., 20.],
                          'per': [0, 1, 0, 1],
                          'line_id': [2, 2, 4, 4]})
    add_to_perioddata(sfrd, flows,
                      flowline_routing=flowline_routing,
                      variable='inflow',
                      line_id_column_in_data='line_id',
                      period_column_in_data='per',
                      variable_column_in_data='Q_avg')
    flows = flows.loc[flows.line_id != 2]
    assert np.allclose(sfrd.period_data['inflow'].values, flows['Q_avg'].values)
    assert np.allclose(sfrd.period_data['per'].values, flows['per'].values)
    assert rd.loc[rd.line_id == seq[-1], 'rno'].values[0] == sfrd.period_data['rno'].values[0]
