import platform
import sys

import numpy as np
import pandas as pd
import pytest

from .test_routing import add_line_sequence
from ..flows import (get_inflow_locations_from_parent_model, get_inflows_from_parent_model,
add_to_perioddata, add_to_segment_data)
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


def test_get_inflow_locations_from_parent_model(outdir, shellmound_grid, shellmound_sfrdata):
    parent_reach_data = 'sfrmaker/test/data/shellmound/merasnwt_sfrlines.shp'
    inset_reach_data = shellmound_sfrdata.reach_data
    inset_sr = shellmound_grid
    df = get_inflow_locations_from_parent_model(parent_reach_data, inset_reach_data,
                                                inset_sr, active_area=None
                                                )
    shellmound_sfrdata.export_cells('{}/shellmound_sfr_cells.shp'.format(outdir))
    shellmound_sfrdata.export_lines('{}/shellmound_sfr_lines.shp'.format(outdir))
    expected_line_ids = {'1000005',  # yazoo
                         '18046670',  # Tallahatchie River
                         '938030409',  # unnamed stream along southern boundary
                         '17955337',  # unnamed stream in NW corner
                         '17955371',  # downstream of "Wild Bill Bayou"
                         '17955445',  # unnamed stream just east of 17955281
                         '17958187',  # big sunflower river
                         '17956547',  # unnamed trib to BSR on west side
                         }
    symmetric_diff = expected_line_ids ^ set(df.line_id)
    assert len(symmetric_diff) == 0


@pytest.mark.skip(reason="still working on this feature")
def test_get_inflows_from_parent_model(shellmound_sfrdata):
    parent_reach_data = 'sfrmaker/test/data/shellmound/merasnwt_sfrlines.shp'
    inset_reach_data = shellmound_sfrdata.reach_data
    result = get_inflows_from_parent_model(parent_reach_data, inset_reach_data,
                                           #mf2005_parent_sfr_outputfile, mf6_parent_sfr_budget_file,
                                           #inset_grid, active_area=None
                                           )

@pytest.fixture(scope='function')
def data(shellmound_sfrdata):
    sfr_routing = shellmound_sfrdata.segment_routing.copy()
    data = {
        'sfrd': shellmound_sfrdata,  # copy.deepcopy(shellmound_sfrdata)
        'rd': shellmound_sfrdata.reach_data.copy()
    }
    line_id = dict(zip(data['rd'].iseg, data['rd'].line_id))

        # routing for source hydrography
    data['flowline_routing'] = {line_id.get(k, 0): line_id.get(v, 0)
                            for k, v in sfr_routing.items()}
    nlines = 4
    data['seq'] = add_line_sequence(data['flowline_routing'], nlines=nlines)
    return data

def test_add_to_perioddata1(data):
    # two inflows applied upstream of model on same path;
    # results should be summed an applied to first line_id in model
    sfrd = data['sfrd']
    flowline_routing = data['flowline_routing']
    rd = data['rd']
    flows = pd.DataFrame({'Q_avg': [100., 10., 200., 20.],
                        'per': [0, 1, 0, 1],
                        'line_id': [2, 2, 4, 4]})
    add_to_perioddata(sfrd, flows,
                    flowline_routing=flowline_routing,
                    variable='inflow',
                    line_id_column='line_id',
                    period_column='per',
                    data_column='Q_avg',
                    one_inflow_per_path=False)
    flows = flows.groupby('per').sum()
    assert np.allclose(sfrd.period_data['inflow'].values, flows['Q_avg'].values)
    assert np.allclose(sfrd.period_data.index.get_level_values(0), 
                    flows.index.values)
    assert rd.loc[rd.line_id == data['seq'][-1], 'rno'].values[0] ==\
        sfrd.period_data.index.levels[1].values[0]

@pytest.mark.skipif((platform.system() == 'Windows') & (sys.version_info[:2] == (3, 12)), 
                    reason="inscrutable error")
def test_add_to_perioddata2(data):
    # two inflows applied upstream of model on different paths;
    # inflows should be applied to the first lines in the model 
    # along their respective paths
    sfrd = data['sfrd']
    flowline_routing = data['flowline_routing']
    rd = data['rd']
    flowline_routing['6'] = '18046688'
    flows = pd.DataFrame({'Q_avg': [100., 10., 200., 20.],
                        'per': [0, 1, 0, 1],
                        'line_id': [6, 6, 4, 4]})
    add_to_perioddata(sfrd, flows,
                    flowline_routing=flowline_routing,
                    variable='inflow',
                    line_id_column='line_id',
                    period_column='per',
                    data_column='Q_avg',
                    one_inflow_per_path=False)
    assert np.allclose(sfrd.period_data['inflow'].sort_values(), flows['Q_avg'].sort_values())
    assert np.allclose(sfrd.period_data.index.get_level_values(0), 
                    sorted(flows['per']))
    lines = {flowline_routing['6'], data['seq'][-1]}
    rnos = rd.loc[rd.line_id.isin(lines), 'rno']
    assert not set(sfrd.period_data.index.levels[1]).difference(rnos)
    assert len(set(sfrd.period_data.index.levels[1])) == 2

def test_add_to_perioddata3(data):
    # two inflows applied along same path in model;
    # one_inflow_per_path=True
    # inflows should be summed and applied to most downstream reach
    # test updating of period_data
    # sfrd.period_data wasn't reset, so should contain
    # updated values of 300 and 30 at reach 354
    # and existing values of 100 and 10 at 395
    sfrd = data['sfrd']
    flowline_routing = data['flowline_routing']
    flows = pd.DataFrame({'Q_avg': [100., 10., 200., 20.],
                        'per': [0, 1, 0, 1],
                        'line_id': [2, 2, 4, 4]})
    add_to_perioddata(sfrd, flows,
                    flowline_routing=flowline_routing,
                    variable='inflow',
                    line_id_column='line_id',
                    period_column='per',
                    data_column='Q_avg',
                    one_inflow_per_path=True)
    assert np.allclose(sfrd.period_data['inflow'].values, [300., 100., 30., 10.])
    # these reach numbers are liable to change
    rnos = [
        sfrd.reach_data.loc[sfrd.reach_data['line_id'] == data['seq'][-1], 'rno'].min(),
        sfrd.reach_data.loc[sfrd.reach_data['line_id'] == '1000005', 'rno'].min()
    ]
    assert np.allclose(sfrd.period_data.index.levels[1], rnos)

def test_add_to_perioddata3(data):
    # two inflows applied along same path in model;
    # one_inflow_per_path=False
    # inflows should be applied to their respective reaches
    sfrd = data['sfrd']
    flowline_routing = data['flowline_routing']
    flows = pd.DataFrame({'Q_avg': [100., 10., 200., 20.],
                        'per': [0, 1, 0, 1],
                        'line_id': ['17955337', '17955337', 
                                    flowline_routing['17955337'], 
                                    flowline_routing['17955337']
                                    ]})
    sfrd._period_data = None
    add_to_perioddata(sfrd, flows,
                    flowline_routing=flowline_routing,
                    variable='inflow',
                    line_id_column='line_id',
                    period_column='per',
                    data_column='Q_avg',
                    one_inflow_per_path=False)
    assert np.allclose(sfrd.period_data['inflow'].sort_index(level='rno').values, 
                    [100., 10., 200., 20.])
    assert np.allclose(sfrd.period_data.index.levels[1], [354, 401])
    
    
def test_add_to_segment_data(shellmound_sfrdata):
    sfrd = shellmound_sfrdata  # copy.deepcopy(shellmound_sfrdata)
    rd = shellmound_sfrdata.reach_data.copy()
    line_id = dict(zip(rd.iseg, rd.line_id))
    segment = {v:k for k, v in line_id.items()}
    sfr_routing = shellmound_sfrdata.segment_routing.copy()

    # routing for source hydrography
    flowline_routing = {line_id.get(k, 0): line_id.get(v, 0)
                        for k, v in sfr_routing.items()}
    nlines = 4
    seq = add_line_sequence(flowline_routing, nlines=nlines)
    line = sfrd.reach_data.line_id.values[0]
    flows = pd.DataFrame({'Q_avg': [100., 10., 200., 20., 51, 52, 53, 54],
                          'per': [0, 1, 0, 1, 2, 3, 4, 5],
                          'line_id': [2, 2, 4, 4, line, line, line, line]})
    sd1 = sfrd.segment_data.copy()
    sd1.index = pd.MultiIndex.from_tuples(zip(sd1.per, sd1.nseg), names=['per', 'nseg'])
    add_to_segment_data(sfrd, flows,
                        flowline_routing=flowline_routing,
                        variable='flow',
                        line_id_column='line_id',
                        period_column='per',
                        data_column='Q_avg')
    sd2 = sfrd.segment_data.copy()
    sd2.index = pd.MultiIndex.from_tuples(zip(sd2.per, sd2.nseg), names=['per', 'nseg'])
    #flows = flows.loc[~flows.line_id.isin([2])]
    #flows = flows.groupby(level=(0, 1)).sum()
    flows['nseg'] = [segment.get(l, segment[seq[-1]]) for l in flows.line_id]
    flow_sums = flows.groupby(['per', 'nseg']).first()
    flow_sums['Q_avg'] = flows.groupby(['per', 'nseg'])['Q_avg'].sum()
    #flows.index = pd.MultiIndex.from_tuples(zip(flows.per, flows.nseg), names=['per', 'nseg'])
    assert np.allclose(sd2.loc[flow_sums.index, 'flow'], flow_sums.Q_avg)
    assert not sd2.isna().any().any()
    pd.testing.assert_frame_equal(sd1.drop('flow', axis=1),
                                  sd2.loc[sd1.index].drop('flow', axis=1),
                                  check_dtype=False)
