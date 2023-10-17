import os
from pathlib import Path
import copy
import io
from packaging import version
import flopy
import numpy as np
import pandas as pd
import pytest
from gisutils import shp2df
import sfrmaker
from sfrmaker.mf5to6 import cellids_to_kij


@pytest.fixture(scope='function')
def period_data(shellmound_sfrdata):
    # shellmound_sfrdata = copy.deepcopy(shellmound_sfrdata)
    perdata = shellmound_sfrdata.period_data
    rd = shellmound_sfrdata.reach_data
    df = pd.DataFrame({'per': [0, 0, 1, 1],
                       'rno': [rd.rno.values[0], rd.rno.values[100]] * 2,
                       'inflow':[100., 1000., 20., 200.]})
    df.set_index(['per', 'rno'], inplace=True)
    perdata['inflow'] = df['inflow']
    return perdata


@pytest.fixture(scope='function')
def shellmound_sfrdata_with_period_data(shellmound_sfrdata, period_data):
    # shellmound_sfrdata = copy.deepcopy(shellmound_sfrdata)
    shellmound_sfrdata._period_data = period_data
    return shellmound_sfrdata


def test_init(tylerforks_sfrdata, tylerforks_model_grid, tylerforks_sfrmaker_grid_from_flopy):
    rd = tylerforks_sfrdata.reach_data.copy()
    sd = tylerforks_sfrdata.segment_data.copy()
    grid = copy.deepcopy(tylerforks_model_grid)
    sfrdata = sfrmaker.SFRData(reach_data=rd, segment_data=sd, grid=grid,
                               isfr=tylerforks_sfrmaker_grid_from_flopy.isfr)
    assert isinstance(sfrdata, sfrmaker.SFRData)
    pd.testing.assert_frame_equal(sfrdata.reach_data, rd,
                                  check_dtype=False,
                                  rtol=0.01
                                  #check_less_precise=2
                                  )
    pd.testing.assert_frame_equal(sfrdata.segment_data, sd,
                                  check_dtype=False)
    assert isinstance(sfrdata.grid, sfrmaker.StructuredGrid)
    assert sfrdata.grid == tylerforks_sfrmaker_grid_from_flopy
    j=2


def test_const(shellmound_sfrdata, sfr_testdata):
    assert shellmound_sfrdata._lenuni == 2  # meters
    assert shellmound_sfrdata._itmuni == 4  # days
    assert shellmound_sfrdata.const == 86400

    sfr_testdata.model_length_units = 'feet'
    sfr_testdata.model_time_units = 'days'
    assert sfr_testdata.const == 86400 * 1.486


def test_empty_period_data(shellmound_sfrdata):
    # shellmound_sfrdata = copy.deepcopy(shellmound_sfrdata)
    perdata = shellmound_sfrdata.period_data
    assert isinstance(perdata, pd.DataFrame)
    assert len(perdata) == 0
    assert set(perdata.columns) == \
           {'iseg', 'evaporation', 'inflow', 'status',
            'stage', 'runoff', 'rainfall',
            'ireach', 'icalc'}
    assert perdata.index.names == ['per', 'rno']


def test_perioddata(shellmound_sfrdata_with_period_data):
    sfrd = shellmound_sfrdata_with_period_data
    perdata = sfrd.period_data
    assert perdata.index.names == ['per', 'rno']
    assert np.array_equal(perdata.index.levels[0], [0, 1])
    assert np.array_equal(perdata.index.levels[1], [1, 101])
    assert np.array_equal(perdata.inflow, [100., 1000., 20., 200.])
    

def test_write_perioddata(shellmound_sfrdata_with_period_data, outdir):
    sfrd = shellmound_sfrdata_with_period_data
    outfile = '{}/test_mf6_sfr_period_data.sfr'.format(outdir)
    sfrd.write_package(outfile, version='mf6')
    with open(outfile) as src:
        rno_values = []
        txt_values = []
        q_values = []
        for line in src:
            if 'begin period' in line.lower():
                for i in range(2):
                    rno, txt, q = next(src).strip().split()
                    rno_values.append(int(rno))
                    txt_values.append(txt)
                    q_values.append(float(q))
    assert np.array_equal(rno_values, sfrd.period_data.index.get_level_values(1))
    assert np.allclose(q_values, sfrd.period_data['inflow'].values)
    assert set(txt_values) == {'inflow'}
    out_tables = Path(outdir / 'tables')
    out_tables.mkdir(exist_ok=True, parents=True)
    sfrd.write_tables(basename=outdir / 'tables/shellmound')
    written_pd = pd.read_csv(outdir / 'tables/shellmound_sfr_period_data.csv')
    assert 'rno' in written_pd.columns and 'per' in written_pd.columns


def test_export_period_data(shellmound_sfrdata_with_period_data, outdir):
    sfrd = shellmound_sfrdata_with_period_data
    outfile = '{}/test_mf6_sfr_period_data_inflow.shp'.format(outdir)
    sfrd.export_period_data(outfile)
    df = shp2df(outfile)
    nodes = dict(zip(sfrd.reach_data.rno, sfrd.reach_data.node))
    pers = [int(c.strip('inflow')) for c in df.columns if 'inflow' in c]
    assert set(pers) == set(sfrd.period_data.index.levels[0])
    assert set(df['rno']) == set(sfrd.period_data.index.levels[1])
    results = pd.concat([df['0inflow'], df['1inflow']]).values
    assert np.allclose(results, sfrd.period_data['inflow'].values)
    assert np.array_equal(df.node.values, np.array([nodes[rno] for rno in df.rno], dtype=int))

    # check that export still works if there are multiple items in a reach
    sfrd._period_data = pd.concat([sfrd.period_data, sfrd.period_data])
    sfrd.export_period_data(outfile)
    df = shp2df(outfile)
    assert np.allclose(sorted(pd.concat([df['0inflow'], df['1inflow']]).values),
                       sorted(sfrd.period_data.groupby(['rno', 'per']).sum().inflow.values))


@pytest.fixture(scope="function")
def mf6sfr(shellmound_sfrdata, shellmound_model):
    return shellmound_sfrdata.create_mf6sfr(model=shellmound_model)


@pytest.fixture(scope="function")
def mf6sfr_outfile(mf6sfr, shellmound_model, outdir):
    m = shellmound_model
    mf6sfr.write()
    outfile = os.path.join(m.model_ws, mf6sfr.filename)
    return outfile


def test_create_mf6sfr(mf6sfr, shellmound_sfrdata, shellmound_model):
    assert isinstance(mf6sfr, flopy.mf6.ModflowGwfsfr)
    assert hasattr(shellmound_model, 'sfr')

    # test that packagedata has same information as sfrmaker.sfrdata
    packagedata = pd.DataFrame(mf6sfr.packagedata.array)
    k, i, j = cellids_to_kij(packagedata['cellid'], drop_inactive=False)
    packagedata['k'] = k
    packagedata['i'] = i
    packagedata['j'] = j
    packagedata.drop('cellid', axis=1, inplace=True)
    for col in packagedata.columns:
        if packagedata[col].dtype == object:
            packagedata[col] = pd.to_numeric(packagedata[col])
    reach_data = shellmound_sfrdata.reach_data
    rno_col = {'rno', 'ifno'}.intersection(packagedata.columns).pop()
    assert np.array_equal(packagedata[rno_col].values + 1, reach_data['rno'].values)
    packagedata['iseg'] = reach_data['iseg']
    segment_data = shellmound_sfrdata.segment_data
    packagedata_groups = packagedata.groupby('iseg').mean()
    for col in packagedata.columns:
        if col in shellmound_sfrdata.mf5names:
            mf2005col = shellmound_sfrdata.mf5names[col]
            if mf2005col in reach_data.columns:
                c1 = packagedata[col].values
                c2 = reach_data[mf2005col].values
                if col == 'rno':
                    c1 += 1
                assert np.allclose(c1, c2, rtol=0.001), \
                    print('values in packagedata.{} != sfrdata.reach_data.{}'.format(col, mf2005col))
            elif mf2005col in segment_data.columns:
                assert np.allclose(packagedata_groups[col].values, segment_data[mf2005col].values, rtol=0.001), \
                    print('values in packagedata.{} != sfrdata.segment_data.{}'.format(col, mf2005col))


def test_flopy_mf6sfr_outfile(mf6sfr, mf6sfr_outfile):
    assert os.path.exists(mf6sfr_outfile)
    # check that the written package data matches in memory package data
    pkdata1 = mf6sfr.packagedata.array.copy()
    mf6sfr.load()
    for col in pkdata1.dtype.names:
        if col == 'cellid':
            continue
        elif pkdata1[col].dtype == object:
            c1 = pd.to_numeric(pkdata1[col])
            c2 = pd.to_numeric(mf6sfr.packagedata.array[col])
        else:
            c1 = pkdata1[col]
            c2 = mf6sfr.packagedata.array[col]
        assert np.allclose(c1, c2)


def test_ibound_representation_of_idomain(shellmound_sfrdata, shellmound_model):
    ibound = shellmound_sfrdata.modflow_sfr2.parent.bas6.ibound.array
    idomain = shellmound_model.dis.idomain.array
    assert np.array_equal(ibound, idomain)


@pytest.mark.xfail(version.parse(flopy.__version__) <= version.parse('3.3.0'),
                   reason="")
def test_write_mf6_package(shellmound_sfrdata, mf6sfr, outdir):
    sfr_package_file = os.path.join(outdir, 'test.package_file.sfr')
    shellmound_sfrdata.write_package(filename=sfr_package_file, version='mf6')
    with open(sfr_package_file) as src:
        for line in src:
            if 'budget' in line.lower() or 'stage' in line.lower():
                _, _, fpath = line.strip().split()
                path, fname = os.path.split(fpath)
                abspath = os.path.abspath(os.path.join(outdir, fname))
                # path should be relative
                assert path == ''
                assert fname.replace('.cbc', '').replace('.stage.bin', '') == \
                       os.path.split(sfr_package_file)[1]
            elif 'unit_conversion' in line.lower():
                _, conversion = line.strip().split()
                assert np.allclose(float(conversion), shellmound_sfrdata.const)
            if 'end options' in line.lower():
                break

    # check for extra line endings
    # (https://github.com/pandas-dev/pandas/issues/25048)
    cols = ['rno', 'k', 'i', 'j', 'rlen', 'rwid', 'rgrd', 'rtp', 'rbth', 'rhk', 'man', 'ncon', 'ustrf', 'ndv']
    with open(sfr_package_file) as src:
        for line in src:
            if 'begin packagedata' in line.lower():
                next(src) # skip the header
                for line in src:
                    if 'end packagedata' in line.lower():
                        break
                    assert len(line.strip().split()) >= len(cols) -2

    # check the actual values
    pdata = pd.DataFrame(mf6sfr.packagedata.array)
    active = [True if cellid != 'none' else False for cellid in pdata['cellid']]
    pdata = pdata.loc[active]
    pdata['k'], pdata['i'], pdata['j'] = cellids_to_kij(pdata.cellid)
    pdata.drop('cellid', axis=1, inplace=True)
    cols = ['rno', 'ifno', 'k', 'i', 'j', 'rlen', 'rwid', 'rgrd', 'rtp', 'rbth', 'rhk', 'man', 'ncon', 'ustrf', 'ndv']
    cols = [c for c in cols if c in pdata.columns]
    pdata = pdata[cols].copy()
    rows = []
    with open(sfr_package_file) as src:
        for line in src:
            if 'begin packagedata' in line.lower():
                for line in src:
                    if 'end packagedata' in line.lower():
                        break
                    if 'none' not in line.lower() and '#' not in line:
                        rows.append(line)
                break
    text = ''.join(rows)
    df = pd.read_csv(io.StringIO(text), delim_whitespace=True,
                         #skiprows=2,
                         header=None,
                         names=cols,
                         index_col=False
                         )
    isna = df.isna().any(axis=1).values.astype(bool)
    df.dropna(axis=0, inplace=True)
    for c in ['k', 'i', 'j', cols[0]]:
        df[c] = df[c].astype(int) - 1
    pd.testing.assert_frame_equal(df,
                                  pdata.loc[~isna].reset_index(drop=True),
                                  check_dtype=False)


@pytest.mark.parametrize('kwargs', [{'rno': 1},  # specified reach(es)
                                    {'segments': 1},  # specified segment(s)
                                    {'line_ids': 17955471},  # specified line numbers in source hydrography
                                    {},  # all segments
                                    ])
def test_to_riv(shellmound_sfrdata, kwargs):
    sfrd = copy.deepcopy(shellmound_sfrdata)

    riv = sfrd.to_riv(**kwargs)

    riv_spd = riv.stress_period_data
    rd = sfrd.reach_data
    assert 'per' in riv_spd.columns
    if len(kwargs) > 0:
        overlap = set(riv_spd.node).intersection(rd.node)
        assert len(overlap) == 0
    assert sfrmaker.checks.valid_nsegs(riv_spd.rno, riv_spd.outreach)
    nseg, outseg = sfrd.segment_data.nseg, sfrd.segment_data.outseg,
    assert sfrmaker.checks.routing_numbering_is_valid(nseg, outseg,
                                                      rd.iseg, rd.ireach,
                                                      rd.rno, rd.outreach,
                                                      increasing_nseg=True)
    # would be nice to have a way to compare stages
    # but some sfr were removed that didn't get converted to RIV
    # (minor reaches collocated with reaches that got converted)


def test_run_diagnostics(sfrdata):
    """Check that flopy diagnostics were run
    (that a .chk output file was produced)."""
    sfrdata.run_diagnostics()
    checkfile = os.path.join(sfrdata.model.model_ws,
                             '{}_SFR.chk'.format(sfrdata.package_name))
    assert os.path.getsize(checkfile) > 0
    
    
def test_reach_asum(sfrdata):
    # note: asums will only be greater than zero if an arbolate sum field 
    # is supplied as input to Lines class
    rd = sfrdata.reach_data.dropna(subset=['asum'], axis=0)
    assert rd.asum.sum() > 0
    assert np.all(rd.asum >= 0)
