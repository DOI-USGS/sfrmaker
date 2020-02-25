import io
import os
import numpy as np
import pandas as pd
import pytest
from gisutils import shp2df


@pytest.fixture
def flux_observation_data():
    data = ('datetime,per,line_id,flow_m3d,site_no,junk,comment\n'
            '1998-04-01,0,17991438,100.,7281600,262,RF model estimate\n'
            '1998-05-01,0,17956091,10.,7288570,204,RF model estimate\n'
            '1998-04-01,1,17991438,200.,7281600,262,base flow separation\n'
            '1998-05-01,1,17956091,20.,7288570,204,base flow separation\n'
            )
    data = pd.read_csv(io.StringIO(data))
    return data


@pytest.fixture
def mf6_observations_file(outdir):
    return os.path.join(outdir, 'test.observations_file.obs')


def check_gage_package(fname, name_file, expected):
    assert os.path.exists(fname)
    assert os.path.exists(name_file)
    data = pd.read_csv(fname, header=None, skiprows=1, delim_whitespace=True)
    pd.testing.assert_frame_equal(data, expected, check_dtype=False)
    with open(name_file) as src:
        for line in src:
            if '.ggo' in line.lower():
                _, unit, fname = line.strip().split()
                assert os.path.exists(fname)


def check_mf6_obs_file(fname, expected):
    with open(fname) as src:
        data = []
        for line in src:
            if '# obsname' in line:
                for line in src:
                    if 'end' in line.lower():
                        break
                    data.append(line.strip().split())
        data = pd.DataFrame(data, columns=['obsname', 'obstype', 'rno'])
        data['rno'] = data.rno.astype(int)
        pd.testing.assert_frame_equal(data, expected, check_dtype=False)


def test_write_gage_package(tylerforks_sfrdata, flux_observation_data, outdir):
    obsdata = flux_observation_data.copy()
    obsdata['line_id'] = tylerforks_sfrdata.reach_data.line_id.unique()[:len(obsdata)]

    gage_package_file = os.path.join(tylerforks_sfrdata.model.model_ws,
                                     tylerforks_sfrdata.package_name + '.gage')
    name_file = os.path.join(tylerforks_sfrdata.model.model_ws,
                             tylerforks_sfrdata.model.namefile)
    obs = tylerforks_sfrdata.add_observations(obsdata,
                                              line_id_column_in_data='line_id',
                                              obsname_column_in_data='site_no'
                                              )
    tylerforks_sfrdata.write_gage_package()
    expected = tylerforks_sfrdata.observations.sort_values(by='obsname')[['iseg', 'ireach']].reset_index(drop=True).copy()
    expected['unit'] = np.arange(228, 228+len(obsdata))
    expected.columns = list(range(3))
    expected[3] = 0

    check_gage_package(gage_package_file, name_file,
                       expected=expected)
    os.remove(gage_package_file)
    tylerforks_sfrdata.write_package()
    check_gage_package(gage_package_file, name_file,
                       expected=expected)


def test_write_mf6_sfr_obsfile(shellmound_sfrdata, flux_observation_data, outdir):
    mf6_observations_file = os.path.join(outdir, 'test.observations_file.obs')
    sfr_package_file = os.path.join(outdir, 'test.package_file.sfr')

    shellmound_sfrdata.observations_file = mf6_observations_file
    obs = shellmound_sfrdata.add_observations(flux_observation_data,
                                              obstype='downstream-flow',
                                              line_id_column_in_data='line_id',
                                              obsname_column_in_data='site_no'
                                              )
    shellmound_sfrdata.write_mf6_sfr_obsfile()
    expected = shellmound_sfrdata.observations[['obsname', 'obstype', 'rno']]
    check_mf6_obs_file(shellmound_sfrdata.observations_file,
                       expected=expected)
    os.remove(shellmound_sfrdata.observations_file)
    shellmound_sfrdata.write_package(filename=sfr_package_file, version='mf6')
    check_mf6_obs_file(shellmound_sfrdata.observations_file,
                       expected=expected)
    with open(sfr_package_file) as src:
        for line in src:
            if 'obs6' in line.lower():
                _, _, fname = line.strip().split()
                assert os.path.exists(fname)
                break


def test_add_observations(shellmound_sfrdata, flux_observation_data, outdir):
    obs = shellmound_sfrdata.add_observations(flux_observation_data,
                                              obstype='downstream-flow',
                                              line_id_column_in_data='line_id',
                                              obsname_column_in_data='site_no'
                                              )
    assert np.all(obs == shellmound_sfrdata._observations)
    assert set(obs.columns) == {'obsname', 'obstype', 'rno', 'iseg', 'ireach'}
    rd = shellmound_sfrdata.reach_data.loc[shellmound_sfrdata.reach_data.ireach == 1]
    rno = dict(zip(rd.line_id, rd.rno))
    assert set(obs.rno) == set([rno[lid] for lid in flux_observation_data.line_id])
    rd = shellmound_sfrdata.reach_data
    iseg_ireach = dict(list(zip(rd.rno, zip(rd.iseg, rd.ireach))))
    for i, r in obs.iterrows():
        assert (r.iseg, r.ireach) == iseg_ireach[r.rno]

    out_shapefile = os.path.join(outdir, 'obs.shp')

    # test shapefile export
    shellmound_sfrdata.export_observations(filename=out_shapefile)
    df = shp2df(out_shapefile)
    pd.testing.assert_frame_equal(df.drop('geometry', axis=1),
                                  shellmound_sfrdata.observations,
                                  check_dtype=False
                                  )

    # test assigning obs from custom reach number column?
    obs = shellmound_sfrdata.add_observations(flux_observation_data,
                                              obstype='downstream-flow',
                                              rno_column_in_data='junk',
                                              obsname_column_in_data='site_no'
                                              )
    assert set(obs.rno) == set(flux_observation_data.junk)
