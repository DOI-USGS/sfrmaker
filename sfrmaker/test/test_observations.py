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
        assert data.equals(expected)


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
    check_mf6_obs_file(shellmound_sfrdata.observations_file,
                       expected=shellmound_sfrdata.observations)
    os.remove(shellmound_sfrdata.observations_file)
    shellmound_sfrdata.write_package(filename=sfr_package_file, version='mf6')
    check_mf6_obs_file(shellmound_sfrdata.observations_file,
                       expected=shellmound_sfrdata.observations)
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
    assert set(obs.columns) == {'obsname', 'obstype', 'rno'}
    rd = shellmound_sfrdata.reach_data.loc[shellmound_sfrdata.reach_data.ireach == 1]
    rno = dict(zip(rd.line_id, rd.rno))
    assert set(obs.rno) == set([rno[lid] for lid in flux_observation_data.line_id])
    out_shapefile = os.path.join(outdir, 'obs.shp')

    # test shapefile export
    shellmound_sfrdata.export_observations(filename=out_shapefile)
    df = shp2df(out_shapefile)
    assert df.drop('geometry', axis=1).equals(shellmound_sfrdata.observations)

    # test assigning obs from custom reach number column?
    obs = shellmound_sfrdata.add_observations(flux_observation_data,
                                              obstype='downstream-flow',
                                              rno_column_in_data='junk',
                                              obsname_column_in_data='site_no'
                                              )
    assert set(obs.rno) == set(flux_observation_data.junk)
