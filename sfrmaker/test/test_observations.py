import io
import os
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, box
import pytest
from gisutils import shp2df, df2shp
from sfrmaker.observations import locate_sites, get_closest_reach


@pytest.fixture
def flux_observation_data():
    data = ('datetime,per,line_id,flow_m3d,site_no,junk,x,y,comment\n'
            '1998-04-01,0,17991438,100.,7281600,262,533280,1192740,RF model estimate\n'
            '1998-05-01,0,17956091,10.,7288570,204,515375.2,1189942.5,RF model estimate\n'
            '1998-04-01,1,17991438,200.,7281600,262,533280,1192740,base flow separation\n'
            '1998-05-01,1,17956091,20.,7288570,204,515375.2,1189942.5,base flow separation\n'
            )
    data = pd.read_csv(io.StringIO(data))
    return data


@pytest.fixture
def mf6_observations_file(outdir):
    return os.path.join(outdir, 'test.observations_file.obs')


def check_gage_package(fname, name_file, expected_data, expected_sites):
    assert os.path.exists(fname)
    assert os.path.exists(name_file)
    data = pd.read_csv(fname, header=None, skiprows=1, delim_whitespace=True)
    pd.testing.assert_frame_equal(data, expected_data, check_dtype=False)
    sites = set()
    with open(name_file) as src:
        for line in src:
            if '.ggo' in line.lower():
                _, unit, fname = line.strip().split()
                site, _ = os.path.splitext(fname)
                sites.add(site)
    assert not any(sites.symmetric_difference(expected_sites))


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
                                              line_id_column='line_id',
                                              obsname_column='site_no'
                                              )
    tylerforks_sfrdata.write_gage_package()
    expected = tylerforks_sfrdata.observations.sort_values(by='obsname')[['iseg', 'ireach']].reset_index(drop=True).copy()
    expected['unit'] = np.arange(250, 250+len(obsdata))
    expected.columns = list(range(3))
    expected[3] = 0
    tylerforks_sfrdata.model.write_name_file()
    check_gage_package(gage_package_file, name_file,
                       expected_data=expected, expected_sites=set(obs.obsname))
    os.remove(gage_package_file)
    tylerforks_sfrdata.write_package()
    check_gage_package(gage_package_file, name_file,
                       expected_data=expected, expected_sites=set(obs.obsname))


def test_write_mf6_sfr_obsfile(shellmound_sfrdata, flux_observation_data, outdir):
    mf6_observations_file = os.path.join(outdir, 'test.observations_file.obs')
    sfr_package_file = os.path.join(outdir, 'test.package_file.sfr')

    shellmound_sfrdata.observations_file = mf6_observations_file
    obs = shellmound_sfrdata.add_observations(flux_observation_data,
                                              obstype='downstream-flow',
                                              line_id_column='line_id',
                                              obsname_column='site_no'
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
                assert os.path.exists(os.path.join(outdir, fname))
                break                                              

def test_add_observations_from_line_ids(shellmound_sfrdata, flux_observation_data, outdir):
    obs = shellmound_sfrdata.add_observations(flux_observation_data,
                                              obstype='downstream-flow',
                                              line_id_column='line_id',
                                              obsname_column='site_no'
                                              )
    assert np.all(obs == shellmound_sfrdata._observations)
    assert set(obs.columns) == {'obsname', 'obstype', 'rno', 'iseg', 'ireach'}
    # get the last reach in each segment
    rd = shellmound_sfrdata.reach_data.sort_values(by=['iseg', 'ireach'], axis=0).groupby('iseg').last()
    rno = dict(zip(rd.line_id, rd.rno))
    assert set(obs.rno) == set([rno[str(lid)] for lid in flux_observation_data.line_id])
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

    # test line_id and obsname same column
    obs = shellmound_sfrdata.add_observations(flux_observation_data,
                                              obstype='downstream-flow',
                                              line_id_column='line_id',
                                              obsname_column='line_id'
                                              )
    assert set(obs.obsname) == set([str(i) for i in flux_observation_data.line_id])

    # test rno_column and obsname same column
    obs = shellmound_sfrdata.add_observations(flux_observation_data,
                                              obstype='downstream-flow',
                                              rno_column='junk',
                                              obsname_column='junk'
                                              )
    assert set(obs.obsname) == set([str(i) for i in flux_observation_data.junk])

    
    # test assigning obs from custom reach number column?
    obs = shellmound_sfrdata.add_observations(flux_observation_data,
                                              obstype='downstream-flow',
                                              rno_column='junk',
                                              obsname_column='site_no'
                                              )
    assert set(obs.rno) == set(flux_observation_data.junk)


def test_add_observations_from_points(shellmound_sfrdata, flux_observation_data, outdir):
    inputs = (flux_observation_data, [flux_observation_data,  # test adding from a list
                                      flux_observation_data])
    for input in inputs:
        obs = shellmound_sfrdata.add_observations(input,
                                                  obstype='downstream-flow',
                                                  x_location_column='x',
                                                  y_location_column='y',
                                                  obsname_column='site_no'
                                                  )
        assert np.all(obs == shellmound_sfrdata._observations)
        assert set(obs.columns) == {'obsname', 'obstype', 'rno', 'iseg', 'ireach'}
        expected = dict(zip(inputs[0].site_no.astype(str), inputs[0].junk))
        returned = dict(zip(obs.obsname, obs.rno))
        assert returned == expected
        rd = shellmound_sfrdata.reach_data
        iseg_ireach = dict(list(zip(rd.rno, zip(rd.iseg, rd.ireach))))
        for i, r in obs.iterrows():
            assert (r.iseg, r.ireach) == iseg_ireach[r.rno]


@pytest.mark.parametrize("x, y, expected", ((515459.9, 1189906.1, 202),
                                            (515375.2, 1189942.5, 204)
                                  )
                         )
def test_get_closest_reach(shellmound_sfrdata, x, y, expected, outdir):
    sfrlines = shellmound_sfrdata.reach_data
    rno, dist = get_closest_reach(x, y, sfrlines,
                                  rno_column='rno')
    assert rno == expected
    assert 0 < dist < shellmound_sfrdata.grid.dx


@pytest.mark.parametrize('reach_id_col', (None, 'rno'))
def test_locate_sites(shellmound_sfrdata, reach_id_col, outdir):

    X, Y, rno = zip(*((515459.9, 1189906.1, 202),
                      (515375.2, 1189942.5, 204)))
    df = pd.DataFrame({'geometry': [Point(x, y) for x, y in zip(X, Y)],
                       'site_no': rno
                       }
                      )
    sites_shapefile = '{}/sites.shp'.format(outdir)
    df2shp(df, sites_shapefile, crs=5070)
    sfrlines_shapefile = '{}/shellmound_lines.shp'.format(outdir)
    shellmound_sfrdata.export_lines(sfrlines_shapefile)
    # test reading sfrlines as a dataframe
    # and sfrlines without a reach number column
    if reach_id_col is None:
        reach_id_col = 'rno'
        sfrlines = gpd.read_file(sfrlines_shapefile)
        sfrlines.drop('rno', axis=1, inplace=True)
        sfrlines_shapefile = sfrlines
    active_area = box(*shellmound_sfrdata.grid.bounds)
    locs = locate_sites(sites_shapefile,
                        sfrlines_shapefile,
                        active_area,
                        keep_columns=None,
                        reach_id_col=reach_id_col,
                        ireach_col='ireach',
                        iseg_col='iseg',
                        site_number_col='site_no',
                        perimeter_buffer=1000,
                        distance_threshold=1600
                         )
    assert np.array_equal(locs.rno.values, locs.site_no.values)
    # check that iseg and ireach columns are in the located sites table
    # (for modflow-2005 style sfr packages)
    assert 'iseg' in locs.columns
    assert 'ireach' in locs.columns
