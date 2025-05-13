import os

import numpy as np
import pandas as pd
import geopandas as gpd
import pytest

import sfrmaker
from sfrmaker.checks import reach_elevations_decrease_downstream
from sfrmaker.gis import get_authority_crs
from gisutils import project
from sfrmaker import Lines
from .test_grid import tyler_forks_grid_from_shapefile


@pytest.fixture(scope='module')
def dem(datapath):
    return os.path.join(datapath, 'tylerforks/dem_26715.tif')


@pytest.mark.parametrize('method', ['cell polygons', 'buffers'])
def test_sample_elevations(dem, tylerforks_sfrdata, datapath, method):
    sfr = tylerforks_sfrdata
    sampled_elevs = sfr.sample_reach_elevations(dem, method=method, smooth=True)
    sfr.reach_data['strtop'] = [sampled_elevs[rno] for rno in sfr.reach_data['rno']]
    assert reach_elevations_decrease_downstream(sfr.reach_data)


@pytest.mark.parametrize('fmt,crs,field_data,errors',(
    ('df', 26715, [692396.7, 5140825.5, 392.], 'raise'),
    pytest.param('df', 26715, [[692396.7, 5140825.5, 392.],
                   [0, 0, 400.]], 'raise', marks=pytest.mark.xfail),
    ('df', 26715, [[692396.7, 5140825.5, 392.],
                   [0, 0, 400.]], 'ignore'),
    ('df', 4326, [-90.497579, 46.39570, 392], 'raise'),
    pytest.param('df', 26715, [0, 0, 0], 'raise', marks=pytest.mark.xfail), 
    ('gdf', 26715, [692396.7, 5140825.5, 392.], 'raise'),
    ('gdf', 4326, [-90.497579, 46.395510, 392], 'raise'),
    ('csv', None, [692396.7, 5140825.5, 392.], 'raise'),
    ('shp', 4326, [-90.497579, 46.395510, 392], 'raise'),
    pytest.param('shp', None, [46.395510,-90.497579, 392], 'raise', marks=pytest.mark.xfail), 
    ('gpkg', 26715, [692396.7, 5140825.5, 392.], 'raise'),
))
def test_sample_elevations_with_field_data(dem, tylerforks_sfrdata, field_data, crs, fmt, 
                                           errors, outdir):
    sfr = tylerforks_sfrdata
    
    field_data_df = pd.DataFrame(np.atleast_2d(field_data),
                                 columns=['x', 'y', 'elevation'])
    def get_gdf(field_data_df):
        return gpd.GeoDataFrame(
            field_data_df, 
            geometry=gpd.points_from_xy(x=field_data_df['x'],
                                        y=field_data_df['y']),
            crs=crs)
        
    if fmt == 'df':
        field_data_input = field_data_df
    elif fmt == 'gdf':
        field_data_input = get_gdf(field_data_df)
    else: 
        field_data_input = outdir / f'elevations.{fmt}'
        gdf = get_gdf(field_data_df)
        if fmt == 'csv':
            gdf.drop('geometry', axis=1).to_csv(field_data_input)
        else:
            gdf.to_file(field_data_input)
    
    sampled_elevs = sfr.sample_reach_elevations(dem, smooth=True, 
                                                elevation_data=field_data_input, 
                                                elevation_data_crs=crs,
                                                elevation_data_errors=errors)
    reach_data = sfr.reach_data.copy()
    reach_data.index = reach_data['rno']
    sfr.reach_data['strtop'] = [sampled_elevs[rno] for rno in sfr.reach_data['rno']]
    from sfrmaker.observations import get_closest_reach
    field_data_gdf = get_gdf(field_data_df)
    if crs is not None and crs != sfr.crs:
        field_data_gdf.to_crs(sfr.crs, inplace=True)
        field_data_df['x'] = [g.x for g in field_data_gdf['geometry']]
        field_data_df['y'] = [g.y for g in field_data_gdf['geometry']]
    reaches, distances = get_closest_reach(field_data_df['x'], field_data_df['y'], sfr.reach_data)
    for i, r in field_data_df.iterrows():
        # arbitrary, but avoids trying to compare elevations that were not included
        # when elevation_data_errors='ignore'
        if distances[i] < 100:  
            assert np.allclose(sampled_elevs[reaches[i]], r['elevation'])
    assert reach_elevations_decrease_downstream(sfr.reach_data)

