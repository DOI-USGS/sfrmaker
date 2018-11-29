import os
import pytest
import pandas as pd
import sfrmaker

@pytest.fixture(scope="module", autouse=True)
def create_tmp():
    if not os.path.isdir('tmp'):
        os.makedirs('tmp')


def test_make_sfr():
    name = 'map_test'
    lns = sfrmaker.lines.from_shapefile('data/{}_flowlines.shp'.format(name),
                                        id_column='COMID',
                                        routing_column='tocomid',
                                        width1_column='width1',
                                        width2_column='width2',
                                        up_elevation_column='elevupsmo',
                                        dn_elevation_column='elevdnsmo',
                                        name_column='GNIS_NAME',
                                        attr_length_units='feet',
                                        attr_height_units='feet',
                                        epsg=5070)
    grd = sfrmaker.StructuredGrid.from_json('data/map_test_grid.json')
    sfr = lns.to_sfr(grid=grd,
                              isfr=None,
                              cull_flowlines_to_active_area=True,
                              one_reach_per_cell=True
                              )

    # add some inflows
    tmp = sfr.segment_data.loc[sfr.segment_data.nseg.isin([17, 18])].sort_values(by='nseg').copy()
    dfs = []
    for i in range(1, 4):
        itmp = tmp.copy()
        itmp['flow'] = [500*i, 1000*i]
        itmp['per'] = i
        dfs.append(itmp)
    sfr.segment_data = sfr.segment_data.append(pd.concat(dfs))

    sfr.export_cells('tmp/{}_cells.shp'.format(name))
    sfr.export_outlets('tmp/{}_outlets.shp'.format(name))
    sfr.export_transient_variable('flow', 'tmp/{}_inlets.shp'.format(name))
    sfr.export_lines('tmp/{}_lines.shp'.format(name))
    sfr.export_routing('tmp/{}_lines.shp'.format(name))
