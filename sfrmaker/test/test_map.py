import pandas as pd
import pytest

import sfrmaker


# TODO: make tests more rigorous

@pytest.fixture(scope="module")
def name_path():
    return 'map_test', 'shellmound'


@pytest.fixture(scope="module")
def lines(testdatapath, name_path):
    name, path = name_path
    lns = sfrmaker.Lines.from_shapefile('{}/{}/flowlines.shp'.format(testdatapath, path),
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
    return lns


@pytest.fixture(scope="module")
def grid(testdatapath, name_path):
    name, path = name_path
    grd = sfrmaker.StructuredGrid.from_json('{0}/{1}/{1}/{1}_grid.json'.format(testdatapath, path))
    return grd


@pytest.fixture(scope="module")
def sfr(lines, grid):
    sfr = lines.to_sfr(grid=grid,
                       isfr=None,
                       cull_flowlines_to_active_area=True,
                       one_reach_per_cell=True
                       )
    return sfr


@pytest.fixture(scope="module")
def sfr_with_inflows(sfr):
    # add some inflows
    tmp = sfr.segment_data.loc[sfr.segment_data.nseg.isin([17, 18])].sort_values(by='nseg').copy()
    dfs = []
    for i in range(1, 4):
        itmp = tmp.copy()
        itmp['flow'] = [500 * i, 1000 * i]
        itmp['per'] = i
        dfs.append(itmp)
    sfr.segment_data = sfr.segment_data.append(pd.concat(dfs))
    return sfr


def test_write_mf2005(sfr, outdir, name_path):
    name, path = name_path
    # write an mf2005 version
    sfr.write_package('{}/{}.sfr'.format(outdir, name))


def test_write_mf6(sfr, outdir, name_path):
    name, path = name_path
    # write a MODFLOW6 version
    sfr.write_package('{}/{}.sfr'.format(outdir, name), version='mf6')


def test_shapefile_export(sfr, outdir, name_path):
    name, path = name_path
    sfr.export_cells('{}/{}_cells.shp'.format(outdir, name))
    sfr.export_outlets('{}/{}_outlets.shp'.format(outdir, name))
    sfr.export_transient_variable('flow', '{}/{}_inlets.shp'.format(outdir, name))
    sfr.export_lines('{}/{}_lines.shp'.format(outdir, name))
    sfr.export_routing('{}/{}_lines.shp'.format(outdir, name))
