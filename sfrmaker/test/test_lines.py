import pytest
import sfrmaker
from sfrmaker.checks import is_to_one
from sfrmaker.nhdplus_utils import load_nhdplus_v2, get_prj_file


@pytest.fixture
def nhdplus_dataframe(datapath):
    pfvaa_files = ['{}/tylerforks/NHDPlus/NHDPlusAttributes/PlusFlowlineVAA.dbf'.format(datapath)]
    plusflow_files = ['{}/tylerforks/NHDPlus/NHDPlusAttributes/PlusFlow.dbf'.format(datapath)]
    elevslope_files = ['{}/tylerforks/NHDPlus/NHDPlusAttributes/elevslope.dbf'.format(datapath)]
    flowlines = ['{}/tylerforks/NHDPlus/NHDSnapshot/Hydrography/NHDFlowline.shp'.format(datapath)]
    df = load_nhdplus_v2(NHDFlowlines=flowlines, PlusFlowlineVAA=pfvaa_files,
                         PlusFlow=plusflow_files, elevslope=elevslope_files,
                         filter='{}/tylerforks/grid.shp'.format(datapath),
                         )
    return df


@pytest.mark.parametrize('one_to_many', (True, False))
def test_original_routing_attribute(nhdplus_dataframe, one_to_many, datapath):
    df = nhdplus_dataframe
    if one_to_many:
        # add some fake distributaries to a few flowlines
        id = df.COMID.values[0]
        inds = df.index.values[:3]
        modified_toids = [[toids[0], id] for toids in df.loc[inds, 'tocomid']] + df.loc[df.index.values[3]:, 'tocomid'].tolist()
        df['tocomid'] = modified_toids
        assert not is_to_one(df.tocomid)
    flowlines = ['{}/tylerforks/NHDPlus/NHDSnapshot/Hydrography/NHDFlowline.shp'.format(datapath)]
    prjfile = get_prj_file(NHDFlowlines=flowlines)

    # convert arbolate sums from km to m
    df['asum2'] = df.ArbolateSu * 1000

    # convert comid end elevations from cm to m
    if 'MAXELEVSMO' in df.columns:
        df['elevup'] = df.MAXELEVSMO / 100.
    if 'MINELEVSMO' in df.columns:
        df['elevdn'] = df.MINELEVSMO / 100.
    lines = sfrmaker.Lines.from_dataframe(df, id_column='COMID',
                                          routing_column='tocomid',
                                          name_column='GNIS_NAME',
                                          attr_length_units='meters',
                                          attr_height_units='meters',
                                          prjfile=prjfile)
    # verify that values (to ids) in original routing dictionary are scalars
    assert is_to_one(lines._original_routing)


def test_write_shapefile(tylerforks_lines_from_NHDPlus, test_data_path):
    lines = tylerforks_lines_from_NHDPlus
    outshp = test_data_path / 'lines.shp'
    lines.write_shapefile(outshp)
    assert outshp.exists()
