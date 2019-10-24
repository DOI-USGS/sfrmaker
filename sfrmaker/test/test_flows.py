from ..flows import get_inflow_locations_from_parent_model


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
