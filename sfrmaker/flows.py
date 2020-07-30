import numpy as np
import pandas as pd
from shapely.geometry import box
import flopy
from sfrmaker.routing import find_path, make_graph
from gisutils import shp2df
from mfexport.budget_output import read_sfr_output
from .fileio import read_tables
from .routing import get_next_id_in_subset
from sfrmaker.fileio import load_modelgrid


def get_inflow_locations_from_parent_model(parent_reach_data, inset_reach_data,
                                           inset_grid, active_area=None
                                           ):
    """Get places in an inset model SFR network where the parent SFR network crosses
    the inset model boundary, using common line ID numbers from parent and inset reach datasets.
    MF2005 or MF6 supported; if either dataset contains only reach numbers (is MODFLOW-6),
    the reach numbers are used as segment numbers, with each segment only having one reach.

    Parameters
    ----------
    parent_reach_data : str (filepath) or DataFrame
        SFR reach data for parent model. Must include columns:
        line_id : int; unique identifier for hydrography line that each reach is based on
        rno : int; unique identifier for each reach. Optional if iseg and ireach columns are included.
        iseg : int; unique identifier for each segment. Optional if rno is included.
        ireach : int; unique identifier for each reach. Optional if rno is included.
        geometry : shapely.geometry object representing location of each reach
    inset_reach_data : str (filepath) or DataFrame
        SFR reach data for inset model. Same columns as parent_reach_data,
        except a geometry column isn't needed. line_id values must correspond to
        same source hydrography as those in parent_reach_data.
    inset_grid : flopy.discretization.StructuredGrid instance describing model grid
        Must be in same coordinate system as geometries in parent_reach_data.
        Required only if active_area is None.
    active_area : shapely.geometry.Polygon object
        Describes the area of the inset model where SFR is applied. Used to find
        inset reaches from parent model. Must be in same coordinate system as
        geometries in parent_reach_data. Required only if inset_grid is None.

    Returns
    -------
    locations : DataFrame
        Columns:
        parent_segment : parent model segment
        parent_reach : parent model reach
        parent_rno : parent model reach number
        line_id : unique identifier for hydrography line that each reach is based on
    """

    # spatial reference instances defining parent and inset grids
    if isinstance(inset_grid, str):
        grid = load_modelgrid(inset_grid)
    elif isinstance(inset_grid, flopy.discretization.grid.Grid):
        grid = inset_grid
    else:
        raise ValueError('Unrecognized input for inset_grid')

    if active_area is None:
        l, r, b, t = grid.extent
        active_area = box(l, b, r, t)

    # parent and inset reach data
    if isinstance(parent_reach_data, str):
        prd = shp2df(parent_reach_data)
    elif isinstance(parent_reach_data, pd.DataFrame):
        prd = parent_reach_data.copy()
    else:
        raise ValueError('Unrecognized input for parent_reach_data')
    if 'rno' in prd.columns and 'iseg' not in prd.columns:
        prd['iseg'] = prd['rno']
        prd['ireach'] = 1
    mustinclude_cols = {'line_id', 'rno', 'iseg', 'ireach', 'geometry'}
    assert len(mustinclude_cols.intersection(prd.columns)) == len(mustinclude_cols)

    if isinstance(inset_reach_data, str):
        if inset_reach_data.endswith('.shp'):
            ird = shp2df(inset_reach_data)
        else:
            ird = pd.read_csv(inset_reach_data)
    elif isinstance(inset_reach_data, pd.DataFrame):
        ird = inset_reach_data.copy()
    else:
        raise ValueError('Unrecognized input for inset_reach_data')
    if 'rno' in ird.columns and 'iseg' not in ird.columns:
        ird['iseg'] = ird['rno']
        ird['ireach'] = 1
    mustinclude_cols = {'line_id', 'rno', 'iseg', 'ireach'}
    assert len(mustinclude_cols.intersection(ird.columns)) == len(mustinclude_cols)

    graph = make_graph(ird.rno.values, ird.outreach.values, one_to_many=False)

    # cull parent reach data to only lines that cross or are just upstream of inset boundary
    buffered = active_area.buffer(5000, cap_style=2)
    close = [g.intersects(buffered) for g in prd.geometry]
    prd = prd.loc[close]
    prd.index = prd.rno
    boundary = active_area.exterior
    inset_line_id_connections = {}  # parent rno: inset line_id
    for i, r in prd.iterrows():
        if r.outreach not in prd.index:
            continue
        downstream_line = prd.loc[r.outreach, 'geometry']
        upstream_line = prd.loc[prd.rno == r.outreach, 'geometry'].values[0]
        intersects = r.geometry.intersects(boundary)
        intersects_downstream = downstream_line.within(active_area)
        # intersects_upstream = upstream_line.within(active_area)
        in_inset_model = r.geometry.within(active_area)
        if intersects_downstream:
            if intersects:
                # if not intersects_upstream: # exclude lines that originated within the model
                #    # lines that cross route to their counterpart in inset model
                inset_line_id_connections[r.rno] = r.line_id
                pass
            elif not in_inset_model:
                # lines that route to a line within the inset model
                # route to that line's inset counterpart
                inset_line_id_connections[r.rno] = prd.loc[r.outreach, 'line_id']
                pass

    prd = prd.loc[prd.rno.isin(inset_line_id_connections.keys())]

    # parent rno lookup
    parent_rno_lookup = {v: k for k, v in inset_line_id_connections.items()}

    # inlet reaches in inset model
    ird = ird.loc[ird.ireach == 1]
    ird = ird.loc[ird.line_id.isin(inset_line_id_connections.values())]

    # for each reach in ird (potential inset inlets)
    # check that there isn't another inlet downstream
    drop_reaches = []
    for i, r in ird.iterrows():
        path = find_path(graph, r.rno)
        another_inlet_downstream = len(set(path[1:]).intersection(set(ird.rno))) > 0
        if another_inlet_downstream:
            drop_reaches.append(r.rno)

    ird = ird.loc[~ird.rno.isin(drop_reaches)]
    # cull parent flows to outlet reaches
    iseg_ireach = zip(prd.iseg, prd.ireach)
    parent_outlet_iseg_ireach = dict(zip(prd.rno, iseg_ireach))

    df = ird[['line_id', 'name', 'rno', 'iseg', 'ireach']].copy()
    df['parent_rno'] = [parent_rno_lookup[lid] for lid in df['line_id']]
    df['parent_iseg'] = [parent_outlet_iseg_ireach[rno][0] for rno in df['parent_rno']]
    df['parent_ireach'] = [parent_outlet_iseg_ireach[rno][1] for rno in df['parent_rno']]
    return df.reset_index(drop=True)


def get_inflows_from_parent_model(parent_reach_data, inset_reach_data,
                                  mf2005_parent_sfr_outputfile, mf6_parent_sfr_budget_file,
                                  inset_grid, active_area=None):
    """Get places in an inset model SFR network where the parent SFR network crosses
    the inset model boundary, using common line ID numbers from parent and inset reach datasets.
    MF2005 or MF6 supported; if either dataset contains only reach numbers (is MODFLOW-6),
    the reach numbers are used as segment numbers, with each segment only having one reach.

    Parameters
    ----------
    parent_reach_data : str (filepath) or DataFrame
        SFR reach data for parent model. Must include columns:
        line_id : int; unique identifier for hydrography line that each reach is based on
        rno : int; unique identifier for each reach. Optional if iseg and ireach columns are included.
        iseg : int; unique identifier for each segment. Optional if rno is included.
        ireach : int; unique identifier for each reach. Optional if rno is included.
        geometry : shapely.geometry object representing location of each reach
    inset_reach_data : str (filepath) or DataFrame
        SFR reach data for inset model. Same columns as parent_reach_data,
        except a geometry column isn't needed. line_id values must correspond to
        same source hydrography as those in parent_reach_data.
    mf2005_parent_sfr_outputfile : str (filepath)
        Modflow-2005 style SFR text file budget output.
    mf6_parent_sfr_budget_file : str (filepath)
        Modflow-6 style SFR binary budget output
    inset_grid : flopy.discretization.StructuredGrid instance describing model grid
        Must be in same coordinate system as geometries in parent_reach_data.
        Required only if active_area is None.
    active_area : shapely.geometry.Polygon object
        Describes the area of the inset model where SFR is applied. Used to find
        inset reaches from parent model. Must be in same coordinate system as
        geometries in parent_reach_data. Required only if inset_grid is None.

    Returns
    -------
    inflows : DataFrame
        Columns:
        parent_segment : parent model segment
        parent_reach : parent model reach
        parent_rno : parent model reach number
        line_id : unique identifier for hydrography line that each reach is based on
    """
    locations = get_inflow_locations_from_parent_model(parent_reach_data=parent_reach_data,
                                                       inset_reach_data=inset_reach_data,
                                                       inset_grid=inset_grid,
                                                       active_area=active_area)
    df = read_sfr_output(mf2005_sfr_outputfile=mf2005_parent_sfr_outputfile,
                         mf6_sfr_stage_file=None,
                         mf6_sfr_budget_file=mf6_parent_sfr_budget_file,
                         model=None)
    j=2


def add_to_perioddata(sfrdata, data, flowline_routing=None,
                      variable='inflow',
                      line_id_column=None,
                      rno_column=None,
                      period_column='per',
                      data_column='Q_avg'):
    """Add data to the period data table (sfrdata.period_data)
    for a MODFLOW-6 style sfrpackage.

    Parameters
    ----------
    sfrdata : sfrmaker.SFRData instance
        SFRData instance with reach_data table attribute. To add observations from x, y coordinates,
        the reach_data table must have a geometry column with LineStrings representing each reach, or
        an sfrlines_shapefile is required. Reach numbers are assumed to be in an 'rno' column.
    data : DataFrame, path to csv file, or list of DataFrames or file paths
        Table with information on the observation sites to be located. Must have
        either reach numbers (rno_column), line_ids (line_id_column),
        or x and y locations (x_column_in_data and y_column_in_data).
    flowline_routing : dict
        Optional dictionary of routing for source hydrography. Only needed
        if locating by line_id, and SFR network is a subset of the full source
        hydrography (i.e. some lines were dropped in the creation of the SFR packge,
        or if the sites are inflow points corresponding to lines outside of the model perimeter).
        In this case, observation points referenced to line_ids that are missing from the SFR
        network are placed at the first reach corresponding to the next downstream line_id
        that is represented in the SFR network. By default, None.
    variable : str, optional
        Modflow-6 period variable (see Modflow-6 Description of Input and Outpu), by default 'inflow'
    line_id_column : str
        Column in data matching observation sites to line_ids in the source hydrography data.
        Either line_id_column or rno_column must be specified. By default, None
    rno_column : str
        Column in data matching observation sites to reach numbers in the SFR network. By default, None.
    period_column : str, optional
        Column with modflow stress period for each inflow value, by default 'per', by default, 'per'.
    data_column : str, optional
        Column with flow values, by default 'Q_avg'

    Returns
    -------
    Updates the sfrdata.perioddata DataFrame.
    """
    sfrd = sfrdata

    # allow input via a list of tables or single table
    data = read_tables(data)

    # cull data to valid periods
    data = data.loc[data[period_column] >= 0].copy()

    # map NHDPlus COMIDs to reach numbers
    if flowline_routing is not None:
        assert line_id_column in data.columns, \
            "Data need an id column so {} locations can be mapped to reach numbers".format(variable)
        # replace ids that are not keys (those outside the network) with zeros
        # (0 is the exit condition for finding paths in get_next_id_in_subset)
        flowline_routing = {k: v if v in flowline_routing.keys() else 0 for k, v in flowline_routing.items()}
        rno_column = 'rno'
        r1 = sfrd.reach_data.loc[sfrd.reach_data.ireach == 1]
        line_id_rno_mapping = dict(zip(r1['line_id'], r1['rno']))
        line_ids = get_next_id_in_subset(r1.line_id, flowline_routing,
                                         data[line_id_column])
        data['line_id_in_model'] = line_ids
        data[rno_column] = [line_id_rno_mapping[lid] for lid in line_ids]
    else:
        assert rno_column in data.columns, \
            "Data to add need reach number or flowline routing information is needed."

    # check for duplicate inflows in same path
    if variable == 'inflow':
        line_ids = set(data[line_id_column])
        drop = set()
        dropped_line_info_file = 'dropped_inflows_locations.csv'
        for lid in line_ids:
            path = find_path(flowline_routing, start=lid)
            duplicated = set(path[1:]).intersection(line_ids)
            if len(duplicated) > 0:
                drop.add(lid)
                txt = ('warning: {}: {} is upstream '
                       'of the following line_ids:\n{}\n'
                       'see {} for details.').format(line_id_column,
                                                     lid, duplicated,
                                                     dropped_line_info_file
                                                     )
                print(txt)
        if len(drop) > 0:
            data.loc[data[line_id_column].isin(drop)].to_csv(dropped_line_info_file, index=False)
            data = data.loc[~data[line_id_column].isin(drop)]

    # add inflows to period_data
    period_data = sfrd.period_data
    period_data['rno'] = data[rno_column]
    period_data['per'] = data[period_column]
    period_data[variable] = data[data_column]
    period_data['specified_line_id'] = data[line_id_column]
    other_columns = [c for c in data.columns
                     if c not in {rno_column, period_column, data_column, line_id_column}]
    for c in other_columns:
        period_data[c] = data[c]


def add_to_segment_data(sfrdata, data, flowline_routing=None,
                        variable='flow',
                        line_id_column=None,
                        segment_column=None,
                        period_column='per',
                        data_column='Q_avg'):
    """Like add_to_perioddata, but for MODFLOW-2005.
    """
    sfrd = sfrdata

    # cull data to valid periods
    data = data.loc[data[period_column] >= 0].copy()

    # map NHDPlus COMIDs to reach numbers
    if flowline_routing is not None:
        assert line_id_column in data.columns, \
            "Data need an id column so {} locations can be mapped to reach numbers".format(variable)
        flowline_routing = {k: v if v in flowline_routing.keys() else 0 for k, v in flowline_routing.items()}
        segment_column = 'segment'
        r1 = sfrd.reach_data.loc[sfrd.reach_data.ireach == 1]
        line_id_iseg_mapping = dict(zip(r1['line_id'], r1['iseg']))
        line_ids = get_next_id_in_subset(r1.line_id, flowline_routing,
                                         data[line_id_column])
        data[segment_column] = [line_id_iseg_mapping[lid] for lid in line_ids]
    else:
        assert segment_column in data.columns, \
            "Data to add need segment number or flowline routing information is needed."

    # check for duplicate inflows in same path
    if variable == 'flow':
        line_ids = set(data[line_id_column])
        drop = set()
        dropped_line_info_file = 'dropped_inflows_locations.csv'
        for lid in line_ids:
            path = find_path(flowline_routing, start=lid)
            duplicated = set(path[1:]).intersection(line_ids)
            if len(duplicated) > 0:
                drop.add(lid)
                txt = ('warning: {}: {} is upstream '
                       'of the following line_ids:\n{}\n'
                       'see {} for details.').format(line_id_column,
                                                     lid, duplicated,
                                                     dropped_line_info_file
                                                     )
                print(txt)
        if len(drop) > 0:
            data.loc[data[line_id_column].isin(drop)].to_csv(dropped_line_info_file, index=False)
            data = data.loc[~data[line_id_column].isin(drop)]

    # rename columns in data to be added to same names as segment_data
    data.rename(columns={period_column: 'per',
                         segment_column: 'nseg',
                         data_column: variable},
                inplace=True)
    # update existing segment data
    sfrd.segment_data.index = pd.MultiIndex.from_tuples(zip(sfrd.segment_data.per, sfrd.segment_data.nseg),
                                                        names=['per', 'nseg'])
    loc = list(zip(data.per, data.nseg))
    data.index = pd.MultiIndex.from_tuples(loc, names=['per', 'nseg'])
    replace = sorted(list(set(data.index).intersection(sfrd.segment_data.index)))
    add = sorted(list(set(data.index).difference(sfrd.segment_data.index)))
    sfrd.segment_data.loc[replace, variable] = data.loc[replace, variable]

    # concat on the added data (create additional rows in segment_data table)
    to_concat = [sfrd.segment_data]
    period_groups = data.loc[add, ['per', 'nseg', variable]].reset_index(drop=True).groupby('per')
    for per, group in period_groups:
        # start with existing data (row) for that segment
        df = sfrd.segment_data.loc[(slice(None, None), group.nseg), :].copy()
        df['per'] = per
        df.index = zip(df.per, df.nseg)
        df[variable] = group[variable].values
        to_concat.append(df)
    sfrd.segment_data = pd.concat(to_concat).reset_index(drop=True)
