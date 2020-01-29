import flopy
import pandas as pd
from shapely.geometry import box

from sfrmaker.routing import find_path, make_graph
from gisutils import shp2df
from .routing import get_next_id_in_subset
from .utils import load_sr


def get_inflow_locations_from_parent_model(parent_reach_data, inset_reach_data,
                                           inset_sr, active_area=None
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
    inset_sr : flopy.utils.reference.SpatialReference instance describing inset model grid
        Must be in same coordinate system as geometries in parent_reach_data.
        Required only if active_area is None.
    active_area : shapely.geometry.Polygon object
        Describes the area of the inset model where SFR is applied. Used to find
        inset reaches from parent model. Must be in same coordinate system as
        geometries in parent_reach_data. Required only if inset_sr is None.

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
    if isinstance(inset_sr, str):
        sr = load_sr(inset_sr)
    elif isinstance(inset_sr, flopy.discretization.grid.Grid):
        sr = inset_sr
    else:
        raise ValueError('Unrecognized input for inset_sr')

    if active_area is None:
        active_area = box(*sr.bounds)

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


def add_to_perioddata(sfrdata, data, flowline_routing=None,
                      variable='inflow',
                      line_id_column_in_data=None,
                      rno_column_in_data=None,
                      period_column_in_data='per',
                      variable_column_in_data='Q_avg'):
    """Add data to the period data table (sfrdata.period_data)
    for a MODFLOW-6 style sfrpackage. Source data

    """
    sfrd = sfrdata

    # cull data to valid periods
    data = data.loc[data[period_column_in_data] >= 0].copy()

    # map NHDPlus COMIDs to reach numbers
    if flowline_routing is not None:
        assert line_id_column_in_data in data.columns, \
            "Data need an id column so {} locations can be mapped to reach numbers".format(variable)
        rno_column_in_data = 'rno'
        r1 = sfrd.reach_data.loc[sfrd.reach_data.ireach == 1]
        line_id_rno_mapping = dict(zip(r1['line_id'], r1['rno']))
        line_ids = get_next_id_in_subset(r1.line_id, flowline_routing,
                                         data[line_id_column_in_data])
        data[rno_column_in_data] = [line_id_rno_mapping[lid] for lid in line_ids]
    else:
        assert rno_column_in_data in data.columns, \
            "Data to add need reach number or flowline routing information is needed."

    # check for duplicate inflows in same path
    if variable == 'inflow':
        line_ids = set(data[line_id_column_in_data])
        drop = set()
        dropped_line_info_file = 'dropped_inflows_locations.csv'
        for lid in line_ids:
            path = find_path(flowline_routing, start=lid)
            duplicated = set(path[1:]).intersection(line_ids)
            if len(duplicated) > 0:
                drop.add(lid)
                txt = ('warning: {}: {} is upstream '
                       'of the following line_ids:\n{}\n'
                       'see {} for details.').format(line_id_column_in_data,
                                                     lid, duplicated,
                                                     dropped_line_info_file
                                                     )
                print(txt)
        if len(drop) > 0:
            data.loc[data[line_id_column_in_data].isin(drop)].to_csv(dropped_line_info_file, index=False)
            data = data.loc[~data[line_id_column_in_data].isin(drop)]

    # add inflows to period_data
    period_data = sfrd.period_data
    period_data['rno'] = data[rno_column_in_data]
    period_data['per'] = data[period_column_in_data]
    period_data[variable] = data[variable_column_in_data]


def add_to_segmentdata(sfrdata, flowline_routing, data,
                       variable='flow',
                       segment_column_in_data='segment',
                       period_column_in_data='per',
                       variable_column_in_data='Q_avg'):
    """Like add_to_perioddata, but for MODFLOW-2005.
    """
    raise NotImplementedError()
