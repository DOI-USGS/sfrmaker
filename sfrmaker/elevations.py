"""Methods related to sampling and smoothing elevations."""
import time
from pathlib import Path
import time
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterstats import zonal_stats
from shapely.geometry import LineString, Point
from gisutils import get_authority_crs, project
from sfrmaker.routing import get_nextupsegs, get_upsegs, make_graph


def get_slopes(streambed_tops, reach_lengths, reach_numbers, outreach_numbers, 
               default_slope=0.001, minimum_slope=0.0001,
                maximum_slope=1.):
    """Compute slopes by reach using values in strtop (streambed top) and rchlen (reach length)
    columns of reach_data. The slope for a reach n is computed as strtop(n) - strtop(n+1) / rchlen(n).
    Slopes for outlet reaches are set equal to a default value (default_slope).
    Populates the slope column in reach_data.

    Parameters
    ----------
    streambed_top : sequence
        Streambed top elevations
    reach_lengths : sequence
        Reach lengths
    reach_numbers : sequence
        Unique identifiers for each reach.
    outreach_numbers : sequence
        Unique identifier of next downtream reach.
    default_slope : float
        Slope value applied to outlet reaches (where water leaves the model).
        Default value is 0.001
    minimum_slope : float
        Assigned to reaches with computed slopes less than this value.
        This ensures that the Manning's equation won't produce unreasonable values of stage
        (in other words, that stage is consistent with assumption that
        streamflow is primarily drive by the streambed gradient).
        Default value is 0.0001.
    maximum_slope : float
        Assigned to reaches with computed slopes more than this value.
        Default value is 1.
    """
    # cast everything to lists to avoid confusion with numpy vs. pandas indexers
    streambed_tops = list(streambed_tops)
    reach_lengths = list(reach_lengths)
    reach_numbers = list(reach_numbers)
    outreach_numbers = list(outreach_numbers)
    assert np.sum(outreach_numbers) > 0, \
        ("outreach_numbers appear to be invalid; make sure outreaches are popluated, "
         "for example by running SFRData.set_outreaches()")
    elev = dict(zip(reach_numbers, streambed_tops))
    dist = dict(zip(reach_numbers, reach_lengths))
    dnelev = {rno: elev[outreach_numbers[i]] if outreach_numbers[i] != 0
              else -9999 for i, rno in enumerate(reach_numbers)}
    slopes = np.array(
        [(elev[rno] - dnelev[rno]) / dist[rno] if dnelev[rno] != -9999 and dist[rno] > 0
            else default_slope for rno in reach_numbers])
    slopes[slopes < minimum_slope] = minimum_slope
    slopes[slopes > maximum_slope] = maximum_slope
    return slopes
        
        
def smooth_elevations(fromids, toids, elevations, start_elevations=None):  # elevup, elevdn):
    """

    Parameters
    ----------
    fromids : sequence of hashables
    toids : sequence of hashables
        Downstream connections of fromids
    elevations : sequence of floats
        Elevation for each edge (line) in a stream network, or if start_elevations
        are specified, the end elevation for each edge.
    start_elevations : sequence of floats, optional
        Start elevation for edge (line) in a stream network.
        By default, None.

    Returns
    -------
    Elevations : dict or tuple
        Dictionary of smoothed edge elevations,
        or smoothed end elevations, start elevations
    """
    # make forward and reverse dictionaries with routing info
    graph = dict(zip(fromids, toids))
    assert 0 in set(graph.values()), 'No outlets in routing network!'
    graph_r = make_graph(toids, fromids)

    # make dictionaries of segment end elevations
    elevations = dict(zip(fromids, elevations))
    if start_elevations is not None:
        elevmax = dict(zip(fromids, start_elevations))

    def get_upseg_levels(seg):
        """Traverse routing network, returning a list of segments
        at each level upstream from the outlets. (level 0 route to seg;
        segments in level 1 route to a segment in level 0, etc.)

        Parameters:
        -----------
        seg : int
            Starting segment number

        Returns
        -------
        all_upsegs : list
            List with list of segments at each level
        """
        upsegs = graph_r[seg].copy()
        all_upsegs = [upsegs]
        for i in range(len(fromids)):
            upsegs = get_nextupsegs(graph_r, upsegs)
            if len(upsegs) > 0:
                all_upsegs.append(upsegs)
            else:
                break
        return all_upsegs

    def reset_elevations(seg):
        """Reset segment elevations above (upsegs) and below (outseg) a node.
        """
        oseg = graph[seg]
        all_upsegs = np.array(list(get_upsegs(graph_r, seg)) + [seg])  # all segments upstream of node
        elevmin_s = np.min([elevations[s] for s in all_upsegs])  # minimum current elevation upstream of node
        oldmin_s = elevations[seg]
        elevs = [elevmin_s, oldmin_s]
        if oseg > 0:  # if segment is not an outlet,
            if start_elevations is not None:
                elevs.append(elevmax[oseg]) # outseg start elevation (already updated)
        # set segment end elevation as min of
        # upstream elevations, current elevation, outseg start elevation
        elevations[seg] = np.min(elevs)
        # if the node is not an outlet, reset the outseg max if the current min is lower
        if oseg > 0:
            if start_elevations is not None:
                next_reach_elev = elevmax[oseg]
                elevmax[graph[seg]] = np.min([elevmin_s, next_reach_elev])
            else:
                next_reach_elev = elevations[oseg]
                elevations[graph[seg]] = np.min([elevmin_s, next_reach_elev])

    print('\nSmoothing elevations...')
    ta = time.time()
    # get list of segments at each level, starting with 0 (outlet)
    segment_levels = get_upseg_levels(0)
    # at each level, reset all of the segment elevations as necessary
    for level in segment_levels:
        for s in level:
            if 0 in level:
                j=2
            reset_elevations(s)
    print("finished in {:.2f}s".format(time.time() - ta))
    if start_elevations is not None:
        return elevations, elevmax
    return elevations


def sample_reach_elevations(sfr_reach_data, dem,
                            sfrmaker_modelgrid=None,
                            model_crs=None,
                            method='buffers',
                            buffer_distance=100,
                            smooth=True,
                            elevation_data=None,
                            elevation_data_crs=None,
                            elevation_data_layer=None,
                            elevation_data_errors='raise',
                            ):
    """Computes zonal statistics on a raster for SFR reaches, using
    either buffer polygons around the reach LineStrings, or the model
    grid cell polygons containing each reach.

    Parameters
    ----------
    sfr_reach_data : pd.DataFrame
        Table of SFR reach information, equivalent to the SFRData.reach_data attribute.
    dem : path to valid raster dataset
        Must be in same Coordinate Reference System as model grid.
    sfrmaker_modelgrid : sfrmaker.grid class instance
    model_crs : obj
        Coordinate reference system for the model. 
        Only needed if sfrmaker_modelgrid does not have an valid `crs` attribute.
        A Python int, dict, str, or :class:`pyproj.crs.CRS` instance
        passed to :meth:`pyproj.crs.CRS.from_user_input`

        Can be any of:
        - PROJ string
        - Dictionary of PROJ parameters
        - PROJ keyword arguments for parameters
        - JSON string with PROJ parameters
        - CRS WKT string
        - An authority string [i.e. 'epsg:4326']
        - An EPSG integer code [i.e. 4326]
        - A tuple of ("auth_name": "auth_code") [i.e ('epsg', '4326')]
        - An object with a `to_wkt` method.
        - A :class:`pyproj.crs.CRS` class 
    method : str; 'buffers' or 'cell polygons'
        If 'buffers', buffers (with flat caps; cap_style=2 in LineString.buffer())
        will be created around the reach LineStrings (geometry column in reach_data).
    buffer_distance : float
        Buffer distance to apply around reach LineStrings, in crs_units.
    smooth : bool
        Run sfrmaker.elevations.smooth_elevations on sampled elevations
        to ensure that they decrease monotonically in the downstream direction
        (default=True).
    elevation_data : str/pathlike or DataFrame
        Optional table of streambed top elevations at x, y locations. These will be mapped to
        individual reaches using the selected method, replacing the DEM-derived elevations.
        The points must fall within the selected buffer distance or cell polygon, otherwise
        an error will be raised. Points can be specified in a shapefile, geopackage, CSV file or
        (Geo)DataFrame, which must contain the following fields:
        
        ========= ===============================================================================
        x         x-coordinate. If a CSV or regular DataFrame is provided, 
                    this must be in the CRS of the model or an elevation_data_crs must be provided.
        y         y-coordinate.
        elevation Elevation in the model units.
        ========= ================================================================================
        
        By default, None.
    elevation_data_layer : None
        Layer name for the elevation data, if the data are being supplied in a multi-layer GeoPackage.
        by default, None.
    elevation_data_crs : obj
        Coordinate reference system of the x,y points in elevation-data. 
        Only needed if the data do not have a valid ESRI projection (.prj) file.
        A Python int, dict, str, or :class:`pyproj.crs.CRS` instance
        passed to :meth:`pyproj.crs.CRS.from_user_input`

        Can be any of:
        - PROJ string
        - Dictionary of PROJ parameters
        - PROJ keyword arguments for parameters
        - JSON string with PROJ parameters
        - CRS WKT string
        - An authority string [i.e. 'epsg:4326']
        - An EPSG integer code [i.e. 4326]
        - A tuple of ("auth_name": "auth_code") [i.e ('epsg', '4326')]
        - An object with a `to_wkt` method.
        - A :class:`pyproj.crs.CRS` class
    elevation_data_errors : {‘ignore’, ‘raise’}
        If ‘ignore’, suppress error when one or more points in elevation data is
        not associated with a stream reach (is out of geographic proximity).
        By default, ‘raise’.
    
    Returns
    -------
    streambed_tops : dict of sampled elevations keyed by reach number
    """
    reach_data = sfr_reach_data.copy()
    reach_data.index = reach_data['rno']

    if model_crs is None:
        if sfrmaker_modelgrid is not None:
            model_crs = sfrmaker_modelgrid.crs
        else:
            raise ValueError("sample_reach_elevations requires a valid model_crs "
                             "or .crs attribute attached to sfrmaker_modelgrid")
        
    # get the CRS and pixel size for the DEM
    with rasterio.open(dem) as src:
        raster_crs = get_authority_crs(src.crs)

        # make sure buffer is large enough for DEM pixel size
        buffer_distance = np.max([np.sqrt(src.res[0] *
                                            src.res[1]) * 1.01,
                                    buffer_distance])

    if method == 'buffers':
        assert isinstance(reach_data.geometry.values[0], LineString), \
            "Need LineString geometries in reach_data.geometry column to use buffer option."
        features = [g.buffer(buffer_distance, join_style='mitre', cap_style='flat') 
                    for g in reach_data.geometry]
        feature_descr = 'buffered LineStrings'
    elif method == 'cell polygons':
        assert sfrmaker_modelgrid is not None, \
            "Need an attached sfrmaker.Grid instance to use cell polygons option."
        # TODO: refactor this to use a Flopy modelgrid
        features = sfrmaker_modelgrid.df.loc[reach_data.node, 'geometry'].tolist()
        feature_descr = method

    # to_crs features if they're not in the same crs
    if raster_crs != model_crs:
        features = project(features,
                            model_crs,
                            raster_crs)

    print(f'running rasterstats.zonal_stats on {feature_descr}...')
    t0 = time.time()
    def get_min(x):
        return np.min(x[x > -1e38])
    results = zonal_stats(features,
                            dem,
                            add_stats={'nanmin': get_min},
                            all_touched=True
                            )
    sampled_elevations = [r['nanmin'] for r in results]
    print("finished in {:.2f}s\n".format(time.time() - t0))

    if all(v is None for v in sampled_elevations):
        raise Exception(f'No {feature_descr} intersected with {dem}. Check projections.')
    if any(v is None for v in sampled_elevations):
        raise Exception(f'Some {feature_descr} not intersected with {dem}. '
                        'Check that DEM covers the area of the stream network.'
                        )
    if elevation_data is not None:
        if isinstance(elevation_data, str) or isinstance(elevation_data, Path):
            measured_elevations = gpd.read_file(elevation_data, layer=elevation_data_layer)
        elif isinstance(elevation_data, pd.DataFrame):
            measured_elevations = elevation_data.copy()
        else:
            raise ValueError(f'Unrecognized input for elevation_data:\n{elevation_data}')
        
        if isinstance(measured_elevations, gpd.GeoDataFrame) and\
            getattr(measured_elevations, 'geometry').is_valid.all():
            #x = [g.x for g in measured_elevations.geometry]
            #y = [g.y for g in measured_elevations.geometry]
            pass
        elif {'x', 'y'}.difference(measured_elevations.columns):
            raise ValueError(f'elevation_data input {elevation_data} needs to be '
                                'a shapefile, geopackage or GeoDataFrame with valid point features,\n'
                                "or 'x' and 'y' columns are needed.")
        else:
            x_coords = measured_elevations['x'].astype(float)
            y_coords = measured_elevations['y'].astype(float)
            measured_elevations['geometry'] = [Point(x, y) for x, y in zip(x_coords, y_coords)]
            measured_elevations = gpd.GeoDataFrame(measured_elevations, crs=elevation_data_crs)
        measured_elevations['elevation'] = measured_elevations['elevation'].astype(float)
        
        if measured_elevations.crs is None:
            if elevation_data_crs is not None:
                measured_elevations.crs = elevation_data_crs
            else:
                measured_elevations.crs = model_crs
        if measured_elevations.crs != model_crs:
            measured_elevations.to_crs(model_crs, inplace=True)
            
        # Get stream buffers or model cell polygons containing each point
        sfr_reach_buffers = gpd.GeoDataFrame({
            'rno': reach_data['rno'].values},
            geometry=features, crs=raster_crs)
        sfr_reach_buffers.to_crs(model_crs, inplace=True)
        # avoid potential geopandas error
        measured_elevations.drop('index_right', axis=1, errors='ignore', inplace=True)
        joined = sfr_reach_buffers.sjoin(measured_elevations, how='right')
        not_joined = ~measured_elevations['geometry'].isin(joined.dropna(subset='rno')['geometry'])
        if not_joined.any():
            if isinstance(elevation_data, str) or isinstance(elevation_data, Path):
                outpath = Path(elevation_data).parent
            else:
                outpath = Path('.')
            outfile = outpath / 'dropped_elevation_points.gpkg'
            measured_elevations.loc[not_joined].to_file(outfile, index=False)
            print(f"wrote {outfile}")
            message = ("The following measured elevations did not intersect an SFR cell, "
                       f"or were not within the buffer distance of {buffer_distance}:\n"
                       f"{measured_elevations.loc[not_joined]}\n"
                       f"see {outfile}")
            if elevation_data_errors == 'raise':
                raise ValueError(message)
            else:
                print(message)
                joined.dropna(subset='rno', axis=0, inplace=True)
        # consolidate multiple elevations within an SFR reach by taking the minimum
        joined['rno'] = joined['rno'].astype(int)
        joined_consolidated = joined.groupby('rno').first()
        joined_consolidated['elevation'] = joined.groupby('rno')['elevation'].min()
        
        sampled_elevations = pd.Series(sampled_elevations, index=reach_data['rno'])
        sampled_elevations.update(joined_consolidated['elevation'])
        sampled_elevations = sampled_elevations.tolist()
            
    if smooth:
        streambed_tops = smooth_elevations(reach_data.rno.tolist(),
                                    reach_data.outreach.tolist(),
                                    sampled_elevations)
    else:
        streambed_tops = dict(zip(reach_data.rno, sampled_elevations))
    return streambed_tops