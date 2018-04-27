import sys
sys.path.append('/Users/aleaf/Documents/GitHub/flopy3')
import numpy as np
import pandas as pd
from shapely.geometry import LineString
import flopy
from .utils import renumber_segments, find_path, make_graph
from .checks import valid_rnos, valid_nsegs, rno_nseg_routing_consistent
from .gis import df2shp

fm = flopy.modflow

class sfrdata:
    """Class for working with a streamflow routing (SFR) dataset,
    where the stream network is discretized into reaches contained
    within individual model cells. Reaches may be grouped into segments,
    with routing between segments specified, and routing between reaches
    within segments based on consecutive numbering (as in MODFLOW-2005).
    In this case, unique identifier numbers will be assigned to each reach
    (in the rno column of the reach_data table), and routing connections
    between rnos will be computed. Alternatively, reaches and their
    routing connections can be specified directly, as in MODFLOW-6. In this
    case, MODFLOW-2005 input will be written with one reach per segment.

    Parameters
    ----------
    reach_data : DataFrame
        Table containing information on the SFR reaches.
    segment_data : DataFrame
        Table containing information on the segments (optional).
    grid : sfrmaker.grid class instance
    model_length_units : str
        'm' or 'ft'
    model_time_units : str
        's': seconds
        'm': minutes
        'h': hours
        'd': days
        'y': years
    enforce_increasing_nsegs : bool
        If True, segment numbering is checked to ensure
        that it only increases downstream, and reset if it doesnt.
    model_name : str
            Base name for writing sfr output.
    kwargs : keyword arguments
        Optional values to assign globally to SFR variables. For example
        icalc=1 would assign all segments an icalc value of 1. For default
        values see the sfrdata.defaults dictionary. Default values can be
        assigned using MODFLOW-2005 or MODFLOW-6 terminology.
    """

    # conversions to MODFLOW6 variable names
    mf6names = {'reachID': 'rno',
                'node': 'cellid',
                'rchlen': 'rlen',
                'width': 'rwid',
                'slope': 'rgrd',
                'strtop': 'rtp',
                'strthick': 'rbth',
                'strhc1': 'rhk',
                'roughch': 'manning',
                'flow': 'inflow',
                'pptsw': 'rainfall',
                'etsw': 'evaporation',
                'runoff': 'runoff',
                'depth1': 'depth1',
                'depth2': 'depth2'}

    mf5names = {v:k for k, v in mf6names.items()}

    # order for columns in reach_data
    rdcols = ['rno', 'node', 'k', 'i', 'j',
              'iseg', 'ireach', 'rchlen', 'width', 'slope',
              'strtop', 'strthick', 'strhc1',
              'thts', 'thti', 'eps', 'uhc',
              'outreach', 'outseg', 'asum', 'line_id', 'name',
              'geometry']

    # order for columns in segment_data
    sdcols = ['per', 'nseg', 'icalc', 'outseg', 'iupseg',
              'iprior', 'nstrpts',
              'flow', 'runoff', 'etsw', 'pptsw',
              'roughch', 'roughbk', 'cdpth', 'fdpth',
              'awdth', 'bwdth',
              'hcond1', 'thickm1', 'elevup', 'width1', 'depth1',
              'thts1', 'thti1', 'eps1', 'uhc1',
              'hcond2', 'thickm2', 'elevdn', 'width2', 'depth2',
              'thts2', 'thti2', 'eps2', 'uhc2']

    dtypes = {'rno': np.int, 'node': np.int, 'k': np.int, 'i': np.int, 'j': np.int,
              'iseg': np.int, 'ireach': np.int, 'outreach': np.int, 'line_id': np.int,
              'per': np.int, 'nseg': np.int, 'icalc': np.int, 'outseg': np.int,
              'iupseg': np.int, 'iprior': np.int, 'nstrpts': np.int,
              'name': np.object, 'geometry': np.object}

    # LENUNI = {"u": 0, "f": 1, "m": 2, "c": 3}
    len_const = {1: 1.486, 2: 1.0, 3: 100.}
    # {"u": 0, "s": 1, "m": 2, "h": 3, "d": 4, "y": 5}
    time_const = {1: 1., 2: 60., 3: 3600., 4: 86400., 5: 31557600.}

    # default values
    defaults = {'icalc': 1,
                'roughch': 0.037,
                'strthick': 1,
                'strhc1': 1,
                }

    def __init__(self, reach_data,
                 segment_data=None, grid=None,
                 model_length_units='ft', model_time_units='d',
                 enforce_increasing_nsegs=True,
                 model_name=None,
                 **kwargs):

        self.reach_data = self._setup_reach_data(reach_data)
        self.segment_data = self._setup_segment_data(segment_data)
        self.grid = grid
        self.structured = self.grid.structured

        # convert any modflow6 kwargs to modflow5
        #for k, v in kwargs.items():
        #    if k in sfrdata.mf6names.values():
        #        kwargs[sfrdata.mf5names[k]] = kwargs.pop(k)
        kwargs = {sfrdata.mf5names[k] if k in sfrdata.mf6names else k:
                      v for k, v in kwargs.items()}
        # apply any defaults not in kwargs
        for k, v in sfrdata.defaults.items():
            if k not in kwargs:
                kwargs[k] = v
        # assign kwargs to reach/segment data
        for k, v in kwargs.items():
            if k in self.rdcols:
                self.reach_data[k] = v
            elif k in self.sdcols:
                self.segment_data[k] = v

        # routing
        self._segment_routing = None  # dictionary of routing connections
        self._reach_routing = None # dictionary of rno routing connections
        self._paths = None  # routing sequence from each segment to outlet

        if not self._valid_nsegs(increasing=enforce_increasing_nsegs):
            self.reset_segments()

        # establish rno routing
        # set_outreaches checks for valid rnos and resets them if not
        # resets out reaches either way using segment data
        # not the ideal logic for MODFLOW 6 case where only
        # rno and connections are supplied
        self.set_outreaches()

        # attached instance of flopy ModflowSfr2 package object
        self._ModflowSfr2 = None
        if model_name is None:
            model_name = 'model'
        self.model_name = model_name

        # units
        self._model_length_units = model_length_units
        self.model_time_units = model_time_units


    @property
    def const(self):
        const = self.len_const[self._lenuni] * \
                self.time_const[self._itmuni]
        return const

    @property
    def _itmuni(self):
        """MODFLOW time units code"""
        d = {"u": 0, "s": 1, "m": 2, "h": 3, "d": 4, "y": 5}
        return d[self.model_time_units]

    @property
    def _lenuni(self):
        """MODFLOW length units code"""
        d = {"ft": 1, "m": 2}
        return d[self.model_length_units]

    @property
    def model_length_units(self):
        """Computational lengths units of numerical model."""
        if self.grid is not None:
            self._model_length_units = self.grid.model_units
        return self._model_length_units

    @property
    def crs_units(self):
        """Length units of the coordinate reference system"""
        return self.grid.crs.length_units

    @property
    def segment_routing(self):
        if self._segment_routing is None or self._routing_changed():
            sd = self.segment_data.groupby('per').get_group(0)
            graph = dict(
                zip(sd.nseg, sd.outseg))
            outlets = set(graph.values()).difference(
                set(graph.keys()))  # including lakes
            graph.update({o: 0 for o in outlets})
            self._routing = graph
        return self._routing

    @property
    def rno_routing(self):
        if self._rno_routing is None or self._routing_changed():
            # enforce valid rnos; and create outreach connections
            # from segment data and sequential reach numbering
            # (ireach values also checked and fixed if necesseary)
            self.set_outreaches()
            rd = self.reach_data
            graph = dict(zip(rd.rno, rd.outreach))
            outlets = set(graph.values()).difference(
                set(graph.keys()))  # including lakes
            graph.update({o: 0 for o in outlets})
            self._rno_routing = graph
        return self._rno_routing

    @property
    def ModflowSfr2(self):
        """Flopy ModflowSfr2 instance."""
        if self._ModflowSfr2 is None:
            self._ModflowSfr2 = self.create_ModflowSfr2()
        return self._ModflowSfr2

    @staticmethod
    def get_empty_reach_data(nreaches=0, default_value=0):
        rd = fm.ModflowSfr2.get_empty_reach_data(nreaches,
                                                 default_value=default_value)
        df = pd.DataFrame(rd)
        for c in sfrdata.rdcols:
            if c not in df.columns:
                df[c] = default_value
            elif c != 'geometry':
                df[c] = df[c].astype(sfrdata.dtypes.get(c, np.float32))
        return df[sfrdata.rdcols]

    def _setup_reach_data(self, reach_data):
        rd = sfrdata.get_empty_reach_data(len(reach_data))
        reach_data.index = range(len(reach_data))
        for c in reach_data.columns:
            rd[c] = reach_data[c].astype(sfrdata.dtypes.get(c, np.float32))
            assert rd[c].dtype == sfrdata.dtypes.get(c, np.float32)
        return rd

    @staticmethod
    def get_empty_segment_data(nsegments=0, default_value=0):
        sd = fm.ModflowSfr2.get_empty_segment_data(nsegments,
                                                   default_value=default_value)
        sd = pd.DataFrame(sd)
        sd['per'] = 0
        for c in sd.columns:
            sd[c] = sd[c].astype(sfrdata.dtypes.get(c, np.float32))
        return sd[sfrdata.sdcols]

    def _setup_segment_data(self, segment_data):
        # if no segment_data was provided
        if segment_data is None:
            # create segment_data from iseg and ireach columns in reach data
            # if
            if valid_nsegs(self.reach_data.iseg, increasing=False) and \
                    self.reach_data.outseg.sum() > 0:
                self.reach_data.sort_values(by=['iseg', 'ireach'], inplace=True)
                nss = self.reach_data.iseg.max()
                sd = sfrdata.get_empty_segment_data(nss)
                routing = dict(zip(self.reach_data.iseg, self.reach_data.outseg))
                sd['nseg'] = range(len(sd))
                sd['outseg'] = [routing[s] for s in sd.nseg]
            # create segment_data from reach routing (one reach per segment)
            else:
                has_rno_routing = self._check_reach_routing()
                assert has_rno_routing, \
                    "Reach data must contain rno column with unique, " \
                    "consecutive reach numbers starting at 1. If no " \
                    "segment_data are supplied, segment routing must be" \
                    "included in iseg and outseg columns, or reach data " \
                    "must contain outreach column with routing connections."
                sd = sfrdata.get_empty_segment_data(len(self.reach_data))
                sd['nseg'] = self.reach_data.rno
                sd['outseg'] = self.reach_data.outreach
        # transfer supplied segment data to default template
        else:
            sd = sfrdata.get_empty_segment_data(len(segment_data))
        if 'per' not in segment_data.columns:
            segment_data['per'] = 0
        segment_data.sort_values(by=['per', 'nseg'], inplace=True)
        segment_data.index = range(len(segment_data))
        for c in segment_data.columns:
            sd[c] = segment_data[c].astype(sfrdata.dtypes.get(c, np.float32))
        return sd

    @property
    def paths(self):
        """Dict listing routing sequence for each segment
        in SFR network."""
        if self._paths is None:
            self._set_paths()
            return self._paths
        if self._routing_changed():
            self._segment_routing = None
            self._set_paths()
        return self._paths

    def _set_paths(self):
        routing = self.segment_routing
        self._paths = {seg: find_path(routing, seg) for seg in routing.keys()}

    def _routing_changed(self):
        sd = self.segment_data.groupby('per').get_group(0)
        rd = self.reach_data
        # check if segment routing in dataframe is consistent with routing dict
        segment_routing = dict(zip(sd.nseg, sd.outseg))
        segment_routing_changed = segment_routing != self._segment_routing

        # check if reach routing in dataframe is consistent with routing dict
        reach_routing = dict(zip(rd.rno, rd.outreach))
        reach_routing_changed = reach_routing != self._reach_routing

        # check if segment and reach routing in dataframe are consistent
        consistent = rno_nseg_routing_consistent(sd.nseg, sd.outseg,
                                                 rd.iseg, rd.ireach,
                                                 rd.rno, rd.outeach)
        # return True if the dataframes changed,
        # or are inconsistent between segments and reach numbers
        return segment_routing_changed & reach_routing_changed & ~consistent

    def repair_outsegs(self):
        """Set any outsegs that are not nsegs or lakes to 0 (outlet status)"""
        isasegment = np.in1d(self.segment_data.outseg,
                             self.segment_data.nseg)
        isasegment = isasegment | (self.segment_data.outseg < 0)
        self.segment_data.loc[~isasegment, 'outseg'] = 0

    def reset_segments(self):
        """Reset the segment numbering so that is consecutive,
        starts at 1 and only increases downstream."""
        r = renumber_segments(self.segment_data.nseg,
                              self.segment_data.outseg)
        self.segment_data['nseg'] = [r[s] for s in self.segment_data.nseg]
        self.segment_data['outseg'] = [r[s] for s in self.segment_data.outseg]
        self.reach_data['iseg'] = [r[s] for s in self.reach_data.iseg]
        self.reach_data['outseg'] = [r[s] for s in self.reach_data.outseg]
        self.segment_data.sort_values(by=['per', 'nseg'], inplace=True)
        self.reach_data.sort_values(by=['iseg', 'ireach'], inplace=True)

    def reset_reaches(self):
        """Ensure that the reaches in each segment are numbered
        consecutively starting at 1."""
        self.reach_data.sort_values(by=['iseg', 'ireach'], inplace=True)
        reach_data = self.reach_data
        segment_data = self.segment_data.groupby('per').get_group(0)
        reach_counts = np.bincount(reach_data.iseg)[1:]
        reach_counts = dict(zip(range(1, len(reach_counts) +1),
                                reach_counts))
        ireach = [list(range(1, reach_counts[s] + 1))
                   for s in segment_data.nseg]
        ireach = np.concatenate(ireach)
        self.reach_data['ireach'] = ireach

    def set_outreaches(self):
        """Determine the outreach for each SFR reach (requires a reachID column in reach_data).
        Uses the segment routing specified for the first stress period to route reaches between segments.
        """
        self.reach_data.sort_values(by=['iseg', 'ireach'], inplace=True)
        self.segment_data.sort_values(by=['per', 'nseg'], inplace=True)
        if not self._valid_rnos():
            self.reach_data['rno'] = np.arange(1, len(self.reach_data) + 1)
        self.reset_reaches()  # ensure that each segment starts with reach 1
        self.repair_outsegs()  # ensure that all outsegs are segments, outlets, or negative (lakes)
        rd = self.reach_data
        outseg = self.segment_routing
        reach1IDs = dict(zip(rd[rd.ireach == 1].iseg,
                             rd[rd.ireach == 1].rno))
        ireach = rd.ireach.values
        iseg = rd.iseg.values
        rno = rd.rno.values
        outreach = []
        for i in range(len(rd)):
            # if at the end of reach data or current segment
            if i + 1 == len(rd) or ireach[i+1] == 1:
                nextseg = outseg[iseg[i]]  # get next segment
                if nextseg > 0:  # current reach is not an outlet
                    nextrchid = reach1IDs[nextseg] # get reach 1 of next segment
                else:
                    nextrchid = 0
            else:  # otherwise, it's the next reachID
                nextrchid = rno[i + 1]
            outreach.append(nextrchid)
        self.reach_data['outreach'] = outreach

    def _valid_rnos(self):
        incols = 'rno' in self.reach_data.columns
        arevalid = valid_rnos(self.reach_data.rno.tolist())
        return incols & arevalid

    def _check_reach_routing(self):
        """Cursory check of reach routing."""
        valid_rnos = self._valid_rnos()
        non_zero_outreaches = 'outreach' in self.reach_data.columns & \
        self.reach_data.outreach.sum() > 0
        return valid_rnos & non_zero_outreaches

    def _valid_nsegs(self, increasing=True):
        return valid_nsegs(self.segment_data.nseg,
                           self.segment_data.outseg,
                           increasing=increasing)

    def create_ModflowSfr2(self, m=None, const=None,
                           isfropt=1, # override flopy default of 0
                           unit_number=None,
                           ipakcb=None, istcb2=None,
                           **kwargs
                           ):

        if const is None:
            const = self.const

        if m is None:
            m = fm.Modflow(model_ws='',
                       structured=self.structured)

        # translate segment data
        # populate MODFLOW 2005 segment variables from reach data
        width1 = self.reach_data.groupby('iseg')['width'].min()
        width1.index = width1.index -1 # convert from iseg to zero-based
        width2 = self.reach_data.groupby('iseg')['width'].max()
        width2.index = width2.index - 1  # convert from iseg to zero-based

        self.segment_data['width2'] = width2
        self.segment_data['width1'] = width1

        assert not np.any(np.isnan(self.segment_data))
        # create record array for each stress period
        sd = self.segment_data.groupby('per')
        sd = {per: sd.get_group(per).drop('per', axis=1).to_records(index=False)
              for per in self.segment_data.per.unique()}

        # translate reach data
        flopy_cols = fm.ModflowSfr2.\
            get_default_reach_dtype(structured=self.structured).names
        columns_not_in_flopy = set(self.reach_data.columns).difference(set(flopy_cols))
        rd = self.reach_data.drop(columns_not_in_flopy, axis=1).copy()
        rd = rd.to_records(index=False)
        nstrm = -len(rd)

        return fm.ModflowSfr2(model=m, nstrm=nstrm, const=const,
                              reach_data=rd, segment_data=sd,
                              isfropt=isfropt, unit_number=unit_number,
                              ipakcb=ipakcb, istcb2=istcb2,
                              **kwargs)

        # create a flopy model attribute that can be updated
        # user should be able to assign dynamically
        # (sfrdata.model = m)
        # if there's no model, one is created from grid
        # grid is still used to write output to make it easier to deal with unstructured
        # ModflowSfr2 property allows dynamic access to flopy methods

    def interpolate_to_reaches(self, segvar1, segvar2, per=0):
        """Interpolate values in datasets 6b and 6c to each reach in stream segment

        Parameters
        ----------
        segvar1 : str
            Column/variable name in segment_data array for representing start of segment
            (e.g. hcond1 for hydraulic conductivity)
            For segments with icalc=2 (specified channel geometry); if width1 is given,
            the eigth distance point (XCPT8) from dataset 6d will be used as the stream width.
            For icalc=3, an abitrary width of 5 is assigned.
            For icalc=4, the mean value for width given in item 6e is used.
        segvar2 : str
            Column/variable name in segment_data array for representing start of segment
            (e.g. hcond2 for hydraulic conductivity)
        per : int
            Stress period with segment data to interpolate

        Returns
        -------
        reach_values : 1D array
            One dimmensional array of interpolated values of same length as reach_data array.
            For example, hcond1 and hcond2 could be entered as inputs to get values for the
            strhc1 (hydraulic conductivity) column in reach_data.

        """
        reach_data = self.reach_data
        segment_data = self.segment_data.groupby('per').get_group(per)
        segment_data.sort_values(by='nseg', inplace=True)
        reach_data.sort_values(by=['iseg', 'ireach'], inplace=True)
        rd_groups = reach_data.groupby('iseg')
        sd_groups = segment_data.groupby('nseg')
        reach_values = []
        for seg in segment_data.nseg:
            reaches = rd_groups.get_group(seg)
            dist = (np.cumsum(reaches.rchlen) - 0.5 * reaches.rchlen).values

            fp = [sd_groups.get_group(seg)[segvar1].values[0],
                  sd_groups.get_group(seg)[segvar2].values[0]]
            xp = [dist[0], dist[-1]]
            reach_values += np.interp(dist, xp, fp).tolist()
        return np.array(reach_values)

    def sample_elevations(self, dem, method='buffers',
                          buffer_distance=None,
                          statistic='min'):
        """Computes zonal statistics on a raster for SFR reaches, using
        either buffer polygons around the reach LineStrings, or the model
        grid cell polygons containing each reach.

        Parameters
        ----------
        dem : path to valid raster dataset
            Must be in same Coordinate Reference System as model grid.
        method : str; 'buffers' or 'cell polygons'
            If 'buffers', buffers (with flat caps; cap_style=2 in LineString.buffer())
            will be created around the reach LineStrings (geometry column in reach_data).
        buffer_distance : float
            Buffer distance to apply around reach LineStrings, in crs_units.
        """
        raise NotImplementedError

    def get_slopes(self, default_slope=0.001, minimum_slope=0.0001,
                   maximum_slope=1.):
        """Compute slopes by reach using values in strtop (streambed top) and rchlen (reach length)
        columns of reach_data. The slope for a reach n is computed as strtop(n+1) - strtop(n) / rchlen(n).
        Slopes for outlet reaches are set equal to a default value (default_slope).
        Populates the slope column in reach_data.

        Parameters
        ----------
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
        rd = self.reach_data
        elev = dict(zip(rd.rno, rd.strtop))
        dist = dict(zip(rd.rno, rd.rchlen))
        dnelev = {rid: elev[rd.outreach.values[i]] if rd.outreach.values[i] != 0
        else -9999 for i, rid in enumerate(rd.rno)}
        slopes = np.array(
            [(elev[i] - dnelev[i]) / dist[i] if dnelev[i] != -9999
             else default_slope for i in rd.rno])
        slopes[slopes < minimum_slope] = minimum_slope
        slopes[slopes > maximum_slope] = maximum_slope
        self.reach_data['slope'] = slopes

    @staticmethod
    def from_package(sfrpackagefile, grid, linework):
        """Read SFR package file

        Parameters
        ----------
        sfrpackagefile : file path
            Modflow-2005 or MODFLOW6 SFR package
        grid : sfrmaker.grid instance
        linework : shapefile path or DataFrame
            Contains linestrings for each reach; must have
            segment and reach, or reach number (rno in MODFLOW 6)
            information.

        Returns
        -------
        sfrdata : sfrmaker.sfrdata instance
        """
        pass

    def write_package(self, filename=None, version='mf2005',
                      **kwargs):
        """Write and SFR package file.

        Parameters
        ----------
        version : str
            'mf2005' or 'mf6'
        """
        if version == 'mf2005':
            if len(kwargs) > 0:
                self.create_ModflowSfr2(**kwargs)

            self.ModflowSfr2.write_file(filename=filename)

        elif version == 'mf6':

            from .mf5to6 import mf6sfr
            # instantiate mf6sfr converter object with mf-nwt model/sfr package from flopy
            sfr6 = mf6sfr(m.sfr)

            # write a MODFLOW 6 file
            sfr6.write_file(outpath=mf6_ws)

    def export_sfrlines(self, filename=None):
        """Export shapefiles of linework"""
        if filename is None:
            filename = self.model_name + '_sfrlines.shp'
        df2shp(self.reach_data, filename, epsg=self.grid.crs.epsg)

    def export_routing(self, filename=None):
        """Export linework shapefile showing all routing connections between SFR reaches.
        A length field containing the distance between connected reaches
        can be used to filter for the longest connections in a GIS.
        """
        if filename is None:
            filename = self.model_name + '_sfr_routing.shp'
        rd = self.reach_data[['node', 'iseg', 'ireach', 'rno', 'outreach']].copy()
        rd.sort_values(by='rno', inplace=True)
        cellgeoms = self.grid.df.loc[rd.node.values, 'geometry']

        # get the cell centers for each reach
        x0 = [g.centroid.x for g in cellgeoms]
        y0 = [g.centroid.y for g in cellgeoms]
        loc = dict(zip(rd.rno, zip(x0, y0)))

        # make lines of the reach connections between cell centers
        geoms = []
        lengths = []
        for i, r in enumerate(rd.rno.values):
            x0, y0 = loc[r]
            outreach = rd.outreach.values[i]
            if outreach == 0:
                x1, y1 = x0, y0
            else:
                x1, y1 = loc[outreach]
            geoms.append(LineString([(x0, y0), (x1, y1)]))
            lengths.append(np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2))
        lengths = np.array(lengths)
        rd['length'] = lengths
        rd['geometry'] = geoms
        df2shp(rd, filename, epsg=self.grid.crs.epsg)