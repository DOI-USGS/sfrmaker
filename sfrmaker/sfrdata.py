import sys
sys.path.append('/Users/aleaf/Documents/GitHub/flopy3')
import numpy as np
import pandas as pd
import flopy
from .utils import renumber_segments, find_path, make_graph
from .checks import valid_rnos, rno_nseg_routing_consistent

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
    model_time_units : str
        's': seconds
        'm': minutes
        'h': hours
        'd': days
        'y': years
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

    def __init__(self, reach_data=None,
                 segment_data=None, grid=None,
                 model_length_units='ft', model_time_units='d',
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
        check that segment numbering is valid
        # call renumber_segments if it isn't, in attached method simliar to self.set_outreaches()
        self.set_outreaches() # establish rno routing
        self._segment_routing = None  # dictionary of routing connections
        self._reach_routing = None # dictionary of rno routing connections
        self._paths = None  # routing sequence from each segment to outlet

        # attached instance of flopy ModflowSfr2 package object
        self._ModflowSfr2 = None

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
    def get_empty_segment_data(nsegments=0):
        sd = fm.ModflowSfr2.get_empty_segment_data(nsegments)
        sd = pd.DataFrame(sd)
        sd['per'] = 0
        for c in sd.columns:
            sd[c] = sd[c].astype(sfrdata.dtypes.get(c, np.float32))
        return sd[sfrdata.sdcols]

    def _setup_segment_data(self, segment_data):
        sd = sfrdata.get_empty_segment_data(len(segment_data))
        segment_data.sort_values(by='nseg', inplace=True)
        segment_data.index = range(len(segment_data))
        for c in segment_data.columns:
            sd[c] = segment_data[c].astype(sfrdata.dtypes.get(c, np.float32))
        return sd

    @property
    def routing(self):
        could base routing off of reach connections in reach_data
        set segment data from that
        but this might just be confusing
        either day, need to also update other dataframe when routing changes
        if self._routing is None or self._routing_changed():
            toid = self.df.toid.values
            # check whether or not routing is
            # many-to-one or one-to-one (no diversions)
            # squeeze it down
            to_one = False
            # if below == True, all toids are scalar or length 1 lists
            to_one = np.isscalar(np.squeeze(toid)[0])
            # if not, try converting any scalars to lists
            if not to_one:
                toid = [[l] if np.isscalar(l) else l for l in toid]
                to_one = np.isscalar(np.squeeze(toid)[0])
            toid = np.squeeze(toid)
            self._routing = make_graph(self.df.id.values, toid,
                                       one_to_many=not to_one)
        return self._routing

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
        self.segment_data.loc[~isasegment, 'outseg'] = 0.

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
        if not self._valid_rnos():
            self.reach_data['rno'] = np.arange(1, len(self.reach_data) + 1)
        self.reset_reaches()  # ensure that each segment starts with reach 1
        self.repair_outsegs()  # ensure that all outsegs are segments, outlets, or negative (lakes)
        rd = self.reach_data
        outseg = self.segment_routing
        reach1IDs = dict(zip(rd[rd.ireach == 1].iseg,
                             rd[rd.ireach == 1].rno))
        outreach = []
        for i in range(len(rd)):
            # if at the end of reach data or current segment
            if i + 1 == len(rd) or rd.ireach[i + 1] == 1:
                nextseg = outseg[rd.iseg[i]]  # get next segment
                if nextseg > 0:  # current reach is not an outlet
                    nextrchid = reach1IDs[
                        nextseg]  # get reach 1 of next segment
                else:
                    nextrchid = 0
            else:  # otherwise, it's the next reachID
                nextrchid = rd.rno[i + 1]
            outreach.append(nextrchid)
        self.reach_data['outreach'] = outreach

    def _valid_rnos(self):
        incols = 'rno' in self.reach_data.columns
        arevalid = valid_rnos(self.reach_data.rno.tolist())
        return incols & arevalid

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


