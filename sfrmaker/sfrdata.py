import sys

sys.path.append('/Users/aleaf/Documents/GitHub/flopy3')
import os
import time
import numpy as np
import pandas as pd
import rasterio
from rasterstats import zonal_stats
from shapely.geometry import LineString
import flopy
from gisutils import df2shp
from sfrmaker.routing import find_path, renumber_segments
from .checks import valid_rnos, valid_nsegs, rno_nseg_routing_consistent
from .elevations import smooth_elevations
from .flows import add_to_perioddata
from .gis import export_reach_data, project, crs
from .grid import StructuredGrid
from .observations import write_mf6_sfr_obsfile, add_observations
from .units import convert_length_units, get_length_units, itmuni_values, lenuni_values
import sfrmaker
from sfrmaker.mf5to6 import segment_data_to_period_data

fm = flopy.modflow
mf6 = flopy.mf6


class SFRData:
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
        'meters' or 'feet'
    model_time_units : str
        's': seconds
        'meters': minutes
        'h': hours
        'd': days
        'y': years
    enforce_increasing_nsegs : bool
        If True, segment numbering is checked to ensure
        that it only increases downstream, and reset if it doesnt.
    package_name : str
            Base name for writing sfr output.
    kwargs : keyword arguments
        Optional values to assign globally to SFR variables. For example
        icalc=1 would assign all segments an icalc value of 1. For default
        values see the sfrdata.defaults dictionary. Default values can be
        assigned using MODFLOW-2005 or MODFLOW-6 terminology.
    """

    # conversions to MODFLOW6 variable names
    mf6names = {'rno': 'rno',
                # 'node': 'cellid',
                'rchlen': 'rlen',
                'width': 'rwid',
                'slope': 'rgrd',
                'strtop': 'rtp',
                'strthick': 'rbth',
                'strhc1': 'rhk',
                'roughch': 'man',
                'flow': 'inflow',
                'pptsw': 'rainfall',
                'etsw': 'evaporation',
                'runoff': 'runoff',
                'depth1': 'depth1',
                'depth2': 'depth2'}

    mf5names = {v: k for k, v in mf6names.items()}

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
    len_const = {0: 1.0, 1: 1.486, 2: 1.0, 3: 100.}
    # {"u": 0, "s": 1, "m": 2, "h": 3, "d": 4, "y": 5}
    time_const = {1: 1., 2: 60., 3: 3600., 4: 86400., 5: 31557600.}

    # default values
    defaults = {'icalc': 1,
                'roughch': 0.037,
                'strthick': 1,
                'strhc1': 1,
                }

    def __init__(self, reach_data=None,
                 segment_data=None, grid=None, sr=None,
                 model=None,
                 isfr=None,
                 model_length_units="undefined", model_time_units='d',
                 enforce_increasing_nsegs=True,
                 package_name=None,
                 **kwargs):

        # attributes
        self._crs = None
        self._period_data = None
        self._observations = None
        self._observations_filename = None

        # convert any modflow6 kwargs to modflow5
        kwargs = {SFRData.mf5names[k] if k in SFRData.mf6names else k:
                      v for k, v in kwargs.items()}
        # update default values (used in segment and reach data setup)
        self.defaults.update(kwargs)

        self.reach_data = self._setup_reach_data(reach_data)
        self.segment_data = self._setup_segment_data(segment_data)
        self.isfropt0_to_1()  # distribute any isfropt=0 segment data to reaches

        if grid is None and sr is not None:
            print('\nCreating grid class instance from flopy SpatialReference...')
            ta = time.time()
            grid = StructuredGrid.from_sr(sr, isfr=isfr)
            print("grid class created in {:.2f}s\n".format(time.time() - ta))
        elif isinstance(grid, flopy.discretization.grid.Grid):
            print('\nCreating grid class instance from flopy modelgrid...')
            ta = time.time()
            grid = StructuredGrid.from_modelgrid(grid, isfr=isfr)
            print("grid class created in {:.2f}s\n".format(time.time() - ta))

        # print grid information to screen
        print(grid)
        self.grid = grid

        # routing
        self._segment_routing = None  # dictionary of routing connections
        self._reach_routing = None  # dictionary of rno routing connections
        self._paths = None  # routing sequence from each segment to outlet

        if not self._valid_nsegs(increasing=enforce_increasing_nsegs):
            self.reset_segments()

        # establish rno routing
        # set_outreaches checks for valid rnos and resets them if not
        # resets out reaches either way using segment data
        # not the ideal logic for MODFLOW 6 case where only
        # rno and connections are supplied
        self.set_outreaches()
        self.get_slopes()

        # units
        self.model_length_units = get_length_units(model_length_units, grid, model)
        self.model_time_units = model_time_units

        self.model = model  # attached flopy model instance
        # self._ModflowSfr2 = None # attached instance of flopy modflow_sfr2 package object
        self.package_name = package_name

    @property
    def const(self):
        const = self.len_const[self._lenuni] * \
                self.time_const[self._itmuni]
        return const

    @property
    def _itmuni(self):
        """MODFLOW time units code"""
        return itmuni_values[self.model_time_units]

    @property
    def _lenuni(self):
        """MODFLOW length units code"""
        return lenuni_values[self.model_length_units]

    @property
    def structured(self):
        if self.grid is not None:
            return self.grid._structured
        else:
            return True

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model):
        self._model = model

        # update the sfr package object as well
        self.create_modflow_sfr2(model)

    @property
    def package_name(self):
        return self._package_name

    @package_name.setter
    def package_name(self, package_name):
        if package_name is None and self.model is not None:
            self._package_name = self.model.name
        elif package_name is None:
            self._package_name = 'model'
        else:
            self._package_name = package_name

    @property
    def crs(self):
        if self._crs is None:
            self._crs = self.grid.crs
        return self._crs

    @property
    def crs_units(self):
        """Length units of the coordinate reference system"""
        return self.crs.length_units

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
    def modflow_sfr2(self):
        """Flopy modflow_sfr2 instance."""
        if self._ModflowSfr2 is None:
            self.create_modflow_sfr2()
        return self._ModflowSfr2

    @classmethod
    def get_empty_reach_data(cls, nreaches=0, default_value=0):
        rd = fm.ModflowSfr2.get_empty_reach_data(nreaches,
                                                 default_value=default_value)
        df = pd.DataFrame(rd)
        for c in SFRData.rdcols:
            if c not in df.columns:
                df[c] = default_value
            elif c != 'geometry':
                df[c] = df[c].astype(cls.dtypes.get(c, np.float32))
        return df[cls.rdcols]

    def _setup_reach_data(self, reach_data):
        rd = SFRData.get_empty_reach_data(len(reach_data))
        reach_data.index = range(len(reach_data))
        for c in reach_data.columns:
            rd[c] = reach_data[c].astype(SFRData.dtypes.get(c, np.float32))
            assert rd[c].dtype == SFRData.dtypes.get(c, np.float32)
        # assign kwargs to reach data
        for k, v in self.defaults.items():
            if k in self.rdcols and k not in reach_data.columns:
                rd[k] = v
        return rd

    @classmethod
    def get_empty_segment_data(cls, nsegments=0, default_value=0):
        sd = fm.ModflowSfr2.get_empty_segment_data(nsegments,
                                                   default_value=default_value)
        sd = pd.DataFrame(sd)
        sd['per'] = 0
        for c in sd.columns:
            sd[c] = sd[c].astype(cls.dtypes.get(c, np.float32))
        return sd[cls.sdcols]

    def _setup_segment_data(self, segment_data):
        # if no segment_data was provided
        if segment_data is None:
            # create segment_data from iseg and ireach columns in reach data
            # if
            if valid_nsegs(self.reach_data.iseg, increasing=False) and \
                    self.reach_data.outseg.sum() > 0:
                self.reach_data.sort_values(by=['iseg', 'ireach'], inplace=True)
                nss = self.reach_data.iseg.max()
                sd = SFRData.get_empty_segment_data(nss)
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
                sd = SFRData.get_empty_segment_data(len(self.reach_data))
                sd['nseg'] = self.reach_data.rno
                sd['outseg'] = self.reach_data.outreach
        # transfer supplied segment data to default template
        else:
            sd = SFRData.get_empty_segment_data(len(segment_data))
        if 'per' not in segment_data.columns:
            segment_data['per'] = 0
        segment_data.sort_values(by=['per', 'nseg'], inplace=True)
        segment_data.index = range(len(segment_data))
        for c in segment_data.columns:
            sd[c] = segment_data[c].astype(SFRData.dtypes.get(c, np.float32))

        # assign defaults to segment data
        for k, v in self.defaults.items():
            if k in self.sdcols and k not in segment_data.columns:
                sd[k] = v

        # add outsegs to reach_data
        routing = dict(zip(sd.nseg, sd.outseg))
        self.reach_data['outseg'] = [routing[s] for s in self.reach_data.iseg]
        return sd

    @property
    def observations(self):
        if self._observations is None:
            self._observations = pd.DataFrame(columns=['obsname', 'obstype', 'rno'])
        return self._observations

    @property
    def observations_file(self):
        if self._observations_filename is None:
            self._observations_filename = self.package_name + '.sfr.obs'
        return self._observations_filename

    @observations_file.setter
    def observations_file(self, observations_filename):
        self._observations_filename = observations_filename

    @property
    def period_data(self):
        if self._period_data is None:
            self._period_data = self._get_period_data()
        if not np.array_equal(self._period_data.index.values,
                              self._period_data.rno.values):
            self._period_data.index = self._period_data.rno
            self._period_data.index.name = None
        return self._period_data

    def _get_period_data(self):
        print('converting segment data to period data...')
        return segment_data_to_period_data(self.segment_data, self.reach_data)

    def add_to_perioddata(self, data, flowline_routing=None,
                          variable='inflow',
                          line_id_column_in_data=None,
                          rno_column_in_data=None,
                          period_column_in_data='per',
                          variable_column_in_data='Q_avg'):
        return add_to_perioddata(self, data, flowline_routing=flowline_routing,
                                 variable=variable,
                                 line_id_column_in_data=line_id_column_in_data,
                                 rno_column_in_data=rno_column_in_data,
                                 period_column_in_data=period_column_in_data,
                                 variable_column_in_data=variable_column_in_data)

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
                                                 rd.rno, rd.outreach)
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
        self.segment_data.index = np.arange(len(self.segment_data))
        assert np.array_equal(self.segment_data.loc[self.segment_data.per == 0, 'nseg'].values,
                              self.segment_data.loc[self.segment_data.per == 0].index.values + 1)
        self.reach_data.sort_values(by=['iseg', 'ireach'], inplace=True)

    def reset_reaches(self):
        """Ensure that the reaches in each segment are numbered
        consecutively starting at 1."""
        self.reach_data.sort_values(by=['iseg', 'ireach'], inplace=True)
        reach_data = self.reach_data
        segment_data = self.segment_data.groupby('per').get_group(0)
        reach_counts = np.bincount(reach_data.iseg)[1:]
        reach_counts = dict(zip(range(1, len(reach_counts) + 1),
                                reach_counts))
        ireach = [list(range(1, reach_counts[s] + 1))
                  for s in segment_data.nseg]
        ireach = np.concatenate(ireach)
        self.reach_data['ireach'] = ireach

    def set_outreaches(self):
        """Determine the outreach for each SFR reach (requires a rno column in reach_data).
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
            if i + 1 == len(rd) or ireach[i + 1] == 1:
                nextseg = outseg[iseg[i]]  # get next segment
                if nextseg > 0:  # current reach is not an outlet
                    nextrchid = reach1IDs[nextseg]  # get reach 1 of next segment
                else:
                    nextrchid = 0
            else:  # otherwise, it's the next rno
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
        sd0 = self.segment_data.loc[self.segment_data.per == 0]
        return valid_nsegs(sd0.nseg,
                           sd0.outseg,
                           increasing=increasing)

    def create_modflow_sfr2(self, model=None, const=None,
                            isfropt=1,  # override flopy default of 0
                            unit_number=None,
                            ipakcb=None, istcb2=None,
                            **kwargs
                            ):

        if const is None:
            const = self.const

        if model is not None and model.version != 'mf6':
            m = model
        # create an flopy mf2005 model instance attached to modflow_sfr2 object
        # this is a parallel model instance to self.model, that is only
        # accessible through modflow_sfr2.parent. As long as this method is
        # called by the @model.setter; this instance should have the same
        # dis and ibound (idomain) as self.model.
        # The modflow_sfr2 instance is used as basis for writing packages,
        # because it includes many features, like checking and exporting,
        # that ModflowGwfsfr doesn't have)
        # ibound in BAS package is used by mf5to6 converter
        # to fill in required "None" values when writing mf6 input
        elif model is None or model.version == 'mf6':
            model_ws = '.' if model is None else model.model_ws
            m = fm.Modflow(model_ws=model_ws, structured=self.structured)
            if model is not None:
                nper = 1
                if 'tdis' in model.simulation.package_key_dict.keys():
                    tdis = model.simulation.package_key_dict['tdis']
                    nper = tdis.nper.array
                if 'dis' in model.package_dict.keys():
                    dis = model.dis
                    botm = dis.botm.array.copy()
                    idomain = dis.idomain.array.copy()
                    nlay, nrow, ncol = botm.shape
                    dis = fm.ModflowDis(m, nlay=nlay, nrow=nrow, ncol=ncol, nper=nper,
                                        top=dis.top.array.copy(),
                                        botm=botm)
                    bas = fm.ModflowBas(m, ibound=idomain)

        # translate segment data
        # populate MODFLOW 2005 segment variables from reach data if they weren't entered
        if self.segment_data[['width1', 'width2']].sum().sum() == 0:
            raise NotImplementedError('Double check indexing below before using this option.')
            width1 = self.reach_data.groupby('iseg')['width'].min().to_dict()
            width2 = self.reach_data.groupby('iseg')['width'].max().to_dict()

            self.segment_data['width1'] = [width1[s] for s in self.segment_data.nseg]
            self.segment_data['width2'] = [width2[s] for s in self.segment_data.nseg]

        assert not np.any(np.isnan(self.segment_data))
        # create record array for each stress period
        sd = self.segment_data.groupby('per')
        sd = {per: sd.get_group(per).drop('per', axis=1).to_records(index=False)
              for per in self.segment_data.per.unique()}

        # translate reach data
        flopy_cols = fm.ModflowSfr2. \
            get_default_reach_dtype(structured=self.structured).names
        columns_not_in_flopy = set(self.reach_data.columns).difference(set(flopy_cols))
        rd = self.reach_data.drop(columns_not_in_flopy, axis=1).copy()
        rd = rd.to_records(index=False)
        nstrm = -len(rd)

        self._ModflowSfr2 = fm.ModflowSfr2(model=m, nstrm=nstrm, const=const,
                                           reach_data=rd, segment_data=sd,
                                           isfropt=isfropt, unit_number=unit_number,
                                           ipakcb=ipakcb, istcb2=istcb2,
                                           **kwargs)

        # if model.version == 'mf6':
        #    self._ModflowSfr2.parent = model
        return self._ModflowSfr2

    def create_mf6sfr(self, model=None, unit_conversion=None,
                      stage_filerecord=None,
                      budget_filerecord=None,
                      flopy_rno_input_is_zero_based=True,
                      **kwargs
                      ):

        if unit_conversion is None:
            unit_conversion = self.const
        if stage_filerecord is None:
            stage_filerecord = '{}.sfr.stage.bin'.format(self.package_name)
        if budget_filerecord is None:
            budget_filerecord = '{}.sfr.cbc'.format(self.package_name)

        if model is not None and model.version == 'mf6':
            m = model
            if 'tdis' in model.simulation.package_key_dict.keys():
                tdis = model.simulation.package_key_dict['tdis']
                nper = tdis.nper.array
        # create an flopy mf2005 model instance attached to modflow_sfr2 object
        # this is a parallel model instance to self.model, that is only
        # accessible through modflow_sfr2.parent. As long as this method is
        # called by the @model.setter; this instance should have the same
        # dis and ibound (idomain) as self.model.
        # The modflow_sfr2 instance is used as basis for writing packages,
        # because it includes many features, like checking and exporting,
        # that ModflowGwfsfr doesn't have)
        # ibound in BAS package is used by mf5to6 converter
        # to fill in required "None" values when writing mf6 input
        elif model is None or model.version != 'mf6':
            # create simulation
            sim = flopy.mf6.MFSimulation(version='mf6', exe_name='mf6',
                                         sim_ws='')
            m = flopy.mf6.ModflowGwf(sim)

        # create an sfrmaker.mf6sfr instance
        from .mf5to6 import mf6sfr

        sfr6 = mf6sfr(self.modflow_sfr2)

        # package data
        # An error occurred when storing data "packagedata" in a recarray.
        # packagedata data is a one or two dimensional list containing the variables
        # "<rno> <cellid> <rlen> <rwid> <rgrd> <rtp> <rbth> <rhk> <man> <ncon> <ustrf> <ndv>"
        # (some variables may be optional, see MF6 documentation)
        packagedata = sfr6.packagedata.copy()
        if self.structured:
            columns = packagedata.drop(['k', 'i', 'j', 'idomain'], axis=1).columns.tolist()
            packagedata['cellid'] = list(zip(packagedata.k,
                                             packagedata.i,
                                             packagedata.j))
            columns.insert(1, 'cellid')

        connectiondata = [(rno, *sfr6.connections[rno])
                          for rno in sfr6.packagedata.rno if rno in sfr6.connections]
        # as of 9/12/2019, flopy.mf6.modflow.ModflowGwfsfr requires zero-based input for rno
        if flopy_rno_input_is_zero_based and self.modflow_sfr2.reach_data['reachID'].min() == 1:
            packagedata['rno'] -= 1
            zero_based_connectiondata = []
            for item in connectiondata:
                zb_items = tuple(i - 1 if i > 0 else i + 1 for i in item)
                zero_based_connectiondata.append(zb_items)
            connectiondata = zero_based_connectiondata
        assert packagedata['rno'].min() == 0
        assert np.min(list(map(np.min, map(np.abs, connectiondata)))) < 1

        # set cellids to None for unconnected reaches or where idomain == 0
        # note: as of 9/12/2019, this results in an error during ModflowGwfsfr package construction
        # with either Nonetype or "NONE"
        # even though MODFLOW 6 requires it
        # unconnected = ~packagedata['rno'].isin(np.array(list(sfr6.connections.keys())) - 1).values
        # inactive = m.dis.idomain.array[packagedata.k.values,
        #                               packagedata.i.values,
        #                               packagedata.j.values]
        # packagedata.loc[unconnected | inactive, 'cellid'] = None
        packagedata = packagedata[columns].values.tolist()

        period_data = None
        if sfr6.period_data is not None:
            # TODO: add method to convert period_data df to MF6 input
            # raise NotImplemented("Support for mf6 period_data input not implemented yet. "
            #                     "Use sfrdata.write_package(version='mf6') instead.")
            pass

        mf6sfr = mf6.ModflowGwfsfr(model=m, unit_conversion=unit_conversion,
                                   stage_filerecord=stage_filerecord,
                                   budget_filerecord=budget_filerecord,
                                   nreaches=len(self.reach_data),
                                   packagedata=packagedata,
                                   connectiondata=connectiondata,
                                   diversions=None,  # TODO: add support for diversions
                                   perioddata=period_data,  # TODO: add support for creating mf6 perioddata input
                                   )
        return mf6sfr

    def add_observations(self, data, flowline_routing=None,
                         obstype=None,
                         line_id_column_in_data=None,
                         rno_column_in_data=None,
                         obstype_column_in_data=None,
                         obsname_column_in_data='site_no'):
        return add_observations(self, data, flowline_routing=flowline_routing,
                                obstype=obstype,
                                line_id_column_in_data=line_id_column_in_data,
                                rno_column_in_data=rno_column_in_data,
                                obstype_column_in_data=obstype_column_in_data,
                                obsname_column_in_data=obsname_column_in_data)

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
        from .utils import interpolate_to_reaches

        reach_data = self.reach_data
        segment_data = self.segment_data.groupby('per').get_group(per)
        segment_data.sort_values(by='nseg', inplace=True)
        reach_data.sort_values(by=['iseg', 'ireach'], inplace=True)

        return interpolate_to_reaches(reach_data, segment_data,
                                      segvar1, segvar2,
                                      reach_data_group_col='iseg',
                                      segment_data_group_col='nseg'
                                      )

    def isfropt0_to_1(self):
        """transfer isfropt=0 segment properties to reaches,
        using linear interpolation.
        """
        snames = {'strtop': ('elevup', 'elevdn'),
                  'strthick': ('thickm1', 'thickm2'),
                  'strhc1': ('hcond1', 'hcond2'),
                  'thts': ('thts1', 'thts2'),
                  'thti': ('thti1', 'thti2'),
                  'eps': ('eps1', 'eps2'),
                  'uhc': ('uhc1', 'uhc2'),
                  }
        sd = self.segment_data.loc[self.segment_data.per == 0]
        for col, sdcols in snames.items():
            if self.reach_data[col].sum() == 0 and \
                    sd[[*sdcols]].values.sum(axis=(0, 1)) != 0.:
                self.reach_data[col] = self.interpolate_to_reaches(*sdcols)

    def sample_reach_elevations(self, dem, dem_z_units=None,
                                method='buffers',
                                buffer_distance=100,
                                statistic='min',
                                smooth=True
                                ):
        """Computes zonal statistics on a raster for SFR reaches, using
        either buffer polygons around the reach LineStrings, or the model
        grid cell polygons containing each reach.

        Parameters
        ----------
        dem : path to valid raster dataset
            Must be in same Coordinate Reference System as model grid.
        dem_z_units : str
            Elevation units for DEM ('feet' or 'meters'). If None, units
            are assumed to be same as model (default).
        method : str; 'buffers' or 'cell polygons'
            If 'buffers', buffers (with flat caps; cap_style=2 in LineString.buffer())
            will be created around the reach LineStrings (geometry column in reach_data).
        buffer_distance : float
            Buffer distance to apply around reach LineStrings, in crs_units.
        statistic : str
            "stats" argument to rasterstats.zonal_stats.
            "min" is recommended for streambed elevations (default).
        smooth : bool
            Run sfrmaker.elevations.smooth_elevations on sampled elevations
            to ensure that they decrease monotonically in the downstream direction
            (default=True).

        Returns
        -------
        elevs : dict of sampled elevations keyed by reach number
        """

        # get the CRS and pixel size for the DEM
        with rasterio.open(dem) as src:
            proj_str = src.crs.to_string()
            epsg = src.crs.to_epsg()
            raster_crs = crs(epsg=epsg, proj_str=proj_str)

            # make sure buffer is large enough for DEM pixel size
            buffer_distance = np.max([np.sqrt(src.res[0] *
                                              src.res[1]) * 1.01,
                                      buffer_distance])

        if method == 'buffers':
            assert isinstance(self.reach_data.geometry[0], LineString), \
                "Need LineString geometries in reach_data.geometry column to use buffer option."
            features = [g.buffer(buffer_distance) for g in self.reach_data.geometry]
            txt = 'buffered LineStrings'
        elif method == 'cell polygons':
            assert self.grid is not None, \
                "Need an attached sfrmaker.Grid instance to use cell polygons option."
            features = self.grid.df.loc[self.reach_data.node, 'geometry'].tolist()
            txt = method

        # reproject features if they're not in the same crs
        if raster_crs != self.crs:
            features = project(features,
                               self.crs.proj_str,
                               raster_crs.proj_str)

        print('running rasterstats.zonal_stats on {}...'.format(txt))
        t0 = time.time()
        results = zonal_stats(features,
                              dem,
                              stats='min')
        elevs = [r['min'] for r in results]
        print("finished in {:.2f}s\n".format(time.time() - t0))

        if all(v is None for v in elevs):
            raise Exception('No {} intersected with {}. Check projections.'.format(txt, dem))
        if any(v is None for v in elevs):
            raise Exception('Some {} not intersected with {}. '
                            'Check that DEM covers the area of the stream network.'
                            '.'.format(txt, dem))

        if smooth:
            elevs = smooth_elevations(self.reach_data.rno.tolist(),
                                      self.reach_data.outreach.tolist(),
                                      elevs)
        else:
            elevs = dict(zip(self.reach_data.rno, elevs))
        return elevs

    def set_streambed_top_elevations_from_dem(self, dem, dem_z_units=None,
                                              method='buffers',
                                              **kwargs):
        """Set streambed top elevations from a DEM raster.
        Runs sfrdata.sample_reach_elevations

        Parameters
        ----------
        dem : path to valid raster dataset
            Must be in same Coordinate Reference System as model grid.
        dem_z_units : str
            Elevation units for DEM ('feet' or 'meters'). If None, units
            are assumed to be same as model (default).
        method : str; 'buffers' or 'cell polygons'
            If 'buffers', buffers (with flat caps; cap_style=2 in LineString.buffer())
            will be created around the reach LineStrings (geometry column in reach_data).
        kwargs : keyword arguments to sfrdata.sample_reach_elevations

        Returns
        -------
        updates strtop column of sfrdata.reach_data
        """
        sampled_elevs = self.sample_reach_elevations(dem=dem, method=method,
                                                     **kwargs)
        if dem_z_units is None:
            dem_z_units = self.model_length_units
        mult = convert_length_units(dem_z_units, self.model_length_units)
        self.reach_data['strtop'] = [sampled_elevs[rno]
                                     for rno in self.reach_data['rno'].values]
        self.reach_data['strtop'] *= mult

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
        assert rd.outreach.sum() > 0, "requires reach routing, must be called after set_outreaches()"
        elev = dict(zip(rd.rno, rd.strtop))
        dist = dict(zip(rd.rno, rd.rchlen))
        dnelev = {rid: elev[rd.outreach.values[i]] if rd.outreach.values[i] != 0
        else -9999 for i, rid in enumerate(rd.rno)}
        slopes = np.array(
            [(elev[i] - dnelev[i]) / dist[i] if dnelev[i] != -9999 and dist[i] > 0
             else default_slope for i in rd.rno])
        slopes[slopes < minimum_slope] = minimum_slope
        slopes[slopes > maximum_slope] = maximum_slope
        self.reach_data['slope'] = slopes

    @classmethod
    def from_package(cls, sfrpackagefile, grid, linework):
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

    @classmethod
    def from_tables(cls, reach_data, segment_data,
                    grid=None, sr=None, isfr=None):
        reach_data = pd.read_csv(reach_data)
        segment_data = pd.read_csv(segment_data)
        return cls(reach_data=reach_data, segment_data=segment_data,
                   grid=grid, sr=sr,
                   isfr=isfr)

    def write_package(self, filename=None, version='mf2005', idomain=None,
                      options=None, write_observations_input=True,
                      **kwargs):
        """Write and SFR package file.

        Parameters
        ----------
        version : str
            'mf2005' or 'mf6'
        """
        # recreate the flopy package object in case it changed
        self.create_modflow_sfr2(model=self.model, **kwargs)
        header_txt = "#SFR package created by SFRmaker v. {}, " \
                     "via FloPy v. {}\n".format(sfrmaker.__version__, flopy.__version__)
        header_txt += "#model length units: {}, model time units: {}".format(self.model_length_units,
                                                                             self.model_time_units)
        self.modflow_sfr2.heading = header_txt

        if filename is None:
            filename = self.modflow_sfr2.fn_path
        if version == 'mf2005':
            self.modflow_sfr2.write_file(filename=filename)

        elif version == 'mf6':

            # instantiate mf6sfr converter object with mf-nwt model/sfr package from flopy
            from .mf5to6 import mf6sfr
            if options is None:
                # save budget and stage output by default
                options = ['save_flows',
                           'BUDGET FILEOUT {}.cbc'.format(filename),
                           'STAGE FILEOUT {}.stage.bin'.format(filename),
                           ]

            if write_observations_input and len(self.observations) > 0:
                obs_input_filename = filename + '.obs'
                self.write_mf6_sfr_obsfile(filename=obs_input_filename)
                options.append('OBS6 FILEIN {}'.format(obs_input_filename))

            sfr6 = mf6sfr(self.modflow_sfr2, period_data=self.period_data,
                          idomain=idomain,
                          options=options)

            # write a MODFLOW 6 file
            sfr6.write_file(filename=filename)

    def write_tables(self, filepath='./sfr'):

        filepath, file_extension = os.path.splitext(filepath)
        self.reach_data.drop('geometry', axis=1).to_csv('{}_reach_data.csv'.format(filepath), index=False)
        self.segment_data.to_csv('{}_segment_data.csv'.format(filepath), index=False)

    def write_mf6_sfr_obsfile(self,
                              filename=None,
                              sfr_output_filename=None):
        if filename is None:
            filename = self.observations_file
        else:
            self.observations_file = filename
        if sfr_output_filename is None:
            sfr_output_filename = filename + '.output.csv'
        return write_mf6_sfr_obsfile(self.observations,
                                     filename,
                                     sfr_output_filename)

    def export_cells(self, filename=None, nodes=None):
        """Export shapefile of model cells with stream reaches."""
        if filename is None:
            filename = self.package_name + '_sfrcells.shp'
        export_reach_data(self.reach_data, self.grid, filename,
                          nodes=nodes, geomtype='polygon')

    def export_outlets(self, filename=None):
        """Export shapefile of model cells with stream reaches."""
        if filename is None:
            filename = self.package_name + '_sfr_outlets.shp'
        nodes = self.reach_data.loc[self.reach_data.outreach == 0, 'node'].values
        export_reach_data(self.reach_data, self.grid, filename,
                          nodes=nodes, geomtype='point')

    def export_lines(self, filename=None):
        """Export shapefile of linework"""
        if filename is None:
            filename = self.package_name + '_sfrlines.shp'
        rd = self.reach_data
        assert 'geometry' in rd.columns and \
               isinstance(rd.geometry.values[0], LineString), \
            "No LineStrings in reach_data.geometry"
        df2shp(rd, filename, epsg=self.grid.crs.epsg,
               proj_str=self.grid.crs.proj_str)

    def export_routing(self, filename=None):
        """Export linework shapefile showing all routing connections between SFR reaches.
        A length field containing the distance between connected reaches
        can be used to filter for the longest connections in a GIS.
        """
        if filename is None:
            filename = self.package_name + '_sfr_routing.shp'
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
        df2shp(rd, filename, epsg=self.grid.crs.epsg,
               proj_str=self.grid.crs.proj_str)

    def export_transient_variable(self, varname, filename=None):
        """Export point shapefile showing locations with
        a given segment_data variable applied. For example, segments
        where streamflow is entering or leaving the upstream end of a stream segment (FLOW)
        or where RUNOFF is applied. Cell centroids of the first reach of segments with
        non-zero terms of varname are exported; values of varname are exported by
        stress period in the attribute fields (e.g. flow0, flow1, flow2... for FLOW
        in stress periods 0, 1, 2...

        Parameters
        ----------
        f : str, filename
        varname : str
            Variable in SFR Package dataset 6a (see SFR package documentation)

        """
        if filename is None:
            filename = self.package_name + '_sfr_{}.shp'.format(varname)

        # if the data are in mf2005 format (by segment)
        sd = self.segment_data.sort_values(by=['per', 'nseg'])

        # pivot the segment data to segments x periods with values of varname
        just_the_values = sd.pivot(index='nseg', columns='per', values=varname)
        hasvalues = np.nansum(just_the_values, axis=1) > 0
        df = just_the_values.loc[hasvalues]
        if len(df) == 0:
            print('No non-zero values of {} to export!'.format(varname))
            return
        # rename the columns to indicate stress periods
        df.columns = ['{}{}'.format(i, varname) for i in range(df.shape[1])]
        segs = df.index

        # join the pivoted values to reach location info
        # for now, follow mf2005 model and assume that variable applies to reach 1
        isseg = np.array([True if s in segs else False for s in self.reach_data.iseg])
        locations = isseg & (self.reach_data.ireach == 1)
        rd = self.reach_data.loc[locations][['node', 'k', 'i', 'j', 'iseg', 'ireach']].copy()
        rd.sort_values(by=['iseg'], inplace=True)
        rd.index = rd.iseg
        assert np.array_equal(rd.index.values, df.index.values)
        rd = rd.join(df)

        export_reach_data(rd, self.grid, filename, geomtype='point')

    def export_period_data(self, filename=None, geomtype='point'):
        """Export point shapefile showing locations of period data
        in a MODFLOW-6 SFR package (e.g. inflows, runoff, etc.)

        Parameters
        ----------
        f : str, filename
        geomtype : str ('point' or 'polygon')
            write the locations as points at the cell centers, or polygons
            of the model cells containing the period data.

        """
        data = self.period_data.dropna(axis=1).sort_values(by=['per', 'rno'])
        nodes = dict(zip(self.reach_data.rno, self.reach_data.node))

        for var in ['evaporation', 'inflow', 'rainfall', 'runoff', 'stage']:
            if var in data.columns:
                # pivot the segment data to segments x periods with values of varname
                aggfunc = 'mean'  # how to aggregate multiple instances of rno/per combinations
                if var in ['inflow', 'runoff']:
                    aggfunc = 'sum'
                df = data.reset_index(drop=True).pivot_table(index='rno', columns='per', values=var,
                                                             aggfunc=aggfunc).reset_index()
                # rename the columns to indicate stress periods
                df.columns = ['rno'] + ['{}{}'.format(i, var) for i in range(df.shape[1] - 1)]
                df['node'] = [nodes[rno] for rno in df['rno']]
                if filename is None:
                    filename = self.package_name + '_sfr_period_data_{}.shp'.format(var)
                elif var not in filename:
                    filename = os.path.splitext(filename)[0] + '_{}.shp'.format(var)
                export_reach_data(df, self.grid, filename, geomtype=geomtype)

    def export_observations(self, filename=None, geomtype='point'):

        data = self.observations
        nodes = dict(zip(self.reach_data.rno, self.reach_data.node))
        data['node'] = [nodes[rno] for rno in data['rno']]
        if filename is None:
            filename = self.observations_file + '.shp'
        export_reach_data(data, self.grid, filename, geomtype=geomtype)