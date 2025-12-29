import os
from pathlib import Path
import time
from packaging import version
import warnings
import yaml
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterstats import zonal_stats
from shapely.geometry import LineString, Point
from gisutils import df2shp, get_authority_crs
from sfrmaker.routing import find_path, renumber_segments
from sfrmaker.checks import valid_rnos, valid_nsegs, rno_nseg_routing_consistent
from sfrmaker.elevations import get_slopes, smooth_elevations, sample_reach_elevations
from sfrmaker.flows import add_to_perioddata, add_to_segment_data
from sfrmaker.gis import export_reach_data, project
from sfrmaker.observations import write_gage_package, write_mf6_sfr_obsfile, add_observations
from sfrmaker.units import convert_length_units, itmuni_values, lenuni_values
from sfrmaker.utils import get_sfr_package_format, get_input_arguments, assign_layers, update
import sfrmaker
from sfrmaker.base import DataPackage
from sfrmaker.mf5to6 import segment_data_to_period_data
from sfrmaker.observations import get_closest_reach
from sfrmaker.reaches import consolidate_reach_conductances
from sfrmaker.rivdata import RivData

try:
    import flopy
    fm = flopy.modflow
    mf6 = flopy.mf6
except:
    flopy = False


class SFRData(DataPackage):
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

    dtypes = {'rno': int, 'node': int, 'k': int, 'i': int, 'j': int,
              'iseg': int, 'ireach': int, 'outreach': int, 'line_id': 'object',
              'per': int, 'nseg': int, 'icalc': int, 'outseg': int,
              'iupseg': int, 'iprior': int, 'nstrpts': int,
              'name': object, 'geometry': object}

    # LENUNI = {"u": 0, "f": 1, "m": 2, "c": 3}
    len_const = {0: 1.0, 1: 1.486, 2: 1.0, 3: 100.}
    # {"u": 0, "s": 1, "m": 2, "h": 3, "d": 4, "y": 5}
    time_const = {1: 1., 2: 60., 3: 3600., 4: 86400., 5: 31557600.}

    # default values
    defaults = {'icalc': 1,
                'roughch': 0.037,
                'strthick': 1,
                'strhc1': 1,
                'istcb2': 223,
                'gage_starting_unit_number': 250
                }
    package_type = 'sfr'

    def __init__(self, reach_data=None,
                 segment_data=None, grid=None,
                 model=None,
                 isfr=None,
                 model_length_units="undefined", model_time_units='days',
                 enforce_increasing_nsegs=True,
                 default_slope=0.001, minimum_slope=0.0001,
                 maximum_slope=1.,
                 package_name='model',
                 **kwargs):
        DataPackage.__init__(self, grid=grid, model=model, isfr=isfr,
                         model_length_units=model_length_units,
                         model_time_units=model_time_units,
                         package_name=package_name)

        # attributes
        self._period_data = None
        self._observations = None
        self._observations_filename = None
        self._minimum_slope = minimum_slope
        self._default_slope = default_slope
        self._maximum_slope = maximum_slope

        # convert any modflow6 kwargs to modflow5
        kwargs = {SFRData.mf5names[k] if k in SFRData.mf6names else k:
                  v for k, v in kwargs.items()}
        # update default values (used in segment and reach data setup)
        self.defaults.update(kwargs)

        self.reach_data = self._setup_reach_data(reach_data)
        self.segment_data = self._setup_segment_data(segment_data)
        self.isfropt0_to_1()  # distribute any isfropt=0 segment data to reaches

        # routing
        self._segment_routing = None  # dictionary of routing connections
        self._rno_routing = None  # dictionary of rno routing connections
        self._paths = None  # routing sequence from each segment to outlet
        self._reach_paths = None  # routing sequence from each reach number to outlet

        if not self._valid_nsegs(increasing=enforce_increasing_nsegs):
            self.reset_segments()

        # establish rno routing
        # set_outreaches checks for valid rnos and resets them if not
        # resets out reaches either way using segment data
        # not the ideal logic for MODFLOW 6 case where only
        # rno and connections are supplied
        self.set_outreaches()
        self.update_slopes(default_slope=default_slope, minimum_slope=minimum_slope,
                        maximum_slope=maximum_slope)

        # have to set the model last, because it also sets up a flopy sfr package instance
        self.model = model  # attached flopy model instance
        # self._ModflowSfr2 = None # attached instance of flopy modflow_sfr2 package object

        # MODFLOW-2005 gages will be assigned sequential unit numbers
        # starting at gage_starting_unit_number
        self.gage_starting_unit_number = self.defaults['gage_starting_unit_number']

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
        if self._package_name is None:
            if self.model is not None:
                self._package_name = self.model.name
            else:
                self._package_name = 'model'
        return self._package_name

    @package_name.setter
    def package_name(self, package_name):
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
        """A `flopy.modflow.mfsfr2.ModflowSfr2` represenation of the sfr dataset."""
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
            values = segment_data[c].astype(SFRData.dtypes.get(c, np.float32))
            # fill any nan values with 0 (same as empty segment_data;
            # for example elevation if it wasn't specified and
            # will be sampled from the DEM)
            values.fillna(0, inplace=True)
            sd[c] = values

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
            self._observations = pd.DataFrame(columns=['obsname', 'obstype', 'rno', 'iseg', 'ireach'])
        return self._observations

    @property
    def observations_file(self):
        if self._observations_filename is None:
            if self.model is not None and self.model.version == 'mf6':
                self._observations_filename = self.package_name + '.sfr.obs'
            else:
                self._observations_filename = self.package_name + '.gage'
        return self._observations_filename

    @observations_file.setter
    def observations_file(self, observations_filename):
        self._observations_filename = observations_filename

    @property
    def period_data(self):
        if self._period_data is None:
            self._period_data = self._get_period_data()
        # index by period and rno
        if not self._period_data.index.names == ['per', 'rno']:
            self._period_data.set_index(['per', 'rno'], inplace=True)
        return self._period_data

    def _get_period_data(self):
        print('converting segment data to period data...')
        return segment_data_to_period_data(self.segment_data, self.reach_data)

    def add_to_perioddata(self, data, flowline_routing=None,
                          variable='inflow',
                          line_id_column=None,
                          rno_column=None,
                          period_column='per',
                          data_column='Q_avg',
                          one_inflow_per_path=False, distribute_flows_to_reaches=False):
        return add_to_perioddata(self, data, flowline_routing=flowline_routing,
                                 variable=variable,
                                 line_id_column=line_id_column,
                                 rno_column=rno_column,
                                 period_column=period_column,
                                 data_column=data_column,
                                 one_inflow_per_path=one_inflow_per_path,
                                 distribute_flows_to_reaches=distribute_flows_to_reaches)

    def add_to_segment_data(self, data, flowline_routing,
                            variable='flow',
                            line_id_column=None,
                            segment_column='segment',
                            period_column='per',
                            data_column='Q_avg',
                            one_inflow_per_path=False
                            ):
        return add_to_segment_data(self, data, flowline_routing,
                                   variable=variable,
                                   line_id_column=line_id_column,
                                   segment_column=segment_column,
                                   period_column=period_column,
                                   data_column=data_column,
                                   one_inflow_per_path=one_inflow_per_path)
    @property
    def paths(self):
        """Dict listing routing sequence for each segment
        in SFR network."""
        if self._paths is None:
            self._set_paths()
            return self._paths
        if self._routing_changed():
            self._reset_routing()
        return self._paths

    @property
    def reach_paths(self):
        """Dict listing routing sequence for each segment
        in SFR network."""
        if self._paths is None:
            self._set_reach_paths()
            return self._reach_paths
        if self._routing_changed():
            self._reset_routing()
        return self._reach_paths

    def _set_paths(self):
        routing = self.segment_routing
        self._paths = {seg: find_path(routing, seg) for seg in routing.keys()}

    def _set_reach_paths(self):
        routing = self.rno_routing
        self._reach_paths = {rno: find_path(routing, rno) for rno in routing.keys()}

    def _reset_routing(self):
        self.reset_reaches()
        self.reset_segments()
        self._segment_routing = None
        self._set_paths()
        self._rno_routing = None
        self._set_reach_paths()

    def _routing_changed(self):
        sd = self.segment_data.groupby('per').get_group(0)
        rd = self.reach_data
        # check if segment routing in dataframe is consistent with routing dict
        segment_routing = dict(zip(sd.nseg, sd.outseg))
        segment_routing_changed = segment_routing != self._segment_routing

        # check if reach routing in dataframe is consistent with routing dict
        reach_routing = dict(zip(rd.rno, rd.outreach))
        reach_routing_changed = reach_routing != self._rno_routing

        # check if segment and reach routing in dataframe are consistent
        consistent = rno_nseg_routing_consistent(sd.nseg, sd.outseg,
                                                 rd.iseg, rd.ireach,
                                                 rd.rno, rd.outreach)
        # return True if the dataframes changed,
        # or are inconsistent between segments and reach numbers
        return segment_routing_changed & reach_routing_changed & ~consistent

    def repair_outsegs(self):
        """Set any outsegs that are not nsegs or lakes to 0 (outlet status)"""
        isasegment = np.isin(self.segment_data.outseg,
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
        try:
            self.reach_data['ireach'] = ireach
        except:
            j=2

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

    def assign_layers(self, adjusted_botm_output_path='.'):
        """Assign model layers to SFR reaches, using the discretzation
        package in the attached model. New botm elevations for the model
        will be written to a text array file, if any streambed bottoms
        are below the model bottom.

        Parameters
        ----------
        adjusted_botm_output_path : str
            Path for writing the text array of adjusted model bottom
            elevations, by default, '.'
        """
        if self.model is not None and hasattr(self.model, 'dis'):
            botm = self.model.dis.botm.array.copy()
            idomain = None
            if self.model.version == 'mf6':
                idomain = self.model.dis.idomain.array.copy()
            else:
                bas6 = getattr(self.model, 'bas6')
                if bas6 is not None:
                    idomain = bas6.ibound.array
                    
            nlay = botm.shape[0] + 1
            layers, new_botm = assign_layers(
                self.reach_data, botm_array=botm, idomain=idomain)
            self.reach_data['k'] = layers
            if new_botm is not None:
                new_bottom_is_same = np.allclose(botm.ravel(), 
                                                 new_botm.ravel(), rtol=0.0001)
                if not new_bottom_is_same:
                    outfiles = []
                    if len(new_botm.shape) == 2:
                        write_layers = [nlay - 1]
                    else:
                        write_layers = list(range(nlay))
                    outpath = Path(adjusted_botm_output_path)
                    for k in write_layers:
                        outfile = outpath / f'{self.package_name}_layer_{k}_new_botm_elevations.dat'
                        np.savetxt(outfile, new_botm[-1], fmt='%.2f')
                        outfiles.append(outfile)
                    msg = ('Sfrmaker pushed some model botm elevations downward '
                        'to accomodate streambed bottoms. New botm elevations '
                        'written to:\n')
                    for f in outfiles:
                        msg += f'{f}\n'
                    print(msg)
        else:
            print('Need a model instance with a discretization package to assign layers. '
                  'A model can be assigned to the SFRData.model attribute.')

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
            m = fm.Modflow(modelname=self.package_name, model_ws=model_ws,
                           structured=self.structured)
            nper = self.segment_data['per'].max() + 1
            if nper > 1:
                dis = fm.ModflowDis(m, nrow=self.grid.nrow, ncol=self.grid.ncol, 
                                    nlay=self.grid.nlay, nper=nper)
                
            if model is not None:
                nper = 1
                if 'tdis' in model.simulation.package_key_dict.keys():
                    tdis = model.simulation.tdis
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
                tdis = model.simulation.tdis
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

        # create an sfrmaker.Mf6SFR instance
        from .mf5to6 import Mf6SFR

        # update the ModflowSfr2 instance in case the SFR dataset has changed
        self.create_modflow_sfr2()
        sfr6 = Mf6SFR(self.modflow_sfr2)

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
            for record in connectiondata:
                zero_based_record = []
                for rno in record:
                    if rno is not None:
                        if rno > 0:
                            rno -= 1
                        elif rno < 0:
                            rno += 1
                    zero_based_record.append(rno)
                #zero_based_record = tuple(i - 1 if i > 0 else i + 1 for i in record)
                #if len(zero_based_record) == 1:
                #    zero_based_record.append(None)
                zero_based_connectiondata.append(zero_based_record)
            connectiondata = zero_based_connectiondata
        assert packagedata['rno'].min() == 0
        #assert np.min(list(map(np.min, map(np.abs, connectiondata)))) < 1

        # set cellids to None for unconnected reaches or where idomain == 0
        # can only do this with flopy versions 3.3.1 and later, otherwise flopy will bomb
        if version.parse(flopy.__version__) > version.parse('3.3.0'):
            unconnected = ~packagedata['rno'].isin(np.array(list(sfr6.connections.keys())) - 1).values
            inactive = m.dis.idomain.array[packagedata.k.values,
                                           packagedata.i.values,
                                           packagedata.j.values] != 1
            packagedata.loc[unconnected | inactive, 'cellid'] = 'none'
        packagedata = packagedata[columns].values.tolist()

        period_data = None
        if sfr6.period_data is not None:
            # TODO: add method to convert period_data df to MF6 input
            # raise NotImplemented("Support for mf6 period_data input not implemented yet. "
            #                     "Use sfrdata.write_package(version='mf6') instead.")
            pass
        auxiliary = None
        if 'auxiliary' in sfr6.options_block:
            auxiliary = sfr6.options_block.split('auxiliary')[1].strip().split('\n')[0].split()
        mf6sfr = mf6.ModflowGwfsfr(model=m, unit_conversion=unit_conversion,
                                   stage_filerecord=stage_filerecord,
                                   budget_filerecord=budget_filerecord,
                                   nreaches=len(self.reach_data),
                                   packagedata=packagedata,
                                   connectiondata=connectiondata,
                                   auxiliary=auxiliary,
                                   diversions=None,  # TODO: add support for diversions
                                   perioddata=period_data,  # TODO: add support for creating mf6 perioddata input
                                   **kwargs)
        return mf6sfr

    def add_observations(self, data, flowline_routing=None,
                         obstype=None, sfrlines_shapefile=None,
                         x_location_column=None,
                         y_location_column=None,
                         line_id_column=None,
                         rno_column=None,
                         obstype_column=None,
                         obsname_column='site_no',
                         gage_starting_unit_number=250):
        self.gage_starting_unit_number = gage_starting_unit_number
        added_obs = add_observations(self, data, flowline_routing=flowline_routing,
                                     obstype=obstype, sfrlines_shapefile=sfrlines_shapefile,
                                     x_location_column=x_location_column,
                                     y_location_column=y_location_column,
                                     line_id_column=line_id_column,
                                     rno_column=rno_column,
                                     obstype_column=obstype_column,
                                     obsname_column=obsname_column)

        # replace any observations that area already in the observations table
        if isinstance(self._observations, pd.DataFrame):
            existing_obs = set(zip(self._observations['obsname'], self._observations['obstype']))
            new_obs = set(zip(added_obs['obsname'], added_obs['obstype']))
            exists_already = {obs[0] for obs in existing_obs.intersection(new_obs)}
            exists_already = self._observations['obsname'].isin(exists_already)
            self._observations = self._observations.loc[~exists_already]
        self._observations = pd.concat([self.observations, added_obs], axis=0).reset_index(drop=True)

        for df in self._observations, added_obs:
            # enforce dtypes (pandas doesn't allow an empty dataframe
            # to be initialized with more than one specified dtype)
            for col in ['rno', 'iseg', 'ireach']:
                df[col] = df[col].astype(int)
        return added_obs

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
        from sfrmaker.reaches import interpolate_to_reaches

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

    def sample_reach_elevations(self, dem,
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
        dem : path to valid raster dataset
            Must be in same Coordinate Reference System as model grid.
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
        return sample_reach_elevations(self.reach_data, dem,
                            sfrmaker_modelgrid=self.grid,
                            model_crs=self.crs,
                            method=method,
                            buffer_distance=buffer_distance,
                            smooth=smooth,
                            elevation_data=elevation_data,
                            elevation_data_crs=elevation_data_crs,
                            elevation_data_layer=elevation_data_layer,
                            elevation_data_errors=elevation_data_errors)


    def set_streambed_top_elevations_from_dem(self, filename, elevation_units=None,
                                              dem=None,
                                              method='buffers',
                                              **kwargs):
        """Set streambed top elevations from a DEM raster.
        Runs sfrdata.sample_reach_elevations

        Parameters
        ----------
        filename : path to valid raster dataset
            Must be in same Coordinate Reference System as model grid.
        elevation_units : str
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
        if dem is not None:
            warnings.warn('set_streambed_top_elevations_from_dem: dem argument is deprecated. '
                          'Use filename instead.', DeprecationWarning)
            filename = dem
        if 'dem_z_units' in kwargs:
            warnings.warn('set_streambed_top_elevations_from_dem: dem_z_units argument is deprecated. '
                          'Use elevation_units instead.', DeprecationWarning)
            elevation_units = kwargs['dem_z_units']
        sampled_elevs = self.sample_reach_elevations(dem=filename, method=method,
                                                     **kwargs)
        if elevation_units is None:
            elevation_units = self.model_length_units
        mult = convert_length_units(elevation_units, self.model_length_units)
        self.reach_data['strtop'] = [sampled_elevs[rno]
                                     for rno in self.reach_data['rno'].values]
        self.reach_data['strtop'] *= mult
        # update the slopes as well
        self.update_slopes()

    def update_slopes(self, default_slope=None, minimum_slope=None,
                   maximum_slope=None):
        """Compute slopes by reach using values in strtop (streambed top) and rchlen (reach length)
        columns of reach_data. The slope for a reach n is computed as strtop(n) - strtop(n+1) / rchlen(n).
        Slopes for outlet reaches are set equal to a default value (default_slope).
        Populates the slope column in reach_data.

        Parameters
        ----------
        default_slope : float
            Slope value applied to outlet reaches (where water leaves the model).
            By default None, in which case '_default_slope' attribute is used.
        minimum_slope : float
            Assigned to reaches with computed slopes less than this value.
            This ensures that the Manning's equation won't produce unreasonable values of stage
            (in other words, that stage is consistent with assumption that
            streamflow is primarily drive by the streambed gradient).
            By default None, in which case '_minimum_slope' attribute is used.
        maximum_slope : float
            Assigned to reaches with computed slopes more than this value.
            By default None, in which case '_maximum_slope' attribute is used.
            
        Notes
        -----
        If default_slope, minimum_slope or maximum_slope are not None, the
        respective attribute (for example '_default_slope') will be set with the supplied value.
        """
        if default_slope is not None:
            self._default_slope = default_slope
        if minimum_slope is not None:
            self._minimum_slope = minimum_slope
        if maximum_slope is not None:
            self._maximum_slope = maximum_slope
        rd = self.reach_data
        slopes = get_slopes(rd['strtop'], rd['rchlen'], rd['rno'], rd['outreach'], 
                            default_slope=self._default_slope, 
                            minimum_slope=self._minimum_slope,
                            maximum_slope=self._maximum_slope)
        self.reach_data['slope'] = slopes

    @classmethod
    def from_package(cls, sfrpackagefile, grid, namefile=None,
                     sim_name=None, model_ws='.',
                     version=None, model_name='model', package_name=None,
                     linework=None):
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
        # todo:  finish SFRData.from_package
        raise NotImplementedError('SFRData.from_package not implemented yet.')

        if version is None:
            version = get_sfr_package_format(sfrpackagefile)
        if package_name is None:
            package_name, _ = os.path.splitext(sfrpackagefile)
        # load the model and SFR package
        if namefile is not None:
            model_ws, namefile = os.path.split(namefile)
            m = fm.Modflow.load(namefile, model_ws=model_ws, load_only=['SFR'])
        elif sim_name is not None:
            sim = flopy.mf6.MFSimulation.load(sim_name, 'mf6', 'mf6', model_ws)
            m = sim.get_model(model_name)
        else:
            if version != 'mf6':
                m = fm.Modflow(model_ws=model_ws)
                sfr = fm.ModflowSfr2.load(sfrpackagefile, model=m)
            else:
                sim = flopy.mf6.MFSimulation(sim_ws=model_ws)
                m = flopy.mf6.ModflowGwf(sim, modelname=model_name,
                                         model_nam_file='{}.nam'.format(model_name))
                sfr = mf6.ModflowGwfsfr.load(sfrpackagefile, model=m)

        if m.version != 'mf6':
            reach_data = pd.DataFrame(m.sfr.reach_data)
            perioddata_list = []
            for per, spd in m.sfr.segment_data.items():
                if spd is not None:
                    spd = spd.copy()
                    spd['per'] = per
                    perioddata_list.append(spd)
            segment_data = pd.concat(perioddata_list)
        else:
            pass
        return cls(reach_data=reach_data, segment_data=segment_data,
                   model=m,
                   grid=grid)

    @classmethod
    def from_tables(cls, reach_data, segment_data,
                    grid=None, isfr=None):
        reach_data = pd.read_csv(reach_data)
        segment_data = pd.read_csv(segment_data)
        return cls(reach_data=reach_data, segment_data=segment_data,
                   grid=grid, isfr=isfr)

    @classmethod
    def from_yaml(cls, config_file, package_name=None, output_path=None,
                  write_output=True):
        """Create an SFRData instance from a yaml-format configuration file.

        Parameters
        ----------
        config_file : str or mapping
            Path to Sfrmaker configuration file in yaml format,
            or a dictionary-like mapping.
        output_path : str
            Where output files will be saved. Default is '.'
            (in default_config.yml).
        package_name : str
            Optional argument for naming the output SFR package; otherwise
            can be specified in configuration file.
        write_output : bool
            Whether or not to write output files, including package file,
            reach and segment data tables and shapefiles.

        Returns
        -------
        sfrdata : sfrmaker.SFRData instance
        """
        # load the configuration file
        wd = os.getcwd()
        if isinstance(config_file, str) or isinstance(config_file, Path):
            # change the cwd to the config file path
            path, config_file = os.path.split(config_file)
            path = path if len(path) > 0 else '.'
            os.chdir(path)
            with open(config_file) as src:
                cfg = yaml.load(src, Loader=yaml.Loader)
        # or accept a mapping as input
        else:
            cfg = config_file

        # read in the default configuration
        defaults_file = os.path.join(os.path.split(__file__)[0],
                                     'default_config.yml')
        with open(defaults_file) as src:
            defaults = yaml.load(src, Loader=yaml.Loader)

        # add defaults to configuration
        #cfg = update(defaults, cfg)

        # set the package_name and output paths
        if package_name is not None:
            _, package_name = os.path.split(package_name)
            package_name, _ = os.path.splitext(package_name)
        elif 'package_name' in cfg:
            _, package_name = os.path.split(cfg['package_name'])
            package_name, _ = os.path.splitext(package_name)
        elif 'simulation' in cfg:
            package_name, _ = os.path.splitext(cfg['simulation']['sim_name'])
        elif 'model' in cfg:
            package_name, _ = os.path.splitext(cfg['model']['namefile'])
        else:
            package_name = 'model'
        if output_path is None:
            output_path = cfg.get('output_path', '.')
        if not os.path.isdir(output_path):
            os.makedirs(output_path)

        # model grids from shapefiles or flopy grids are supported
        grid = None
        if 'modelgrid' in cfg:
            grid_cfg = cfg['modelgrid']
            if 'shapefile' in grid_cfg:
                grid_kwargs = get_input_arguments(grid_cfg, sfrmaker.StructuredGrid.from_shapefile)
                grid = sfrmaker.StructuredGrid.from_shapefile(**grid_kwargs)
            elif not flopy:
                raise ImportError("Specifying a modelgrid with xoffset, yoffset, etc. requires flopy.")
            else:
                # convert delr/delc to ndarrays
                for delrc, ncells in {'delr': 'ncol', 'delc': 'nrow'}.items():
                    if np.isscalar(grid_cfg[delrc]):
                        if ncells not in grid_cfg:
                            msg = ("Scalar value for {} requires "
                                   "specification of {}.".format(delrc, ncells))
                            raise KeyError(msg)
                        grid_cfg[delrc] = np.array([grid_cfg[delrc]] * grid_cfg[ncells])
                grid_cfg['xoff'] = grid_cfg.get('xoffset', grid_cfg.get('xoff', 0.))
                grid_cfg['yoff'] = grid_cfg.get('yoffset', grid_cfg.get('yoff', 0.))
                grid_kwargs = get_input_arguments(grid_cfg, flopy.discretization.StructuredGrid)
                grid = flopy.discretization.StructuredGrid(**grid_kwargs)

        # create a Lines instance
        no_lines_msg = "A flowlines: block must be specified in configuration file."
        if 'flowlines' not in cfg:
            raise KeyError(no_lines_msg)
        lines_config = cfg['flowlines']
        # custom hydrography option
        if 'filename' in lines_config:
            lines_config['shapefile'] = lines_config['filename']
            lines_kwargs = get_input_arguments(lines_config, sfrmaker.Lines.from_shapefile)
            lines = sfrmaker.Lines.from_shapefile(**lines_kwargs)
        # nhdplus option
        else:
            # renames to deal with case issues
            renames = {'nhdplus_paths': 'NHDPlus_paths',
                       'nhdflowlines': 'NHDFlowlines',
                       'plusflowlinevaa': 'PlusFlowlineVAA',
                       'plusflow': 'PlusFlow'}
            lines_config = {renames.get(k, k): v for k, v in lines_config.items()}
            lines_kwargs = get_input_arguments(lines_config, sfrmaker.Lines.from_nhdplus_v2)
            lines = sfrmaker.Lines.from_nhdplus_v2(**lines_kwargs)

        # load a model if there is one
        # modflow 6
        model_version = cfg.get('package_version', 'mf6')
        no_flopy_msg = "Specifying a model requires flopy."
        if 'simulation' in cfg:
            if not flopy:
                raise ImportError(no_flopy_msg)
            sim_kwargs = get_input_arguments(cfg['simulation'], flopy.mf6.MFSimulation)
            sim = mf6.MFSimulation.load(**sim_kwargs)
            model = sim.get_model(model_name=cfg['model']['modelname'])
            model_version = model.version
            if grid is None:
                grid = model.modelgrid
        # modflow-2005
        elif 'model' in cfg:
            if not flopy:
                raise ImportError(no_flopy_msg)
            cfg['model']['f'] = cfg['model']['namefile']
            model_kwargs = get_input_arguments(cfg['model'], flopy.modflow.Modflow.load)
            model = flopy.modflow.Modflow.load(**model_kwargs)
            model_version = model.version
            if grid is None:
                grid = model.modelgrid
        else:
            model = None

        # create an SFRData instance
        to_sfr_kwargs = cfg.copy()
        to_sfr_kwargs.update(cfg.get('options', {}))
        to_sfr_kwargs = get_input_arguments(to_sfr_kwargs, sfrmaker.Lines.to_sfr)
        if cfg.get('active_area') is not None:
            warnings.warn('the active_area: block is deprecated. '
                          'Use active_area: argument in options: block instead.', DeprecationWarning)
            to_sfr_kwargs['active_area'] = cfg.get('active_area', {}).get('filename')
        to_sfr_kwargs['model'] = model
        to_sfr_kwargs['grid'] = grid
        to_sfr_kwargs['package_name'] = package_name
        # to_sfr() populates isfr from model if neither active_area or isfr are argued
        sfrdata = lines.to_sfr(**to_sfr_kwargs)
        sfrdata._tables_path = os.path.join(output_path, sfrdata._tables_path)
        sfrdata._shapefiles_path = os.path.join(output_path, sfrdata._shapefiles_path)

        # setup elevations
        if cfg.get('options', cfg).get('set_streambed_top_elevations_from_dem', False):
            warnings.warn(('The set_streambed_top_elevations_from_dem argument '
                           'is deprecated. This option is now activated by including a dem: '
                           'block in the configuration file.'), DeprecationWarning)
        if 'dem' in cfg:
            error_msg = ("set_streambed_top_elevations_from_dem=True "
                         "requires a dem: block.")
            #if 'dem' not in cfg:
            #    raise KeyError(error_msg)
            # get valid arguments for set_streambed_top_elevations_from_dem
            # and sample_reach_elevations,
            # which accepts kwargs from set_streambed_top_elevations_from_dem
            dem_kwargs = get_input_arguments(cfg['dem'],
                                             sfrmaker.SFRData.sample_reach_elevations)
            dem_kwargs2 = get_input_arguments(cfg['dem'],
                                              sfrmaker.SFRData.set_streambed_top_elevations_from_dem)
            dem_kwargs.update(dem_kwargs2)
            sfrdata.set_streambed_top_elevations_from_dem(**dem_kwargs)
        else:
            sfrdata.reach_data['strtop'] = sfrdata.interpolate_to_reaches('elevup', 'elevdn')
            sfrdata.update_slopes()

        # assign layers to the sfr reaches
        if model is not None:
            sfrdata.assign_layers(adjusted_botm_output_path=output_path)

        # option to convert reaches to the River Package
        if 'to_riv' in cfg:
            to_riv_args = get_input_arguments(cfg['to_riv'], SFRData.to_riv)
            rivdata = sfrdata.to_riv(**to_riv_args)
            sfrdata.setup_riv(rivdata)
            rivdata.write_table()
            rivdata.write_shapefiles()
            print('to_riv option: output table and shapefiles written, but'
                  'writing of RIV package not implemented yet. Use flopy to'
                  'write the RIV package.')

        # add inflows
        if 'inflows' in cfg:
            inflows_input = cfg['inflows']
            # resample inflows to model stress periods
            inflows_input['id_column'] = inflows_input['line_id_column']
            inflows_by_stress_period = pd.read_csv(inflows_input['filename'])
            inflows_by_stress_period[inflows_input['id_column']] =\
                inflows_by_stress_period[inflows_input['id_column']].astype(int).astype(str)

            # check if all inflow sites are included in sfr network
            missing_sites = set(inflows_by_stress_period[inflows_input['id_column']]). \
                difference(lines._original_routing.keys())
            if any(missing_sites):
                inflows_routing_input = cfg.get('inflows_routing')
                if inflows_routing_input is None:
                    raise KeyError(('inflow sites {} are not within the model sfr network. '
                                    'Please supply an inflows_routing: block'.format(missing_sites)))
                routing = pd.read_csv(inflows_routing_input['filename'])
                routing = dict(zip(routing[inflows_routing_input['id_column']],
                                   routing[inflows_routing_input['routing_column']]))
            else:
                routing = lines._original_routing
            missing_sites = any(set(inflows_by_stress_period[inflows_input['id_column']]). \
                                difference(routing.keys())),
            if any(missing_sites):
                raise KeyError(('Inflow sites {} not found in {}'.format(missing_sites,
                                                                         inflows_routing_input['filename'])))

            # add resampled inflows to SFR package
            inflows_input['data'] = inflows_by_stress_period
            inflows_input['flowline_routing'] = routing
            if model_version == 'mf6':
                inflows_input['variable'] = 'inflow'
                method = sfrdata.add_to_perioddata
            else:
                method = sfrdata.add_to_segment_data
            kwargs = get_input_arguments(inflows_input.copy(), method)
            method(**kwargs)

        # add observations
        if 'observations' in cfg:
            if not model:
                pass #  raise KeyError('Setup of observations input results a model: block')
            key = 'filename' if 'filename' in cfg['observations'] else 'filenames'
            cfg['observations']['data'] = cfg['observations'][key]
            kwargs = get_input_arguments(cfg['observations'].copy(), sfrdata.add_observations)
            obsdata = sfrdata.add_observations(**kwargs)

        if write_output:
            if model is None:
                package_file_path = os.path.join(output_path, package_name + '.sfr')
            else:
                package_file_path = os.path.join(model.model_ws, package_name + '.sfr')
            sfrdata.write_package(package_file_path,
                                  version=model_version)
            sfrdata.write_tables()
            sfrdata.write_shapefiles()

        # change the cwd back
        os.chdir(wd)
        return sfrdata

    def to_riv(self, segments=None, rno=None, line_ids=None, drop_in_sfr=True):
        """Cast one or more reaches to a RivData instance,
        which can then be written as input to the MODFLOW RIV package.

        Parameters
        ----------
        segments : sequence of ints
            Convert the listed segments, and all downstream segments,
            to the RIV package. If None, all segments are converted (default).
        rno : sequence of ints
            Convert the listed reach numbers (nro), and all downstream reaches,
            to the RIV package. If None, all reaches are converted (default).
        line_ids : sequence of ints or strings
            Convert the segments corresponding to the line_ids 
            (unique identifiers for LineString features in the source hydrography), 
            and all downstream segments to the RIV package. If None, all reaches 
            are converted (default).
        drop_in_sfr : bool
            Whether or not to remove the converted segments from the SFR package
            (default True)

        Returns
        -------
        riv : SFRmaker.RivData instance
        """

        # get the downstream segments to be converted
        loc = np.array([True] * len(self.reach_data))
        if segments is None and rno is None and line_ids is None:
            #loc = slice(None, None)  # all reaches
            # if all reaches are being converted,
            # skip dropping the reaches from the sfr package
            drop_in_sfr = False
        if segments is not None:
            if np.isscalar(segments):
                segments = [segments]
            loc = loc & self.reach_data.iseg.isin(segments)
        if line_ids is not None:
            if np.isscalar(line_ids):
                line_ids = [line_ids]
            line_ids = [str(lid) for lid in line_ids]
            loc = loc & self.reach_data.line_id.isin(line_ids)
        if rno is not None:
            if np.isscalar(rno):
                rno = [rno]
            loc = loc & self.reach_data.rno.isin(rno)
        reaches = self.reach_data.loc[loc, 'rno'].tolist()
        to_riv_reaches = set()
        for rno in reaches:
            # skip reaches that have already been considered
            if rno not in to_riv_reaches:
                to_riv_reaches.update(set(self.reach_paths[rno]))

        # subset the RIV reaches from reach_data;
        # populate RIV input
        is_riv_reach = self.reach_data.rno.isin(to_riv_reaches)
        df = self.reach_data.loc[is_riv_reach].copy()

        # consolidate the reaches to 1 per cell
        df = consolidate_reach_conductances(df, keep_only_dominant=True)

        # add in a column for the stress period
        # even though transfer of any specified stage changes
        # aren't implemented yet
        df['per'] = 0

        df['cond'] = df['Cond_sum']

        if any(self.segment_data.depth1 > 0):
            # find period with most data
            per = np.argmax(self.segment_data.groupby('per').count().nseg)
            reach_depths = self.interpolate_to_reaches('depth1', 'depth2', per=per)
            df['stage'] = df['strtop'] + reach_depths
        else:
            df['stage'] = df['strtop']
        df['rbot'] = df['strtop'] - df['strthick']
        cols = ['per', 'rno', 'node', 'k', 'i', 'j', 'cond', 'stage', 'rbot',
                'outreach', 'asum', 'line_id', 'name', 'geometry']
        cols = [c for c in cols if c in df.columns]
        riv_data = df[cols].copy()

        # retain routing information between reaches;
        # but reset the numbering since the SFR number will be reset anways
        # (after the RIV reaches are removed from the SFR dataset)
        new_rnos = renumber_segments(riv_data['rno'], riv_data['outreach'])
        riv_data['rno'] = [new_rnos[outreach] for outreach in riv_data['rno']]
        riv_data['outreach'] = [new_rnos[outreach] for outreach in riv_data['outreach']]

        riv = RivData(stress_period_data=riv_data, grid=self.grid,
                      model=self.model, model_length_units=self.model_length_units,
                      model_time_units=self.model_time_units,
                      package_name=self.package_name,)

        if drop_in_sfr:
            riv_reaches = self.reach_data.rno.isin(to_riv_reaches)
            riv_nodes = self.reach_data.loc[riv_reaches, 'node']
            is_riv_node = self.reach_data.node.isin(riv_nodes)

            # drop sfr from all model cells with a riv reach
            self.reach_data = self.reach_data.loc[~is_riv_node].copy()

            # reset riv outreaches (that are no longer in SFR network)
            # to zero (exit for sfr network)
            riv_outreaches = set(self.reach_data.outreach).difference(self.reach_data.rno)
            outreach_is_riv = self.reach_data.outreach.isin(riv_outreaches)
            self.reach_data.loc[outreach_is_riv, 'outreach'] = 0

            # drop segments that were completely converted to riv
            riv_segments = set(self.segment_data.nseg).difference(set(self.reach_data.iseg))
            # reset outseg numbers in reach_data to zero
            outseg_is_riv = self.reach_data.outseg.isin(riv_segments)
            self.reach_data.loc[outseg_is_riv, 'outseg'] = 0
            # drop these segments from segment_data
            is_riv_segment = self.segment_data.nseg.isin(riv_segments)
            self.segment_data = self.segment_data.loc[~is_riv_segment].copy()
            outseg_is_riv = self.segment_data.outseg.isin(riv_segments)
            self.segment_data.loc[outseg_is_riv, 'outseg'] = 0
            # reset the numbering to be consecutive
            self._reset_routing()
        return riv

    def run_diagnostics(self, checkfile=None, **kwargs):
        """Run the Flopy SFR diagnostic suite.

        Parameters
        ----------
        checkfile : str
            Path of file to write results to, by default, results
            are written to {}_SFR.chk, where {} is the
            SFRData.package_name attribute.
        kwargs : keyword arguments to flopy.modflow.ModflowSfr2.check()

        """
        if flopy:
            print('\nRunning Flopy v. {} diagnostics...'.format(flopy.__version__))
            self.create_modflow_sfr2(model=self.model)
            if checkfile is None:
                checkfile = '{}_SFR.chk'.format(self.package_name)
            self._ModflowSfr2.check(checkfile, **kwargs)
        print('wrote {}'.format(checkfile))

    def write_package(self, filename=None, version='mf2005', idomain=None,
                      options=None, run_diagnostics=True,
                      write_observations_input=True,
                      external_files_path=None, gage_starting_unit_number=None,
                      **kwargs):
        """Write an SFR package input file.

        Parameters
        ----------
        filename : str, optional
            File name for SFR package. By default None, in which
            case the filename of the attached :py:class:`flopy.modflow.mfsfr2.ModflowSfr2` instance
            (``SFRData.modflow_sfr2.fn_path``) is used.
        version : str, optional, {'mf2005', 'mfnwt', 'mf6'}
            MODFLOW version for the SFR package, by default 'mf2005'
        idomain : ndarray, optional
            3D numpy array designating active cells (idomain==1).
            SFR reaches in inactive cells will be written with 'none' in the cellid field.
            by default None
        options : list, optional
            List of strings to write to the MODFLOW-6 SFR options block. For example::

                options=['save_flows',
                         'BUDGET FILEOUT model.sfr.cbc',
                         'STAGE FILEOUT model.sfr.stage.bin']

            * An appropriate unit_conversion is written by default.
            * Cell-by-cell budget and stage output files are also added by default if not specified.
            * If the :attr:`SFRData.observations` is populated and `write_observations_input=True`,
              the observation input file is also added.
            * All files added by default are placed in the same path as the SFR package file. If another
              location is desired, the file(s) need to be specified explicitly as strings in the options
              list.

            See MODFLOW-6 documentation for other options.
            By default None.
        run_diagnostics : bool, optional
            Option to run the :ref:`diagnostic checks <Running diagnostics>`.
            by default True
        write_observations_input : bool, optional
            Option to write input to the the MODFLOW-2005 gage package or MODFLOW-6 Observations Process.
            Requires attached observations, as added through the :py:meth:`SFRData.add_observations` method.
            by default True
        external_files_path : str, optional
            Path for writing an external file for packagedata, relative to the location of the SFR package file.
            If specified, an open/close statement referencing the file is written to the packagedata block.
            By default, None (packagedata table is written to the SFR package file)
        gage_starting_unit_number : int, optional
            Starting unit number for gage output files, 
            by default None
        """        
        print('SFRmaker v. {}'.format(sfrmaker.__version__))
        # run the flopy SFR diagnostics
        if run_diagnostics:
            self.run_diagnostics()

        # recreate the flopy package object in case it changed
        self.create_modflow_sfr2(model=self.model, **kwargs)
        header_txt = "#SFR package created by SFRmaker v. {}, " \
                     "via FloPy v. {}\n".format(sfrmaker.__version__, flopy.__version__)
        header_txt += "#model length units: {}, model time units: {}".format(self.model_length_units,
                                                                             self.model_time_units)
        self.modflow_sfr2.heading = header_txt

        if filename is None:
            filename = self.modflow_sfr2.fn_path

        if gage_starting_unit_number is None:
            gage_starting_unit_number = self.gage_starting_unit_number
        if version in {'mf2005', 'mfnwt'}:

            if write_observations_input and len(self.observations) > 0:
                gage_package_filename = os.path.splitext(filename)[0] + '.gage'
                self.write_gage_package(filename=gage_package_filename,
                                        gage_starting_unit_number=gage_starting_unit_number)

            self.modflow_sfr2.write_file(filename=filename)
            print('wrote {}.'.format(filename))

        elif version == 'mf6':

            # instantiate Mf6SFR converter object with mf-nwt model/sfr package from flopy
            from .mf5to6 import Mf6SFR
            if options is None:
                # save budget and stage output by default
                # if options weren't specified, assume the output files
                # should be in the same location as the SFR package file
                # (do the same for obs below)
                _, just_the_filename = os.path.split(filename)
                options = ['save_flows',
                           'BUDGET FILEOUT {}.cbc'.format(just_the_filename),
                           'STAGE FILEOUT {}.stage.bin'.format(just_the_filename),
                           ]
            if write_observations_input and len(self.observations) > 0:
                if 'obs6 filein' not in ''.join(options).lower():
                    obs_input_filename = filename + '.obs'
                    options.append('OBS6 FILEIN {}'.format(just_the_filename + '.obs'))
                else:
                    for entry in options:
                        if 'obs6 filein' in entry.lower():
                            break
                    _, _, obs_input_filename = entry.split()
                self.write_mf6_sfr_obsfile(filename=obs_input_filename)


            sfr6 = Mf6SFR(SFRData=self, period_data=self.period_data,
                          idomain=idomain,
                          options=options)

            # write a MODFLOW 6 file
            sfr6.write_file(filename=filename, external_files_path=external_files_path)

    def write_tables(self, basename=None):
        """Write :py:attr:`~SFRData.reach_data`, :py:attr:`~SFRData.segment_data`,
        and :py:attr:`~SFRData.period_data` (if populated) to csv files. 

        Parameters
        ----------
        basename : str, optional
            Base name for csv files, by default None, in which case
            :py:attr:`SFRData.package_name` is used.
        """
        if basename is None:
            output_path = self._tables_path
            if not os.path.isdir(output_path):
                os.makedirs(output_path)
            basename = self.package_name
        else:
            output_path, basename = os.path.split(basename)
            basename, _ = os.path.splitext(basename)
        reach_data_file = os.path.normpath('{}/{}_sfr_reach_data.csv'.format(output_path, basename))
        self.reach_data.drop('geometry', axis=1).to_csv(reach_data_file, index=False)
        print('wrote {}'.format(reach_data_file))
        segment_data_file = os.path.normpath('{}/{}_sfr_segment_data.csv'.format(output_path, basename))
        self.segment_data.to_csv(segment_data_file, index=False)
        print('wrote {}'.format(segment_data_file))

        if self.period_data is not None and len(self.period_data) > 0:
            pd_file = os.path.normpath('{}/{}_sfr_period_data.csv'.format(output_path, basename))
            self.period_data.dropna(axis=1, how='all').to_csv(pd_file)
            print('wrote {}'.format(pd_file))

    def write_gage_package(self, filename=None, sitename_col='obsname',
                           gage_package_unit=25,
                           gage_starting_unit_number=None):
        """Write observation input for the MODFLOW-2005 Gage Package.

        Parameters
        ----------
        filename : str, optional
            Gage package file, by default None, in which case
            :py:attr:`SFRData.observations_file` is used.
        sitename_col : str
            Unique name or number for each gage site.
            By default, 'obsname', which is the default field in the 
            :attr:`SFRData.observations` table containing either 
            site numbers or site number-variable pairs.
        gage_package_unit : int, optional
            Unit number for Gage Package, by default 25
        gage_starting_unit_number : int, optional
            Starting unit number for gage output files, 
            by default None

        Returns
        -------
        gag : :py:class:`flopy.modflow.mfgage.ModflowGage` instance
        """        
        if filename is None:
            filename = self.observations_file
        else:
            self.observations_file = filename
        if gage_starting_unit_number is None:
            gage_starting_unit_number = self.gage_starting_unit_number
        gage_namfile_entries_file = filename + '.namefile_entries'
        model = self.model
        if model is None:
            model = self.modflow_sfr2.parent
            model.model_ws = os.path.split(filename)[0]
        return write_gage_package(self.observations,
                                  gage_package_filename=filename,
                                  gage_namfile_entries_file=gage_namfile_entries_file,
                                  model=model,
                                  sitename_col=sitename_col,
                                  gage_package_unit=gage_package_unit,
                                  start_gage_unit=gage_starting_unit_number)

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

    def export_outlets(self, filename=None):
        """Export shapefile of model cells with stream reaches."""
        if filename is None:
            filename = self.package_name + '_sfr_outlets.shp'
        nodes = self.reach_data.loc[self.reach_data.outreach == 0, 'node'].values
        export_reach_data(self.reach_data, self.grid, filename,
                          nodes=nodes, geomtype='point')

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
        df2shp(rd, filename, crs=self.grid.crs)

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

    def export_observations(self, filename=None, geomtype='point'):

        data = self.observations
        if len(data) == 0:
            print('No observations to export!')
            return
        nodes = dict(zip(self.reach_data.rno, self.reach_data.node))
        data['node'] = [nodes[rno] for rno in data['rno']]
        if filename is None:
            filename = self.observations_file + '.shp'
        export_reach_data(data, self.grid, filename, geomtype=geomtype)

