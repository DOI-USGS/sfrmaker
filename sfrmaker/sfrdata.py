import sys
sys.path.append('/Users/aleaf/Documents/GitHub/flopy3')
import numpy as np
import pandas as pd
import flopy
fm = flopy.modflow

class sfrdata:
    """Class for working with a streamflow routing (SFR) dataset,
    where the stream network is discretized into reaches contained
    within individual model cells. Reaches may be grouped into segments,
    with routing between segments specified, and routing between reaches
    within segments based on consecutive numbering (as in MODFLOW-2005).
    In this case, unique identifier numbers will be assigned to each reach
    (in the reachID column of the reach_data table), and routing connections
    between reachIDs will be computed. Alternatively, reaches and their
    routing connections can be specified directly, as in MODFLOW-6. In this
    case, MODFLOW-2005 input will be written with one reach per segment.

    Parameters
    ----------
    reach_data : DataFrame
        Table containing information on the SFR reaches.
    segment_data : DataFrame
        Table containing information on the segments (optional).
    grid : sfrmaker.grid class instance
    """
    # conversions from reach_data columns to
    # flopy (ModflowSfr2.reach_data) column names
    mf5names = {'length': 'rchlen'}

    # conversions to MODFLOW6 variable names
    mf6names = {'reachID': 'rno',
                'node': 'cellid',
                'length': 'rlen',
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

    # order for columns in reach_data
    rdcols = ['reachID', 'node', 'k', 'i', 'j',
              'segment', 'reach', 'length', 'width', 'slope',
              'strtop', 'strthick', 'strhc1', 'roughch',
              'thts', 'thti', 'eps', 'uhc',
              'outreach', 'outseg', 'line_id',
              'asum', 'geometry']
    # order for columns in segment_data
    sdcols = ['nper', 'nseg', 'icalc', 'outseg', 'iupseg',
              'iprior', 'nstrpts',
              'flow', 'runoff', 'etsw', 'pptsw',
              'roughch', 'roughbk', 'cdpth', 'fdpth',
              'awdth', 'bwdth',
              'hcond1', 'thickm1', 'elevup', 'width1', 'depth1',
              'thts1', 'thti1', 'eps1', 'uhc1',
              'hcond2', 'thickm2', 'elevdn', 'width2', 'depth2',
              'thts2', 'thti2', 'eps2', 'uhc2']

    def __init__(self, reach_data=None,
                 segment_data=None, grid=None, **kwargs):

        self.reach_data = reach_data
        self.segment_data = segment_data
        self.grid = grid
        self.structured = self.grid.structured

        # assign and default values from kwargs
        for k, v in kwargs.items():
            if k in self.rdcols:
                self.reach_data[k] = kwargs.pop(k)
            elif k in self.sdcols:
                self.segment_data[k] = kwargs.pop(k)

        # establish reachID routing
        self.set_outreaches()

    @property
    def model_units(self):
        """Computational lengths units of numerical model."""
        return self.grid.model_units

    @property
    def crs_units(self):
        """Length units of the coordinate reference system"""
        return self.grid.crs.length_units

    @property
    def segment_routing(self):
        sd = self.segment_data.groupby('per').get_group(0)
        graph = dict(
            zip(sd.nseg, sd.outseg))
        outlets = set(graph.values()).difference(
            set(graph.keys()))  # including lakes
        graph.update({o: 0 for o in outlets})
        return graph

    def repair_outsegs(self):
        isasegment = np.in1d(self.segment_data.outseg,
                             self.segment_data.nseg)
        isasegment = isasegment | (self.segment_data.outseg < 0)
        self.segment_data.loc[~isasegment, 'outseg'] = 0.

    def reset_reaches(self):
        """Ensure that the reaches in each segment are numbered
        consecutively starting at 1."""
        self.reach_data.sort_values(by=['iseg', 'ireach'])
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
        self.reach_data.sort_values(by=['iseg', 'ireach'])
        if not self._valid_reachIDs():
            self.reach_data['reachID'] = np.arange(1, len(self.reach_data) + 1)
        self.reset_reaches()  # ensure that each segment starts with reach 1
        self.repair_outsegs()  # ensure that all outsegs are segments, outlets, or negative (lakes)
        rd = self.reach_data
        outseg = self.segment_routing
        reach1IDs = dict(zip(rd[rd.ireach == 1].iseg,
                             rd[rd.ireach == 1].reachID))
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
                nextrchid = rd.reachID[i + 1]
            outreach.append(nextrchid)
        self.reach_data['outreach'] = outreach

    def _valid_reachIDs(self):
        incols = 'reachID' in self.reach_data.columns
        sorted_reaches = sorted(self.reach_data.reachID.tolist())
        consecutive = np.diff(sorted_reaches).sum() \
                      == len(self.reach_data) - 1
        onebased = sorted_reaches.min() == 1
        return incols & consecutive & onebased

    def to_flopy(self):
        m = fm.Modflow(model_ws='',
                       structured=self.structured)

        nstrm = -len(self.reach_data)
        sd = self.segment_data.groupby('per')
        rd = fm.ModflowSfr2.get_empty_reach_data(abs(nstrm))
        sd = {per: sd.get_group(per)
              for per in self.segment_data.per.unique()}
        nss = len(sd[0])
        sfr = fm.ModflowSfr2(nstrm=nstrm, nss=nss)

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
    def from_package(sfrpackagefile, grid):
        pass


