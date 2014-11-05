# program to plot SFR segment profiles

import os
import numpy as np
import discomb_utilities as disutil
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import shutil


import fiona
from shapely.geometry import Polygon, mapping
import pandas as pd




class plot_elevation_profiles:
    # takes information from classes in SFR_classes.py and formats for plotting
    def __init__(self, SFRdata):

        self.SFRdata = SFRdata

        self.elevs_by_cellnum = dict()

        # read in the DIS file for structured grids
        if self.SFRdata.gridtype == 'structured':
            self.read_DIS()


    def read_DIS(self):
        DX, DY, NLAY, self.NROW, self.NCOL, i = disutil.read_meta_data(self.SFRdata.MFdis)

        # get layer tops/bottoms
        self.layer_elevs = np.zeros((NLAY+1, self.NROW, self.NCOL))
        for c in range(NLAY + 1):
            tmp, i = disutil.read_nrow_ncol_vals(self.SFRdata.MFdis, self.NROW, self.NCOL, 'float', i)
            self.layer_elevs[c, :, :] = tmp

        # make dictionary of model top elevations by cellnum
        for c in range(self.NCOL):
            for r in range(self.NROW):
                cellnum = r*self.NCOL + c + 1
                self.elevs_by_cellnum[cellnum] = self.layer_elevs[0, r, c]


    # function to reshape distances and elevations to plot actual cell elevations (flat tops)
    def cells2vertices(self, elevations, distance=False):
        if distance:
            reshaped = [0]
        else:
            reshaped = []
        for i in range(len(elevations)):
            reshaped.append(elevations[i])
            reshaped.append(elevations[i])
        if distance:
            reshaped = reshaped[:-1] # trim last, since last seg distances denotes end of last reach
        return reshaped


    def plot_Mat1_profiles(self, outpdf='Segment_profiles.pdf', interval=1, units='ft', add_profiles={}):
        '''
        Plot elevation profiles from Mat1
        default is to plot streambed top with model top;
        additional profiles can be added with 'add_profiles' argument (they must be added to Mat1 first as additional columns)

        outpdf: filename for multipage PDF
        interval: plot every <interval> segments
        units: (string) for plot lables
        add_profiles: (dict) {<profile name>: <Mat 1 column name>}
        '''
        print 'plotting SFR segment profiles from {}...'.format(self.SFRdata.MAT1)
        try:
            import pandas as pd
        except:
            raise Exception("Requires pandas module, which doesn't seem to be installed")

        m1 = pd.read_csv(self.SFRdata.MAT1)

        # assign node numbers if they aren't in Mat1
        try:
            m1[self.SFRdata.node_attribute]
        except:
            m1[self.SFRdata.node_attribute] = (self.NCOL * (m1['row'] - 1) + m1['column']).astype('int')

        # add model top elevations to dictionary using node numbers
        m1['model_top'] = [self.elevs_by_cellnum[c] for c in m1[self.SFRdata.node_attribute]]

        segs2plot = sorted(np.unique(m1['segment']))[::interval]

        pdf = PdfPages(outpdf)
        for seg in segs2plot:
            print '\r{}'.format(seg),
            df = m1[m1['segment'] == seg].sort('reach')
            df['distance'] = np.cumsum(df['length_in_cell']) - 0.5 * df['length_in_cell'].values
            '''
            dist = self.cells2vertices(df['distance'].tolist(), distance=True)
            sbtops = self.cells2vertices(df['top_streambed'].tolist())
            modeltops = self.cells2vertices(df['model_top'].tolist())
            '''
            dist = df['distance']
            sbtops = df['top_streambed']
            modeltops = df['model_top']

            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.plot(dist, modeltops, label='Model top', lw=0.5)
            plt.plot(dist, sbtops, label='Streambed top', lw=0.5)

            for p in add_profiles.iterkeys():
                #elevs = self.cells2vertices(df[add_profiles[p]].tolist())
                elevs = df[add_profiles[p]]
                plt.plot(dist, elevs, label=p, lw=0.5)

            ax.set_xlabel('Distance along SFR segment')
            ax.set_ylabel('Elevation, {}'.format(units))
            ax.set_title('SFR segment {}'.format(seg))
            ax.set_ylim(int(np.floor(ax.get_ylim()[0])), int(np.ceil(ax.get_ylim()[1])))
            plt.legend()
            pdf.savefig()
            plt.close()
        pdf.close()
        print 'Done, saved to {}'.format(outpdf)



    def get_comid_plotting_info(self, FragIDdata, COMIDdata, SFRdata, interval=False):

        if not interval:
            interval = SFRdata.profile_plot_interval
        self.segs2plot = sorted(FragIDdata.COMID_orderedFragID.keys())[::interval]
        self.seg_dist_dict = dict()
        self.L1top_elev_dict = dict()
        self.seg_elev_fromNHD_dict = dict()
        self.seg_elev_fromContours_dict = dict()
        self.seg_elev_fromDEM_dict = dict()
        self.profiles = [self.L1top_elev_dict, self.seg_elev_fromNHD_dict]
        self.profile_names = ['model top', 'NHDPlus']

        if SFRdata.calculated_contour_elevs:
            self.profiles.append(self.seg_elev_fromContours_dict)
            self.profile_names.append('topographic contours')
        if SFRdata.calculated_DEM_elevs:
            self.profiles.append(self.seg_elev_fromDEM_dict)
            self.profile_names.append('DEM')

        for seg in self.segs2plot:
            distances = []
            elevs_fromNHD = []
            elevs_fromContours = []
            elevs_fromDEM = []
            L1top_top_elevs = []
            dist = 0
            for fid in FragIDdata.COMID_orderedFragID[seg]:
                dist += FragIDdata.allFragIDs[fid].lengthft
                distances.append(dist)
                mean_elev_fromNHD = 0.5 * (FragIDdata.allFragIDs[fid].NHDPlus_elev_max +
                                               FragIDdata.allFragIDs[fid].NHDPlus_elev_min)
                if SFRdata.calculated_contour_elevs:
                    mean_elev_fromContours = 0.5 * (FragIDdata.allFragIDs[fid].interpolated_contour_elev_max +
                                                FragIDdata.allFragIDs[fid].interpolated_contour_elev_min)
                if SFRdata.calculated_DEM_elevs:
                    mean_elev_fromDEM = 0.5 * (FragIDdata.allFragIDs[fid].smoothed_DEM_elev_max +
                                           FragIDdata.allFragIDs[fid].smoothed_DEM_elev_min)
                elevs_fromNHD.append(mean_elev_fromNHD)
                if SFRdata.calculated_contour_elevs:
                    elevs_fromContours.append(mean_elev_fromContours)
                if SFRdata.calculated_DEM_elevs:
                    elevs_fromDEM.append(mean_elev_fromDEM)
                cellnum = FragIDdata.allFragIDs[fid].cellnum
                L1top_top_elevs.append(self.elevs_by_cellnum[cellnum])

            self.seg_dist_dict[seg] = distances
            self.seg_elev_fromNHD_dict[seg] = elevs_fromNHD
            if SFRdata.calculated_contour_elevs:
                self.seg_elev_fromContours_dict[seg] = elevs_fromContours
            if SFRdata.calculated_DEM_elevs:
                self.seg_elev_fromDEM_dict[seg] = elevs_fromDEM
            self.L1top_elev_dict[seg] = L1top_top_elevs


    def get_segment_plotting_info(self, SFRSegsAll, interval=False):

        if not interval:
            interval = self.SFRdata.profile_plot_interval
        seglist = sorted(list(SFRSegsAll.allSegs.keys()))
        self.segs2plot = seglist[::self.SFRdata.profile_plot_interval]
        self.seg_dist_dict = dict()
        self.L1top_elev_dict = dict()
        self.seg_elevs_dict = dict()
        self.profiles = [self.L1top_elev_dict, self.seg_elevs_dict]
        self.profile_names = ['model top', 'streambed top']

        for cseg in seglist:
            reachlist = sorted(SFRSegsAll.allSegs[cseg].seg_reaches)
            curr_reaches = SFRSegsAll.allSegs[cseg].seg_reaches

            distances = []
            L1top_top_elevs = []
            elevs = []
            dist = 0
            for creach in reachlist:
                dist += curr_reaches[creach].eff_length
                distances.append(dist)
                elev = curr_reaches[creach].elevreach
                r, c = curr_reaches[creach].row, curr_reaches[creach].column
                cellnum = (r - 1) * self.NCOL + c
                L1top_top_elevs.append(self.elevs_by_cellnum[cellnum])
                elevs.append(elev)

            self.seg_dist_dict[cseg] = distances
            self.L1top_elev_dict[cseg] = L1top_top_elevs
            self.seg_elevs_dict[cseg] = elevs

    def check4elevation_issues(self, FragIDdata, COMIDdata, SFRdata, SFRSegsAll):
        print "\n\nChecking for 'floating' streambed elevations..."

        # function to build dictionaries of elevations for each segment/COMID and write to output
        def output_elevation_comparison(SFRdata, outfile, header):
            print "Checking by {}".format(header.split(",")[0])
            ofp = open(os.path.join(SFRdata.working_dir, outfile), 'w')
            ofp.write(header)
            for seg in self.segs2plot:
                tops = self.L1top_elev_dict[seg]
                streambed_elevs = self.seg_elevs_dict[seg]
                for i in range(len(tops)):
                    ofp.write('{0},{1},{2},{3}\n'.format(seg, i+1, tops[i], streambed_elevs[i]))

        # check elevations against model top by COMID
        outfile = 'streambed_model_top_comparison.txt'
        header = 'COMID,fragment,streambedtop,modeltop\n'
        self.get_comid_plotting_info(FragIDdata, COMIDdata, SFRdata, interval=1)
        if SFRdata.calculated_contour_elevs:
            self.seg_elevs_dict = self.seg_elev_fromContours_dict
        if SFRdata.calculated_DEM_elevs:
            self.seg_elevs_dict = self.seg_elev_fromDEM_dict
        else:
            self.seg_elevs_dict = self.seg_elev_fromNHD_dict
        output_elevation_comparison(SFRdata, outfile, header)

        # check elevations against model top by segment
        outfile = 'streambed_model_top_comparison.txt'
        header = 'COMID,fragment,streambedtop,modeltop\n'
        self.get_segment_plotting_info(SFRSegsAll, interval=1)
        output_elevation_comparison(SFRdata, outfile, header)


    def plot_profiles(self, pdffile, **kwargs):

        # segs2plot= list of segments to plot
        # seg_distdict= list of distances along segments
        # profiles= list of dictionaries containing profiles for each segment
        # profilenames= list of names, one for each type of profile e.g., model top, STOP post-fix_w_DEM, etc.
        # pdffile= name for output pdf

        try:
            Bottomsdict = kwargs['Bottoms']
            Bottoms = True
        except KeyError:
            Bottoms = False
        try:
            plot_slopes = kwargs['plot_slopes']
            Slopesdict = kwargs['slopes']
            Reach_lengthsdict = kwargs['reach_lengths']
        except KeyError:
            plot_slopes = False

        # function to reshape distances and elevations to plot actual cell elevations
        def reshape_seglist(seg_dict, distance):
            if distance:
                seg_list = [0]
            else:
                seg_list = []
            for i in range(len(seg_dict)):
                seg_list.append(seg_dict[i])
                seg_list.append(seg_dict[i])
            if distance:
                seg_list = seg_list[:-1] # trim last, since last seg distances denotes end of last reach
            return seg_list

        # set name for plotting COMIDs (NHD) or segments (post-NHD)
        try:
            self.seg_elevs_dict
            streamunit = "segment"
        except:
            streamunit = "COMID"

        pdf = PdfPages(pdffile)
        print "\nsaving plots of selected {0}s to {1}".format(streamunit, pdffile)
        knt = 0
        for seg in self.segs2plot:
            knt += 1
            print "\r{0}: {1} ({2} of {3})".format(streamunit, seg, knt, len(self.segs2plot)),
            # reshape distances and elevations to plot actual cell elevations
            seg_distances = reshape_seglist(self.seg_dist_dict[seg], True)
            profiles2plot = []
            for i in range(len(self.profiles)):
                profile = reshape_seglist(self.profiles[i][seg], False)
                profiles2plot.append(profile)

            if Bottoms:
                seg_Bots = reshape_seglist(Bottomsdict[seg], False)
                profiles2plot.append(seg_Bots)
                self.profile_names.append('model bottom')

            if plot_slopes:
                slopes = reshape_seglist(Slopesdict[seg], False)
                reachlengths = reshape_seglist(Reach_lengthsdict[seg], False)


            fig = plt.figure()
            if plot_slopes:
                ((ax1, ax2)) = fig.add_subplot(2, 1, sharex=True, sharey=False)
            else:
                ax1 = fig.add_subplot(1, 1, 1)
            ax1.grid(True)
            colors = ['b', 'g', 'r', 'k']

            for i in range(len(profiles2plot)):
                ax1.plot(seg_distances, profiles2plot[i], color=colors[i], label=self.profile_names[i])

            handles, labels = ax1.get_legend_handles_labels()
            ax1.legend(handles, labels, loc='best')
            ax1.set_title('Streambed profile for {0} {1}'.format(streamunit, seg))
            plt.xlabel('distance along {0} (ft.)'.format(streamunit))
            ax1.set_ylabel('Elevation (ft)')

            # adjust limits to make all profiles visible
            ymax, ymin = np.max(profiles2plot), np.min(profiles2plot)
            ax1.set_ylim([ymin-10, ymax+10])

            # plot segment slopes if desired
            if plot_slopes:
                ax2.grid(True)
                ax2.plot(seg_distances, slopes, color='0.75', label='streambed slopes')
                ax2.set_ylabel('Streambed slope')
                ax3 = ax2.twinx()
                ax3.plot(self.seg_dist_dict[seg], reachlengths, 'b', label='reach length')
                ax3.set_ylabel('reach length (ft)')
                handles, labels = ax2.get_legend_handles_labels()
                ax2.legend(handles, labels, fontsize=6)
                ax3.legend(loc=0)
            pdf.savefig(fig)
            plt.close(fig)
        pdf.close()
        plt.close('all')


class plot_streamflows:
    # plots simulated flows over the SFR network
    def __init__(self, DISfile, streams_shp, SFR_out, node_num_attribute):
        self.streams_shp = streams_shp
        self.SFR_out = SFR_out
        self.flow_by_cellnum = dict()
        self.seg_rch_by_cellnum = dict()
        self.loss_by_cellnum = dict()
        self.overland_by_cellnum = dict()
        self.state_by_cellnum = dict()
        self.stage_by_cellnum = dict()
        self.depth_by_cellnum = dict()
        self.DISfile = DISfile
        self.outpath = os.path.split(SFR_out)[0]
        self.node_num_attribute = node_num_attribute
        if len(self.outpath) == 0:
            self.outpath = os.getcwd()

    def join_SFR_out2streams(self, use_arcpy=True):

        # get model info
        try:
            DX, DY, NLAY, NROW, NCOL, i = disutil.read_meta_data(self.DISfile)
        except:
            raise IOError("Cannot read MODFLOW DIS file {0}".format(self.DISfile))

        print "\naggregating flow information by cellnum..."
        indata = np.genfromtxt(self.SFR_out, skiprows=8, dtype=None)
        for line in indata:
            r, c = line[1], line[2]
            cellnum = (r-1)*NCOL + c
            seg_rch = "{0} {1}; ".format(line[3], line[4])
            flow = 0.5 * (line[5] + line[7])
            loss = float(line[6])
            overland = float(line[8])
            stage = float(line[11])
            depth = float(line[12])

            try:
                existingflow = self.flow_by_cellnum[cellnum]
                seg_rch_info = self.seg_rch_by_cellnum[cellnum]
            except KeyError:
                existingflow = 0
                seg_rch_info = 'segs  rchs: '

            # determine state
            if flow == 0:
                state = 'dry'
            elif loss > 0:
                state = 'loosing'
            elif loss < 0:
                state = 'gaining'
            else:
                print 'Stream reach in cell {} has flow, but no interaction with aquifer.'.format(cellnum)

            self.flow_by_cellnum[cellnum] = existingflow + flow
            self.seg_rch_by_cellnum[cellnum] = seg_rch_info + seg_rch
            self.loss_by_cellnum[cellnum] = loss
            self.state_by_cellnum[cellnum] = state
            self.overland_by_cellnum[cellnum] = overland
            self.stage_by_cellnum[cellnum] = stage
            self.depth_by_cellnum[cellnum] = depth

        # write to temporary output file
        ofp = open(os.path.join(self.outpath, 'temp.csv'), 'w')
        ofp.write('{},row,column,seg_reach,flow,loss,overland,state,stage,depth\n'.format(self.node_num_attribute))
        for cn in self.flow_by_cellnum.keys():
            ofp.write('{0},{1},{2},"{3}",{4:.6e},{5},{6},{7},{8},{9}\n'.format(cn, 1, 1,
                                                                   self.seg_rch_by_cellnum[cn],
                                                                   self.flow_by_cellnum[cn],
                                                                   self.loss_by_cellnum[cn],
                                                                   self.overland_by_cellnum[cn],
                                                                   self.state_by_cellnum[cn],
                                                                   self.stage_by_cellnum[cn],
                                                                   self.depth_by_cellnum[cn]))
        ofp.close()

        outfile = os.path.join(self.outpath, "{0}.shp".format(self.SFR_out[:-4]))
        if use_arcpy:
            try:
                import arcpy
                import SFR_arcpy
            except:
                print 'module arcpy not found!'

            # make feature/table layers
            arcpy.env.workspace = self.outpath
            arcpy.env.overwriteOutput = True
            arcpy.CopyFeatures_management(self.streams_shp, self.streams_shp[:-4]+'_backup.shp')
            arcpy.MakeFeatureLayer_management(self.streams_shp[:-4]+'_backup.shp', "streams")
            arcpy.CopyRows_management(os.path.join(self.outpath, 'temp.csv'), os.path.join(self.outpath, 'temp.dbf'))


            # drop all fields except for cellnum from stream linework
            Fields = arcpy.ListFields("streams")
            Fields = [f.name for f in Fields if f.name not in ["FID", "Shape", self.node_num_attribute]]
            if len(Fields) > 0:
                arcpy.DeleteField_management("streams", Fields)

            SFR_arcpy.general_join(outfile, "streams", self.node_num_attribute, "temp.dbf", self.node_num_attribute, keep_common=True)

        else:
            import sys
            sys.path.append('../../GIS_utils')
            try:
                import GISops
            except:
                print 'GIS_utils.GISops not found!'
            GISops.join_csv2shp(self.streams_shp, self.node_num_attribute, os.path.join(self.outpath, 'temp.csv'), self.node_num_attribute, outfile, how='inner')


class SFRshapefile:
    '''
    Allow for a shapefile to be constructed independent of other shapefiles or main SFR program
    (for usg, will need to bring in another file with node geometric information)
    '''

    def __init__(self, SFRdata, xll=0, yll=0, mult=1, outshp=None, prj=None):

        self.SFRdata = SFRdata
        self.xll = xll
        self.yll = yll
        self.mult = mult # convert model units to GIS units

        if not outshp:
            self.outshp = self.SFRdata.GISSHP
        else:
            self.outshp = outshp

        self.elevs_by_cellnum = {}
        self.prj = prj

        # read in Mat 1 and 2
        self.m1 = pd.read_csv(self.SFRdata.MAT1).sort(['segment', 'reach'])
        self.m2 = pd.read_csv(self.SFRdata.MAT2)
        self.m2.index = self.m2.segment
        self.upsegs = pd.Series(self.m2.segment, index=self.m2.outseg).to_dict()
        self.nSFRcells = len(self.m1)

        # read information for structured grids
        if self.SFRdata.gridtype == 'structured':
            self.read_DIS()

            # add cell geometries to Mat1
            self.m1['geometry'] = [self.cell_geometry_from_rc(self.m1.ix[i, 'row'], self.m1.ix[i, 'column'])
                                   for i in self.m1.index]

            # add node numbers to Mat1
            self.m1[self.SFRdata.node_attribute] = (self.NCOL * (self.m1['row'] - 1) + self.m1['column']).astype('int')

        # add model tops to Mat1
        self.m1['model_top'] = [self.elevs_by_cellnum[c] for c in self.m1[self.SFRdata.node_attribute]]


    def read_DIS(self):
        self.DX, self.DY, self.NLAY, self.NROW, self.NCOL, i = disutil.read_meta_data(self.SFRdata.MFdis)

        self.DX *= self.mult
        self.DY = (self.DY * self.mult)[::-1]
        self.delx = np.append(np.diff(self.DX), np.diff(self.DX)[-1]) # add another spacing on because the vector doesn't include distal edge of grid
        self.dely = np.append(np.diff(self.DY), np.diff(self.DY)[-1])

        # get layer tops/bottoms
        self.layer_elevs = np.zeros((self.NLAY+1, self.NROW, self.NCOL))
        for c in range(self.NLAY + 1):
            tmp, i = disutil.read_nrow_ncol_vals(self.SFRdata.MFdis, self.NROW, self.NCOL, 'float', i)
            self.layer_elevs[c, :, :] = tmp

        # make dictionary of model top elevations by cellnum
        for c in range(self.NCOL):
            for r in range(self.NROW):
                cellnum = r*self.NCOL + c + 1
                self.elevs_by_cellnum[cellnum] = self.layer_elevs[0, r, c]


    def cell_geometry_from_rc(self, r, c):

        # calculate vertices for parent cell
        x0 = self.DX[c-1] + self.xll
        x1 = x0 + self.delx[c-1]
        y0 = self.DY[r-1] + self.yll
        y1 = y0 - self.dely[r-1]

        return Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)])


    def build(self):
        '''
        Build the shapefile
        '''

        schema = {'geometry': 'Polygon',
                      'properties': {self.SFRdata.node_attribute: 'int',
                                     'row': 'int',
                                     'column': 'int',
                                     'layer': 'int',
                                     'segment': 'int',
                                     'reach': 'int',
                                     'outseg': 'int',
                                     'upseg': 'int',
                                     'outlet': 'int',
                                     'sb_elev': 'float',
                                     'modeltop': 'float',
                                     'sb_slope': 'float',
                                     'width': 'float'}}

        with fiona.collection(self.outshp, "w", "ESRI Shapefile", schema) as output:

            knt = 0
            for i in range(self.nSFRcells):

                print "\r{:d}%".format(100 * knt / self.nSFRcells),
                knt += 1

                segment = int(self.m1.ix[i, 'segment'])
                outseg = int(self.m2.ix[segment, 'outseg'])

                # handle headwaters (no upseg)
                try:
                    self.upsegs[segment]
                except KeyError:
                    self.upsegs[segment] = 0

                # shapefiles are incompatible with int64.
                # but apparently they are compatible with float64 ARGH!

                # python int() worked on Mac, whereas .astype('int32') or np.int32() failed
                output.write({'properties': {self.SFRdata.node_attribute:
                                                 int(self.m1.ix[i, self.SFRdata.node_attribute]),
                                             'row': int(self.m1.ix[i, 'row']),
                                             'column': int(self.m1.ix[i, 'column']),
                                             'layer': int(self.m1.ix[i, 'layer']),
                                             'segment': segment,
                                             'reach': int(self.m1.ix[i, 'reach']),
                                             'outseg': outseg,
                                             'upseg': int(self.upsegs[segment]),
                                             'outlet': int(self.m1.ix[i, 'Outlet']),
                                             'sb_elev': self.m1.ix[i, 'top_streambed'],
                                             'modeltop': self.m1.ix[i, 'model_top'],
                                             'sb_slope': self.m1.ix[i, 'bed_slope'],
                                             'width': self.m1.ix[i, 'width_in_cell'],
                },
                              'geometry': mapping(self.m1.ix[i, 'geometry'])})

        if self.prj:
            shutil.copyfile(self.prj, self.outshp[:-4] + '.prj')