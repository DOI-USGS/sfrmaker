__author__ = 'Fienen, Reeves, Leaf - USGS'

import xml.etree.ElementTree as ET
import arcpy
import os
import numpy as np
import SFR_arcpy
import time
import datetime
import discomb_utilities as disutil #  utility to read in a dis file
import pickle
import math
import gzip
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import shutil

'''
debugging functions
'''
def savetmp(savedict):
    for carg in savedict:
        outfilename = '{0:s}###.pklz'.format(carg)
        print "Pickling down {0:s} to file {1:}".format(carg, outfilename)
        ofp = gzip.open(outfilename, 'wb')
        pickle.dump(savedict[carg], ofp)
        ofp.close()

def loadtmp(savedict):
    a = dict()
    for carg in savedict:
        infilename = '{0:s}###.pklz'.format(carg)
        print "loading {0:s} from file {1:}".format(carg, infilename)
        ofp = gzip.open(infilename, 'rb')
        a[carg] = (pickle.load(ofp))
        ofp.close()
    return a

'''
end debugging functions
'''

class SFRInput:
    """
    the SFRInput class holds all data from the XML-based input file
    """
    def __init__(self, infile):
        try:
            inpardat = ET.parse(infile)
        except:
            raise(InputFileMissing(infile))

        inpars = inpardat.getroot()

        try:
            self.working_dir = inpars.findall('.//working_dir')[0].text
        except IndexError:
            self.working_dir = os.getcwd()

        if not os.path.exists(self.working_dir):
            os.makedirs(self.working_dir)
        # note: compute_zonal is not needed, except for by the DEM_elevs_by_cellnum method,
        # which has been more or less replaced by DEM_elevs_by_FragID
        #self.compute_zonal = self.tf2flag(inpars.findall('.//compute_zonal')[0].text)
        self.preproc = self.tf2flag(inpars.findall('.//preproc')[0].text)
        if not self.preproc:
            print "Running without pre-processing (preproc set to False)\n"
        self.reach_cutoff = float(inpars.findall('.//reach_cutoff')[0].text)
        self.rfact = float(inpars.findall('.//rfact')[0].text)
        self.Lowerbot = self.tf2flag(inpars.findall('.//Lowerbot')[0].text)
        self.buff = float(inpars.findall('.//buff')[0].text)
        self.minimum_slope = float(inpars.findall('.//minimum_slope')[0].text)
        self.tpl = self.tf2flag(inpars.findall('.//tpl')[0].text)
        self.MFgrid = inpars.findall('.//MFgrid')[0].text
        self.MFdomain = inpars.findall('.//MFdomain')[0].text
        self.MFdis = inpars.findall('.//MFdis')[0].text
        self.DEM = inpars.findall('.//DEM')[0].text
        self.intersect = self.add_path(inpars.findall('.//intersect')[0].text)
        self.intersect_points = self.add_path(inpars.findall('.//intersect_points')[0].text)
        self.rivers_table = inpars.findall('.//rivers_table')[0].text
        self.PlusflowVAA = inpars.findall('.//PlusflowVAA')[0].text
        self.Elevslope = inpars.findall('.//Elevslope')[0].text
        self.Flowlines_unclipped = inpars.findall('.//Flowlines_unclipped')[0].text
        self.arcpy_path = inpars.findall('.//arcpy_path')[0].text
        self.FLOW = inpars.findall('.//FLOW')[0].text
        self.FTab = inpars.findall('.//FTab')[0].text
        self.Flowlines = self.add_path(inpars.findall('.//Flowlines')[0].text)
        self.ELEV = self.add_path(inpars.findall('.//ELEV')[0].text)
        self.CELLS = self.add_path(inpars.findall('.//CELLS')[0].text)
        self.CELLS_DISS = self.add_path(inpars.findall('.//CELLS_DISS')[0].text)
        self.NHD = self.add_path(inpars.findall('.//NHD')[0].text)
        self.OUT = self.add_path(inpars.findall('.//OUT')[0].text)
        self.MAT1 = self.add_path(inpars.findall('.//MAT1')[0].text)
        self.MAT2 = self.add_path(inpars.findall('.//MAT2')[0].text)
        self.WIDTH = self.add_path(inpars.findall('.//WIDTH')[0].text)
        self.MULT = self.add_path(inpars.findall('.//MULT')[0].text)
        self.ELEVcontours = inpars.findall('.//ELEVcontours')[0].text
        self.Routes = self.add_path(inpars.findall('.//Routes')[0].text)
        self.Contours_intersect = self.add_path(inpars.findall('.//Contours_intersect')[0].text)
        self.Contours_intersect_distances = self.add_path(inpars.findall('.//Contours_intersect_distances')[0].text)
        self.RCH = self.add_path(inpars.findall('.//RCH')[0].text)
        self.nsfrpar = int(inpars.findall('.//nsfrpar')[0].text)
        self.nparseg = int(inpars.findall('.//nparseg')[0].text)
        self.const = float(inpars.findall('.//const')[0].text)
        self.dleak = float(inpars.findall('.//dleak')[0].text)
        self.nstrail = int(inpars.findall('.//nstrail')[0].text)
        self.isuzn = int(inpars.findall('.//isuzn')[0].text)
        self.nsfrsets = int(inpars.findall('.//nsfrsets')[0].text)
        self.istcb1 = int(inpars.findall('.//istcb1')[0].text)
        self.istcb2 = int(inpars.findall('.//istcb2')[0].text)
        self.isfropt = int(inpars.findall('.//isfropt')[0].text)
        self.bedK = float(inpars.findall('.//bedK')[0].text)
        self.bedKmin = float(inpars.findall('.//bedKmin')[0].text)
        self.bedthick = float(inpars.findall('.//bedthick')[0].text)
        self.icalc = int(inpars.findall('.//icalc')[0].text)
        self.nstrpts = int(inpars.findall('.//nstrpts')[0].text)
        self.iprior = int(inpars.findall('.//iprior')[0].text)
        self.flow = float(inpars.findall('.//flow')[0].text)
        self.runoff = float(inpars.findall('.//runoff')[0].text)
        self.etsw = float(inpars.findall('.//etsw')[0].text)
        self.pptsw = float(inpars.findall('.//pptsw')[0].text)
        self.roughch = float(inpars.findall('.//roughch')[0].text)
        self.roughbk = float(inpars.findall('.//roughbk')[0].text)
        self.cdepth = float(inpars.findall('.//cdepth')[0].text)
        self.fdepth = float(inpars.findall('.//fdepth')[0].text)
        self.awdth = float(inpars.findall('.//awdth')[0].text)
        self.bwdth = float(inpars.findall('.//bwdth')[0].text)
        self.thickm1 = float(inpars.findall('.//thickm1')[0].text)
        self.thickm2 = float(inpars.findall('.//thickm2')[0].text)
        self.Hc1fact = float(inpars.findall('.//Hc1fact')[0].text)
        self.Hc2fact = float(inpars.findall('.//Hc2fact')[0].text)
        self.stream_depth = float(inpars.findall('.//stream_depth')[0].text)
        self.GISSHP = self.add_path(inpars.findall('.//GISSHP')[0].text)
        self.elevflag = inpars.findall('.//elevflag')[0].text


        self.calculated_DEM_elevs = False
        self.calculated_contour_elevs = False

        # read in model information
        self.DX, self.DY, self.NLAY, self.NROW, self.NCOL, i = disutil.read_meta_data(self.MFdis)

        # make backup copy of MAT 1, if it exists
        MAT1backup = "{0}_backup".format(self.MAT1)
        if os.path.isfile(self.MAT1):
            try:
                shutil.copyfile(self.MAT1, MAT1backup)
            except IOError:
                print "Tried to make a backup copy of {0} but got an error.".format(self.MAT1)
        try:
            self.plotflag = self.tf2flag(inpars.findall('.//plotflag')[0].text)
        except:
            self.plotflag = True
        try:
            self.eps = float(inpars.findall('.//eps')[0].text)
        except:
            self.eps = 1.0000001e-02  # default value used if not in the input file
        try:
            self.min_elev = float(inpars.findall('.//minimum_elevation')[0].text)
        except:
            self.min_elev = -999999.0  # default value used if not in the input file

        # conversion for vertical length units between NHDPlus v2. (usually cm) and model
        try:
            self.z_conversion = float(inpars.findall('.//z_conversion')[0].text)
        except:
            self.z_conversion = 1.0/(2.54 *12)  # default value used if not in the input file

        # conversion for vertical length units between DEM and model
        try:
            self.DEM_z_conversion = float(inpars.findall('.//DEM_z_conversion')[0].text)
        except:
            self.DEM_z_conversion = 1.0  # default value used if not in the input file

        #cutoff to check stream length in cell against fraction of cell dimension
        #if the stream length is less than cutoff*side length, the piece of stream is dropped
        try:
            self.cutoff = float(inpars.findall('.//cutoff')[0].text)
        except:
            self.cutoff = 0.0

        try:
            self.distanceTol = float(inpars.findall('.//distanceTol')[0].text)
        except:
            self.distanceTol = 3*np.min(self.DX)

        try:
            self.profile_plot_interval = int(inpars.findall('.//profile_plot_interval')[0].text)
        except:
            self.profile_plot_interval = 20

        try:
            self.elev_comp = self.tf2flag(inpars.findall('.//compare_elevation_methods')[0].text)

        except:
            self.elev_comp = False

        # PFJ: The exception assignment to a string when no field is present in the input XML file is problematic because
        # it causes the code to skip the checking analysis near line 1265, where 'cellnum' had previously been generated.
        # I've added a boolean assigment here, and added a "self.node_attribute = 'cellnum' " assigment around line 1295.
        # Of course, now I had to add a check for a newly generated cellnum field in the shapefile (following try/except)

        # if no node attribute is specified, will default to "old way," where it looks for a "cellnum" column,
        # tests if it has unique values, otherwise attempts to make one from row and column columns
        try:
            self.node_attribute = inpars.findall('.//node_attribute')[0].text
        except:
            self.node_attribute = False
        # PFJ: now need to add a test to see if 'cellnum' has been added to self.MFgrid from a prior run with preproc
        fields = arcpy.ListFields(self.MFgrid)
        for field in fields:
            if str(field.name) == "cellnum":
                self.node_attribute = 'cellnum'


        #read the Fcode-Fstring table and save it into a dictionary, Fstring
        descrips = arcpy.SearchCursor(self.FTab)
        self.Fstring = dict()
        for description in descrips:
            Fcodevalue = int(description.FCode)
            if not Fcodevalue in self.Fstring:
                self.Fstring[Fcodevalue]=description.Descriptio
        del descrips

        # initialize the arcpy environment
        #arcpy.env.workspace = os.getcwd()
        arcpy.env.workspace = self.working_dir
        arcpy.env.overwriteOutput = True
        arcpy.env.qualifiedFieldNames = False
        # Check out any necessary arcpy licenses
        arcpy.CheckOutExtension("spatial")

    def tf2flag(self, intxt):
        # converts text written in XML file to True or False flag
        if intxt.lower() == 'true':
            return True
        else:
            return False

    def add_path(self, infile):

        # if there is no path in front of the input file, add one
        if len(os.path.split(infile)[0]) == 0:
            fullpath = os.path.join(self.working_dir, infile)

        # otherwise, leave the path as is
        else:
            fullpath = infile
        return fullpath


class FragIDprops(object):
    """
    Properties for each COMID
    """
    '''
    __slots__ = ['comid', 'startx', 'starty', 'endx', 'endy', 'FragID',
                 'maxsmoothelev', 'minsmoothelev', 'lengthft',
                 'cellnum', 'elev', 'sidelength',
                 'start_has_end', 'end_has_start', 'contour_elev', 'contour_elev_distance', 'segelevinfo',
                 'NHDPlus_elev_min', 'NHDPlus_elev_max', 'NHDPlus_elev_mean', 'NHDPlus_slope',
                 'interpolated_contour_elev_min', 'interpolated_contour_elev_max', 'interpolated_contour_elev_mean',
                 'interpolated_contour_slope', 'DEM_elev_min', 'DEM_elev_max', 'DEM_elev_mean', 'DEM_slope',
                 'smoothed_DEM_elev_min', 'smoothed_DEM_elev_max', 'smoothed_DEM_elev_mean', 'smoothed_DEM_slope']
    '''
    #  using __slots__ makes it required to declare properties of the object here in place
    #  and saves significant memory
    def __init__(self, comid, startx, starty, endx, endy, FragID,
                 maxsmoothelev, minsmoothelev, lengthft, cellnum, contour_elev, contour_elev_distance,
                 NHDPlus_elev_min, NHDPlus_elev_max, NHDPlus_elev_mean, NHDPlus_slope,
                 interpolated_contour_elev_min, interpolated_contour_elev_max, interpolated_contour_elev_mean,
                 interpolated_contour_slope, DEM_elev_min, DEM_elev_max, DEM_elev_mean, DEM_slope,
                 smoothed_DEM_elev_min, smoothed_DEM_elev_max, smoothed_DEM_elev_mean, smoothed_DEM_slope):
        self.comid = comid
        self.startx = startx
        self.starty = starty
        self.endx = endx
        self.endy = endy
        self.FragID = FragID
        self.maxsmoothelev = maxsmoothelev
        self.minsmoothelev = minsmoothelev
        self.lengthft = lengthft
        self.cellnum = cellnum
        self.contour_elev = contour_elev
        self.contour_elev_distance = contour_elev_distance
        self.NHDPlus_elev_min = NHDPlus_elev_min
        self.NHDPlus_elev_max = NHDPlus_elev_max
        self.NHDPlus_elev_mean = NHDPlus_elev_mean
        self.NHDPlus_slope = NHDPlus_slope
        self.interpolated_contour_elev_min = interpolated_contour_elev_min
        self.interpolated_contour_elev_max = interpolated_contour_elev_max
        self.interpolated_contour_elev_mean = interpolated_contour_elev_mean
        self.interpolated_contour_slope = interpolated_contour_slope
        self.DEM_elev_min = DEM_elev_min
        self.DEM_elev_max = DEM_elev_max
        self.DEM_elev_mean = DEM_elev_mean
        self.DEM_slope = DEM_slope
        self.smoothed_DEM_elev_min = smoothed_DEM_elev_min
        self.smoothed_DEM_elev_max = smoothed_DEM_elev_max
        self.smoothed_DEM_elev_mean = smoothed_DEM_elev_mean
        self.smoothed_DEM_slope = smoothed_DEM_slope
        #self.slope = slope


class COMIDprops(object):
    """
    routing information by COMIDs
    """
    '''
    __slots__ = ['from_comid', 'to_comid', 'hydrosequence', 'uphydrosequence',
                 'downhydrosequence', 'levelpathID','stream_order','arbolate_sum',
                 'est_width','reachcode','Fcode']
    '''
    def __init__(self):
        self.from_comid = list()
        self.to_comid = list()
        self.hydrosequence = None
        self.uphydrosequence = None
        self.downhydrosequence = None
        self.levelpathID = None
        self.stream_order = None
        self.arbolate_sum = None
        self.est_width = None
        self.reachcode = None
        self.Fcode = None


class LevelPathIDprops(object):
    """
    routing of LevelPathIDs
    """
    '''
    __slots__ = ["down_levelpathID", "ordered_cellnums", "ordered_FragIDs", "ordered_hydrosequence"]
    '''
    def __init__(self):
        self.down_levelpathID = None
        self.ordered_hydrosequence = list()
        self.ordered_cellnums = list()  # NB - this is unique even though duplicates may have existed due to meanders
        self.ordered_FragIDs = list()


class LevelPathIDpropsAll:
    def __init__(self):
        self.allids = dict()
        self.level_ordered = list()
        self.levelpath_FragID = dict()

    def return_cutoffs(self, FragIDdata, CELLdata, SFRdata):
        for lpID in self.level_ordered:
            #check to see if individual reachlengths are less than cutoff
            #prescibed by sidelength*cutoff
            rmlist=[]
            for FragID in self.levelpath_FragID[lpID]:
                reachlength=FragIDdata.allFragIDs[FragID].lengthft
                cellnum = FragIDdata.allFragIDs[FragID].cellnum
                if reachlength < np.min(SFRdata.DX)*SFRdata.cutoff:
                    rmlist.append(FragID)
            #if any were too short remove from levelpath list of FragIDs
            newlist = [FragID for FragID in self.levelpath_FragID[lpID] if FragID not in rmlist]
            self.levelpath_FragID[lpID] = newlist
        rmlist = []
        #if any of the levelpath list of FragIDs is now empty, remove that levelpathID
        for lpID in self.level_ordered:
            if len(self.levelpath_FragID[lpID]) == 0:
                rmlist.append(lpID)
        newlist = [lpID for lpID in self.level_ordered if lpID not in rmlist]
        self.level_ordered = newlist


class CellProps(object):
    """
    class for cell objects
    """
    '''
    __slots__ = ['delx', 'dely', 'centroid','sidelength', 'row', 'column', 'fromCell']
    '''
    def __init__(self, delx, dely, centroid, sidelength, row, column):
        self.delx = delx
        self.dely = dely
        self.centroid = centroid
        self.sidelength = sidelength
        self.row = row
        self.column = column
        self.fromCell = []


class CellPropsAll:
    def __init__(self):
        self.allcells = dict()
        self.centroids = dict()

    def populate_cells(self, SFRdata):
        #use the CELLS shapefile to get the length of the cell sides, used
        #to weed out short river reaches.  If the length in the cell is
        #less than sidelength * the input 'cutoff', also get row and column
        #for each cellnumber


        cells = arcpy.da.SearchCursor(SFRdata.CELLS, ['SHAPE@XY', 'delx', 'dely', SFRdata.node_attribute, 'row', 'column'])

        try:
            cells.fields
        except RuntimeError:
            print "At least one of the attribute names in {} is invalid. If grid is structured, the attributes " \
                  "containing row and column information should be called 'row' and 'column' respectively.".format(SFRdata.MFgrid)
        for cell in cells:
            cellnum = int(cell[3])
            dx = float(cell[1])
            dy = float(cell[2])
            centroid = cell[0]
            minside = min([dx, dy])
            row = int(cell[4])
            column = int(cell[5])
            self.allcells[cellnum] = CellProps(dx, dy, centroid, minside, row, column)


    def return_centroids(self):
        # returns a dictionary, keyed by cellnum, with cell centroids
        for ccell in self.allcells:
            self.centroids[ccell] = self.allcells[ccell].centroid


class COMIDPropsAll:
    def __init__(self):
        self.allcomids = dict()
        self.hydrosequence_sorted = list()
        self.hydrosequence_comids = dict()

    def populate_routing(self, SFRdata, FragIDdata, LevelPathdata, CELLdata):
        """
        Read the COMID routing information from the SFRdata.FLOW file
        """
        print ('Reading in routing information from {0:s}'.format(SFRdata.FLOW))
        # open the SFRdata.FLOW file as read-only (using SearchCursor)

        CLIP = np.loadtxt(os.path.join(SFRdata.working_dir, 'boundaryClipsRouting.txt'), skiprows=1, delimiter=',', dtype=int)

        if len(CLIP) == 0:
            print "Warning! boundaryClipsRouting.txt is empty. This implies that there are no NHD streams crossing" \
                "the boundary of the SFR area specified by ()".format(SFRdata.MFdomain)
            CLIP = np.zeros((1, 2))
        for ccomid in FragIDdata.allcomids:
            self.allcomids[ccomid] = COMIDprops()

        # if there are no boundary clips to route, skip it (IndexError happens otherwise)
        knt = 0
        ncomids = len(FragIDdata.allcomids)
        nfragIDs = len(FragIDdata.allFragIDs)
        print "this may take awhile..."
        with arcpy.da.SearchCursor(SFRdata.FLOW, ("FROMCOMID", "TOCOMID")) as cursor:
            for crow in cursor:
                if int(crow[0]) in FragIDdata.allcomids:
                    if (crow[0]) in CLIP[:, 0]:
                        self.allcomids[crow[0]].to_comid.append(999999)
                    else:
                        self.allcomids[crow[0]].to_comid.append(int(crow[1]))

                if int(crow[1]) in FragIDdata.allcomids:
                    if crow[1] in CLIP[:, 1]:
                        self.allcomids[crow[1]].from_comid.append(999999)
                    else:
                        self.allcomids[crow[1]].from_comid.append(int(crow[0]))

        del crow, cursor

        print ('\nReading in routing information from {0:s}'.format(SFRdata.PlusflowVAA))
        comidseen = list()
        with arcpy.da.SearchCursor(SFRdata.PlusflowVAA,
                                   ("ComID", "Hydroseq", "uphydroseq", "dnhydroseq",
                                   "ReachCode","StreamOrde","ArbolateSu","Fcode", "levelpathI")) as cursor:
            knt = 0
            for crow in cursor:
                comid = int(crow[0])

                if int(comid) in FragIDdata.allcomids:
                    # moved these assignments here so they aren't done unless comid in allcomids
                    hydrosequence = int(crow[1])
                    uphydrosequence = int(crow[2])
                    downhydrosequence = int(crow[3])
                    reachcode = crow[4]
                    stream_order = int(crow[5])
                    arbolate_sum = float(crow[6])
                    Fcode = int(crow[7])
                    levelpathid = int(crow[8])
                    self.allcomids[crow[0]].hydrosequence = hydrosequence
                    self.allcomids[crow[0]].uphydrosequence = uphydrosequence
                    self.allcomids[crow[0]].downhydrosequence = downhydrosequence
                    self.allcomids[crow[0]].levelpathID = levelpathid
                    LevelPathdata.level_ordered.append(levelpathid)
                    self.allcomids[crow[0]].reachcode = reachcode
                    self.allcomids[crow[0]].stream_order = stream_order
                    self.allcomids[crow[0]].arbolate_sum = arbolate_sum
                    self.allcomids[crow[0]].Fcode = Fcode
                    #estimate the width
                    self.allcomids[crow[0]].est_width = widthcorrelation(arbolate_sum)
                    comidseen.append(comid)
                    knt += 1
                    print "\r{:d}%".format(100 * knt / ncomids),

            # find unique levelpathIDs
            LevelPathdata.level_ordered = sorted(list(set(LevelPathdata.level_ordered)), reverse=True)

            for clevelpathid in LevelPathdata.level_ordered:
                LevelPathdata.allids[clevelpathid] = LevelPathIDprops()
                LevelPathdata.levelpath_FragID[clevelpathid] = []


            # assign levelpathID routing
        del crow, cursor

        knt = 0
        with arcpy.da.SearchCursor(SFRdata.PlusflowVAA,
                                   ("ComID", "LevelPathI", "DnLevelPat")) as cursor:
            for crow in cursor:
                comid = int(crow[0])
                levelpathid = int(crow[1])
                downlevelpathid = int(crow[2])
                if levelpathid in LevelPathdata.level_ordered:
                    if downlevelpathid != levelpathid:
                        LevelPathdata.allids[levelpathid].down_levelpathID = downlevelpathid
                    ### KLUDGE for Missori River
                    if levelpathid == 550000017:
                        LevelPathdata.allids[levelpathid].down_levelpathID = downlevelpathid
                    if comid in FragIDdata.comid_FragID:
                        LevelPathdata.levelpath_FragID[levelpathid].extend(FragIDdata.comid_FragID[comid])
                        knt += 1
                    print "\r{:d}%".format(100 * knt / ncomids),

        comid_missing = list(set(FragIDdata.allcomids).difference(comidseen))
        if len(comid_missing) > 0:
            print "WARNING! the following COMIDs are missing from \n{0:s}".format('\n'.join(map(str(comid_missing))))

        #populate the cell-to-cell routing dictionary, fromCell, using the routing information.
        #this dictionary returns a list of cells and is keyed by cellnum- the list are the next downstream
        #cell from the key.  It is a list because multiple stream segments could be going through a cell
        #and connecting to different downstream cells

        for ccomid in FragIDdata.COMID_orderedFragID:
            if not comid in FragIDdata.noelev:
                ordfragid = FragIDdata.COMID_orderedFragID[ccomid]
                for i in range(0,len(ordfragid)-1):
                    cfid = ordfragid[i]
                    cfidp1 = ordfragid[i+1]
                    fmcell = FragIDdata.allFragIDs[cfid].cellnum
                    tocell = FragIDdata.allFragIDs[cfidp1].cellnum
                    if not tocell == fmcell:
                        if fmcell not in CELLdata.allcells:
                            CELLdata.allcells[frmcell].fromCell = tocell
                        else:
                            CELLdata.allcells[fmcell].fromCell.append(tocell)
                lastcell = FragIDdata.allFragIDs[ordfragid[-1]].cellnum
                nextcomid = self.allcomids[ccomid].to_comid[0]
                if nextcomid in FragIDdata.COMID_orderedFragID:
                    nextordfragid = FragIDdata.COMID_orderedFragID[nextcomid]
                    nextfragid = nextordfragid[0]
                    if nextfragid in FragIDdata.allFragIDs:
                        nextcell = FragIDdata.allFragIDs[nextfragid].cellnum
                        if lastcell not in CELLdata.allcells:
                            CELLdata.allcells[lastcell].fromCell = nextcell
                        else:
                            CELLdata.allcells[lastcell].fromCell.append(nextcell)


    def return_hydrosequence_comid(self):
        """
        return a dictionary of hydrosequence linked up with COMIDs
        """
        # first get a unique set of hydrosequences
        for ccomid in self.allcomids:
            self.hydrosequence_sorted.append(self.allcomids[ccomid].hydrosequence)
        self.hydrosequence_sorted.sort(reverse=True)
        for ccomid in self.allcomids:
            self.hydrosequence_comids[self.allcomids[ccomid].hydrosequence] = ccomid


class COMIDPropsForIntersect:
    """
    Properties for each COMID
    """
    def __init__(self, comid, inout, startx, starty, endx, endy, clmaxel,
                 clminel, cllen, lenkm):
        self.comid = comid
        self.inout = inout
        self.newstartx = startx
        self.newstarty = starty
        self.newnewendx = endx
        self.newendy = endy
        self.newlength = cllen
        self.slope = clmaxel - clminel
        if inout == 'OUT':
            self.newmaxel = round(clmaxel)
            self.newminel = round(clmaxel - self.slope * cllen / lenkm)
        elif inout == 'IN':
            clippedlength = lenkm - cllen
            self.newminel = round(clminel)
            self.newmaxel = round(clmaxel - self.slope * clippedlength / lenkm)


class FragIDPropsAll:
    def __init__(self):
        self.allFragIDs = dict()
        self.allcomids = list()  # comprehensive list of all comids
        self.unique_cells = list()  # list of unique cellnum values in the grid/streams intersection
        self.maxelev = None
        self.minelev = None
        self.noelev = dict()
        self.COMID_orderedFragID = dict()
        self.comid_FragID = None
        self.cellnum_FragID = None
        self.cellnum_LevelPathID = None

    def return_FragID_comid_list(self):
        """
        Return a dict of lists of comids linked with each FragID
        """
        self.comid_FragID = {ccomid: [] for ccomid in self.allcomids}
        allFragIDs = self.allFragIDs.keys()
        for cFragID in allFragIDs:
            self.comid_FragID[self.allFragIDs[cFragID].comid].append(cFragID)

    def return_smoothelev_comid(self, comid):
        self.maxelev = -np.inf
        self.minelev = np.inf
        FragIDs = self.comid_FragID[comid]
        for cid in FragIDs:
            if self.allFragIDs[cid].maxsmoothelev > self.maxelev:
                self.maxelev = self.allFragIDs[cid].maxsmoothelev
            if self.allFragIDs[cid].minsmoothelev < self.minelev:
                self.minelev = self.allFragIDs[cid].minsmoothelev

    def return_unique_cells(self):
        for cid in self.allFragIDs.keys():
            self.unique_cells.append(self.allFragIDs[cid].cellnum)
        self.unique_cells = set(self.unique_cells)

    def return_cellnum_FragID(self):
        """
        method to return a dictionary keyed by cellnum returning a list FragIDs - used
        to relate the properties in the FragID object to the cellnumber
        """
        self.cellnum_FragID = {cellnum: [] for cellnum in self.unique_cells}
        for cFragID in self.allFragIDs.iterkeys():
            self.cellnum_FragID[self.allFragIDs[cFragID].cellnum].append(cFragID)

    def return_cellnum_LevelPathID(self, LevelPathdata):
        """
        method to return a list of levelpathIDs within a given cell; keyed by
        cellnum and returning a list of LevelPathIDs - especially used where properties are
        accumulated into reaches for a given cell to distinguish between different segments
        """
        self.cellnum_LevelPathID = {cellnum : [] for cellnum in self.unique_cells}
        #each FragID has only one level path, invert the levelpath/fragID dictionary
        #accounting for the value of the dict is a list of FragIDs
        fragid_lpid=dict()
        for lpID, fragidlist in LevelPathdata.levelpath_FragID.iteritems():
            for frg in fragidlist:
                fragid_lpid[frg]=lpID

        for localcell, fragidlist in self.cellnum_FragID.items():
            for cfragid in fragidlist:
                self.cellnum_LevelPathID[localcell].append(fragid_lpid[cfragid])

        #go over the localcells once more and make the list unique
        for localcell in self.cellnum_LevelPathID.iterkeys():
            self.cellnum_LevelPathID[localcell] = list(set(self.cellnum_LevelPathID[localcell]))




    def populate(self, SFRdata):
        """
        read in the main COMID-related properties from the intersect file
        """
        segments = arcpy.SearchCursor(SFRdata.intersect)

        for seg in segments:
            FragID = int(seg.FragID)
            self.allFragIDs[FragID] = FragIDprops(
                int(seg.COMID),
                float(seg.X_start),
                float(seg.Y_start),
                float(seg.X_end),
                float(seg.Y_end),
                int(seg.FragID),
                float(seg.MAXELEVSMO)*SFRdata.z_conversion,  # UNIT CONVERSION
                float(seg.MINELEVSMO)*SFRdata.z_conversion,  # UNIT CONVERSION
                float(seg.LengthFt),
                seg.cellnum,
                list(), list(), None, None, None, None, None, None, None, None,
                None, None, None, None, None, None, None, None)


            self.allcomids.append(int(seg.COMID))
        self.allcomids = list(set(self.allcomids))
        self.return_unique_cells()
        self.return_cellnum_FragID()




    def populate_elevations(self, SFRdata):
        """
        Read elevation information, per COMID, from the SFRdata.ELEV file
        """


        arcpy.MakeFeatureLayer_management(SFRdata.intersect, "grid_temp")
        SFR_arcpy.general_join(SFRdata.ELEV, "grid_temp", "FragID", SFRdata.rivers_table, "OLDFragID", keep_common=False)

        with arcpy.da.SearchCursor(SFRdata.ELEV, ("FragID", "ELEVAVE")) as cursor:
            for crow in cursor:
                self.allFragIDs[int(crow[0])].elev = float(crow[1])


class SFRReachProps(object):
    """
    class containing an object for each reach
    """
    '''
    __slots__ = ['cellnum','eff_length','eff_width','eff_slope','elevreach','bedthick','bedK']
    '''
    def __init__(self, cellnum, eff_length, eff_width, eff_slope, elevreach, bedthick, bedK):
        self.cellnum = cellnum
        self.row = None
        self.column = None
        self.eff_length = eff_length
        self.eff_width = eff_width
        self.eff_slope = eff_slope
        self.elevreach = elevreach
        self.bedthick = bedthick
        self.bedK = bedK


class SFRSegmentProps(object):
    """
    class object to hold final SFR segment information
    """
    '''
    __slots__ = ['seg_cells','icalc','iupseg','outseg','runoff','etsw','pptsw',
    'seg_reaches','seg_label','levelpathID','inseg', 'roughch', 'startreach_xy', 'endreach_xy']
    '''
    def __init__(self):
        self.seg_cells = list()
        self.outseg = None
        self.endreach_xy = None
        self.startreach_xy = None
        self.icalc = None
        self.iupseg = 0    #iupseg default zero right now, diversions could be added
        self.outseg = None
        self.runoff = None
        self.etsw = None
        self.pptsw = None
        self.roughch = None
        self.seg_reaches = dict()   #reach object with properties from SFRReachProps keyed by rch_number
        self.seg_label = None  #segment label from confluence step
        self.levelpathID = None          #levelpathID corresponding to the segment
        self.inseg = list()     #list of segments connecting to the current segment,
                                 #convenience list for looking at connections


class SFRSegmentsAll:
    """
    class that makes up a dictionary of SFR objects
    and contains methods to accumulate properties of the
    same levelpathID within a call, find confluences, subdivide
    levelpathIDs and establish the final SFR segments
    """
    def __init__(self):
        self.allSegs = dict()              # dictionary of segment objects keyed by final segment number
        self.confluences = defaultdict(list)    # dictionary of lists keyed by levelpathID
        self.segment_levelpaths = dict()
        self.levelpaths_segments = dict()

    def divide_at_confluences(self, LevelPathdata, FragIDdata, COMIDdata, CELLdata, SFRdata):
        #establish provisional segment numbers from downstream (descending)
        #list of levelpathIDs
        provSFRseg = dict()
        for i, lpID in enumerate(LevelPathdata.level_ordered):
            provSFRseg[lpID] = i+1

        CELLdata.return_centroids()
        LP_check_log = open(os.path.join(SFRdata.working_dir, 'LevelPath_check_log.txt'), 'w')
        print 'finding confluences for provisional SFR segments'
        for clevelpathid in LevelPathdata.level_ordered:
            nextlevelpath = LevelPathdata.allids[clevelpathid].down_levelpathID
            cellist = LevelPathdata.allids[clevelpathid].ordered_cellnums
            if len(cellist) == 0:
                print 'check levelpathID {0}, no cells are assigned'.format(clevelpathid)
                LP_check_log.write('check {}, no cells are assigned\n'.format(clevelpathid))
            else:
                if nextlevelpath == 0:
                    outseg = 0
                else:
                    if nextlevelpath in provSFRseg:
                        outseg = provSFRseg[nextlevelpath]
                    else:
                        outseg = int(0)
                isegend = cellist[-1]
                if outseg > 0:
                    outid = nextlevelpath
                    confl = -1    # flag in case an end is not found
                    for cell in LevelPathdata.allids[outid].ordered_cellnums:
                        if cell == isegend:
                            confl = cell
                else:
                    #no downstream levelpath
                    outid = 0
                    confl = 0
                #no end found
                if confl == -1:
                    lastfragid = LevelPathdata.allids[clevelpathid].ordered_FragIDs[-1]
                    lastcomid = FragIDdata.allFragIDs[lastfragid].comid
                    #check down_hydrosequence
                    if COMIDdata.allcomids[lastcomid].downhydrosequence == 0:
                        confl = -2     #its a downstream end within the model
                #put confluence cells into lists, iseg end is a confluence for iseg
                self.confluences[clevelpathid].append(isegend)
                #if a confluence was found, put it on the confluence list
                if confl > 0:
                    self.confluences[outid].append(confl)
                if confl == -1:
                    print 'check downstream connection for levelpathID {0}'.format(clevelpathid)
                    LP_check_log.write('check {}, possible issue with downstream connection\n'.format(clevelpathid))
                #put it in downstream order with no repeats
                self.confluences[clevelpathid]=[cell for cell in LevelPathdata.allids[clevelpathid].ordered_cellnums \
                    if cell in set(self.confluences[clevelpathid])]
        #now use the list of confluences for each levelpathid to subdivide and
        #assign cells to each subsection - these become the final SFR segments
        subseg = dict()
        subconfl = defaultdict(list)
        subcell = dict()
        subordered_cells = dict()
        conf_count = 1
        CHK2 = open('check2.out', 'w')
        CHK2.write("LevelpathID, provisionalSegment, Sublabel, finalseg, cells\n")
        for clevelpathid in LevelPathdata.level_ordered:
            iseg = provSFRseg[clevelpathid]
            numconfls = len(self.confluences[clevelpathid])
            #break up segments with multiple confluences
            strt = 0
            for confl in range(0, numconfls):
                sublabel = str(iseg) + '-' + str(confl)
                subseg[sublabel] = conf_count
                subconfl[iseg].append(sublabel)
                subcell[sublabel] = self.confluences[clevelpathid][confl]
                conf_count += 1
                #build  cell lists for each subsection of a provisional segment
                #find the index that matches the cell given by conflunces[id][confl]
                lpcells = LevelPathdata.allids[clevelpathid].ordered_cellnums
                endindx = lpcells.index(self.confluences[clevelpathid][confl])
                subordered_cells[sublabel] = lpcells[strt:endindx+1]
                CHK2.write('{0}, {1}, {2}, {3}:'.format(clevelpathid, iseg, sublabel, subseg[sublabel]))
                CHK2.write(','.join(map(str, subordered_cells[sublabel]))+'\n')
                strt = endindx + 1
        CHK2.close()
        #use confluence information to build final SFR segments and reaches
        for clevelpathid, iseg in provSFRseg.iteritems():
            for i in range(0, len(subconfl[iseg])):
                labl = subconfl[iseg][i]
                self.allSegs[subseg[labl]] = SFRSegmentProps()
                self.allSegs[subseg[labl]].seg_cells = subordered_cells[labl]
                self.allSegs[subseg[labl]].startreach_xy = CELLdata.centroids[
                    self.allSegs[subseg[labl]].seg_cells[0]
                ]
                self.allSegs[subseg[labl]].endreach_xy = CELLdata.centroids[
                    self.allSegs[subseg[labl]].seg_cells[-1]
                ]
                self.allSegs[subseg[labl]].seg_label = labl
                self.allSegs[subseg[labl]].levelpathID = clevelpathid
                lastcell = subordered_cells[labl][-1]
                nextlabl = 999999  # initialize in case nextlevelpath is not found
                if i < len(subconfl[iseg])-1:
                    nextlevelpath = clevelpathid
                    nextlabl = subconfl[iseg][i+1]
                else:
                    nextlevelpath = LevelPathdata.allids[clevelpathid].down_levelpathID
                    if nextlevelpath in provSFRseg and nextlevelpath > 0:
                        try:
                            nextlabl = subconfl[provSFRseg[nextlevelpath]][0]
                        except:
                            nextlabl = 999999

                if nextlevelpath == 0:
                    outseg = 0
                elif nextlabl == 999999:
                    outseg = 999999
                else:
                    outseg = subseg[nextlabl]

                # now check for circular routing (typically happens when only one reach in a segment)subseg[labl]
                # this repeats much of the logic of the preceding 10 or so lines
                if subseg[labl] == outseg:
                    circular_condition = True
                    circular_test_counter = 0
                    circ_thresh = 10  # number of iterations of circular test loop below
                    while circular_condition and circular_test_counter < circ_thresh:
                        circular_test_counter += 1
                        if nextlevelpath > 0:
                            nextlevelpath = LevelPathdata.allids[nextlevelpath].down_levelpathID
                            if nextlevelpath in provSFRseg and nextlevelpath > 0:
                                nextlabl = subconfl[provSFRseg[nextlevelpath]][0]
                            else:
                                nextlabl = 999999

                        if nextlevelpath == 0:
                            outseg = 0
                        elif nextlabl == 999999:
                            outseg = 999999
                        else:
                            outseg = subseg[nextlabl]
                        if subseg[labl] != outseg:
                            circular_condition = False
                    if circular_condition == True and circular_test_counter >= circ_thresh:
                        outseg = 999999

                self.allSegs[subseg[labl]].outseg = outseg

        # go through the segments and make sure the outseg values are correct, based on being closest to the
        # furthest upstream reach of the assumed outseg. If this is not correct, reassign the outseg

        self.return_segs_levelpaths()
        self.return_levelpaths_segs()

        for seg in self.allSegs.iterkeys():
            outseg = self.allSegs[seg].outseg
            if outseg > 0 and outseg < 999999:
                # find the levelpath ids for current seg and downseg
                levelpathid_outseg = self.segment_levelpaths[outseg]
                x_out = list()
                y_out = list()

                for clpid in self.levelpaths_segments[levelpathid_outseg]:
                    x_out.append(self.allSegs[clpid].startreach_xy[0])
                    y_out.append(self.allSegs[clpid].startreach_xy[1])
                x_out = np.array(x_out)
                y_out = np.array(y_out)
                x, y = self.allSegs[seg].endreach_xy
                dist = np.sqrt((x-x_out)**2 + (y-y_out)**2)
                minloc = np.squeeze(np.where(dist == np.min(dist)))
                if len(np.atleast_1d(minloc)) == 1:
                    xxxx = 5
                if minloc.ndim > 0: # Only apply the following to 1-d or greater arrays
                    minloc = minloc[0] # If 2 segs have same dist, just chose the first one to prevend error...may need to test this
                # if closest segment is the current segment, choose next closest (to prevent circular routing!)
                if self.levelpaths_segments[levelpathid_outseg][minloc] == seg:
                    minloc = np.argsort(dist)[1]
                # if current outseg is already the closest, keep it
                if self.allSegs[seg].outseg != self.levelpaths_segments[levelpathid_outseg][minloc]:
                    print 'Changing outseg from {0:d} to {1:d}'.format(
                        self.allSegs[seg].outseg,
                        self.levelpaths_segments[levelpathid_outseg][minloc])
                    self.allSegs[seg].outseg = self.levelpaths_segments[levelpathid_outseg][minloc]
        print "\n"

        # make a correction to find outsegs equal to 999999 that should really be 0
        # in other words, find cells that are on the model boundary

        # first find which cells border the model boundary
        if arcpy.Exists('boundary_cells.shp'):
            arcpy.Delete_management('boundary_cells.shp')
        arcpy.SpatialJoin_analysis(SFRdata.CELLS_DISS,
                                   SFRdata.MFdomain[:-4] + "_outline.shp",
                                   "boundary_cells",
                                   "JOIN_ONE_TO_MANY",
                                   "KEEP_COMMON",
                                   "BOUNDARY_TOUCHES")
        arcpy.MakeFeatureLayer_management("boundary_cells.shp", "boundary_cells")

        SFR_arcpy.general_join("boundary_cells2.shp",
                               "boundary_cells",
                               "TARGET_FID",
                               SFRdata.CELLS_DISS,
                               "FID")

        # make a list of the cellnums for the border cells
        self.boundary_cellnums = list()
        #with arcpy.da.SearchCursor("boundary_cells2.shp", 'cellnum') as cells:
        with arcpy.da.SearchCursor("boundary_cells2.shp", SFRdata.node_attribute) as cells:
            for cell in cells:
                self.boundary_cellnums.append(int(cell[0]))

        #make a list of what segments are connected to each using outseg information
        #also load in other segment information from SFRdata (XML file information)

        for seg in self.allSegs.iterkeys():
            outseg = self.allSegs[seg].outseg
            if outseg > 0 and outseg < 999999:
                #if outseg is > 0 and not 999999, check that the first reach
                #of the outseg is within a tolerance of the last reach of seg
                #problems occured in NACP model when downstream portions of
                #esturary streams were cut from the model and tributaries were
                #mistakenly routed back up to the remaining headwater.  If the
                #cell that seg is routed to is too far away, set the outseg to zero (boundary).
                npoutpt = np.array(self.allSegs[seg].endreach_xy)
                npinpt = np.array(self.allSegs[outseg].startreach_xy)
                dist = np.sqrt((npoutpt[0]-npinpt[0])**2 + (npoutpt[1]-npinpt[1])**2)
                #tol = CELLdata.allcells[self.allSegs[seg].seg_cells[-1]].sidelength* np.sqrt(2)* 2
                if dist > SFRdata.distanceTol:
                    self.allSegs[seg].outseg = 0
                else:
                    self.allSegs[outseg].inseg.append(seg)
            self.allSegs[seg].etsw = SFRdata.etsw
            self.allSegs[seg].pptsw = SFRdata.pptsw
            self.allSegs[seg].runoff = SFRdata.runoff
            self.allSegs[seg].roughch = SFRdata.roughch
            self.allSegs[seg].icalc = SFRdata.icalc




            # finally shift cases where routing is to 999999 AND cells is on the boundary
            # to make outseg = 0 (outlet)
            if outseg == 999999:
                if self.allSegs[seg].seg_cells[-1] in self.boundary_cellnums:
                    self.allSegs[seg].outseg = 0

    def return_segs_levelpaths(self):
        # return a dictionary, keyed by segment, of levelpathIDs
        for cseg in self.allSegs:
            self.segment_levelpaths[cseg] = self.allSegs[cseg].levelpathID

    def return_levelpaths_segs(self):
        all_lpids = list()
        for cseg in self.allSegs:
            all_lpids.append(self.allSegs[cseg].levelpathID)
        all_lpids = list(set(all_lpids))
        self.levelpaths_segments = {clpid: [] for clpid in all_lpids}
        # return a dictionary, keyed by segment, of levelpathIDs
        for cseg in self.allSegs:
            self.levelpaths_segments[self.allSegs[cseg].levelpathID].append(cseg)


    def accumulate_same_levelpathID(self, LevelPathdata, COMIDdata, FragIDdata, SFRdata, CELLdata):
        """
        method to add lengths and weight widths and slopes
        for parts of a stream that have the same levelpathID
        within a cell and within the same SFR segment

        """

        #use flag from SFRdata to determine which elevation to use
        #for the reaches

        elevflag = SFRdata.elevflag
        if elevflag == 'DEM':
            elevattr = 'DEM_elev_mean'
            slopeattr = 'DEM_slope'
        elif elevflag == 'smoothed_DEM':
            elevattr = 'smoothed_DEM_elev_mean'
            slopeattr = 'smoothed_DEM_slope'
        elif elevflag == 'elevation_contours':
            elevattr = 'interpolated_contour_elev_mean'
            slopeattr = 'interpolated_contour_slope'
        elif elevflag == 'NHDPlus':
            elevattr = 'NHDPlus_elev_mean'
            slopeattr = 'NHDPlus_slope'
        else:
            raise(BadElevChoice(elevflag))

        #elevattr = 'maxsmoothelev'
        #slopeattr = 'slope'
        #CHK = open('check.out','w')
        #CHK.write('segment, cell, lpID ** segment, cell, lpID\n')
        for segment in self.allSegs.iterkeys():
            rch = 0
            seglpid = self.allSegs[segment].levelpathID
            #get the levelpath ID for the segment of interest
            #and the first and last cell; see if either overlap
            #with the proceeding or next segment and then check to see
            #if the levelpathIDs are the same
            '''
            firstcell = self.allSegs[segment].seg_cells[0]
            lastcell = self.allSegs[segment].seg_cells[-1]
            outseg = self.allSegs[segment].outseg
            for inseg in self.allSegs[segment].inseg:
                if firstcell == self.allSegs[inseg].seg_cells[-1]:
                    CHK.write('topoverlap {0}, {1}, {2} ** {3}, {4}, {5}'.format(segment,
                                                                                 firstcell,
                                                                                 seglpid,
                                                                                 inseg,
                                                                                 self.allSegs[inseg].seg_cells[-1],
                                                                                 self.allSegs[inseg].levelpathID))
                    if self.allSegs[inseg].levelpathID == seglpid:
                        CHK.write(" overlap split needed!\n")
                    else:
                        CHK.write("\n")
            if outseg > 0 and outseg < 999999:
                if lastcell == self.allSegs[outseg].seg_cells[0]:
                    CHK.write('endoverlap  {0}, {1}, {2} ** {3}, {4}, {5}\n'.format(segment,
                                                                                    lastcell,
                                                                                    seglpid,
                                                                                    outseg,
                                                                                    self.allSegs[outseg].seg_cells[0],
                                                                                    self.allSegs[outseg].levelpathID))
                    if self.allSegs[outseg].levelpathID == seglpid:
                        CHK.write(" overlap split needed!\n")
                    else:
                        CHK.write("\n")
            '''
            for localcell in self.allSegs[segment].seg_cells:
                rch += 1        #seg_cells attribute is ordered...
                tt = 0.
                ww = 0.
                el = 0.
                ws = 0.
                #loop over fids and combine/assign to reach if levelpathID matches
                #the levelpath ID for the segment of interest; other levelpathIDs will
                #be combined and assigned during their segment - doesn't subdivide
                #at confluences...
                knt = 0
                for cFragID in FragIDdata.cellnum_FragID[localcell]:
                    ccomid = FragIDdata.allFragIDs[cFragID].comid
                    clpID = COMIDdata.allcomids[ccomid].levelpathID
                    if clpID == seglpid:
                        knt = knt+1
                        tt += FragIDdata.allFragIDs[cFragID].lengthft
                        ww += COMIDdata.allcomids[ccomid].est_width * FragIDdata.allFragIDs[cFragID].lengthft
                        try: #To handle segments that were not attributed from the DEM --> short-term workaround
                            el += getattr(FragIDdata.allFragIDs[cFragID], elevattr)
                            ws += getattr(FragIDdata.allFragIDs[cFragID], slopeattr)*FragIDdata.allFragIDs[cFragID].lengthft
                        except TypeError:
                            el += 9999
                            ws += 9999


                eff_length = tt
                eff_slope = 99999.
                if tt > 0:
                    eff_width = ww/tt
                    eff_slope = ws/tt
                else:
                    eff_width = 99999.
                    eff_slope = 99999.
                if knt > 0:
                    elevreach = el/knt
                else:
                    elevreach = 99999.
                self.allSegs[segment].seg_reaches[rch] = SFRReachProps(localcell, eff_length, eff_width,
                                                            eff_slope, elevreach, SFRdata.bedthick,
                                                            SFRdata.bedK)
                #if its a structured grid, add row and column to seg_reaches
                #the if for unstructured will be added in the future
                self.allSegs[segment].seg_reaches[rch].row = CELLdata.allcells[localcell].row
                self.allSegs[segment].seg_reaches[rch].column = CELLdata.allcells[localcell].column



        #CHK.close()


class SFRpreproc:
    def __init__(self, SFRdata):
        self.indata = SFRdata
        ts = time.time()
        st_time = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        self.ofp = open(os.path.join(self.indata.working_dir, 'SFR_preproc.log'), 'w')
        self.ofp.write('SFR_preproc log.')
        self.ofp.write('\n' + '#' * 25 + '\nStart Time: {0:s}\n'.format(st_time) + '#' * 25 + '\n')


    def clip_and_join_attributes(self, SFRoperations):
        indat = self.indata

        print "Clip original NHD flowlines to model domain..."
        # this creates a new file so original dataset untouched
        arcpy.Clip_analysis(indat.Flowlines_unclipped,
                            indat.MFdomain,
                            indat.Flowlines)

        print "Make a lines only version of the model domain outline..."
        arcpy.FeatureToLine_management(indat.MFdomain, indat.MFdomain[:-4] + '_outline.shp')

        print "joining PlusflowVAA and Elevslope tables to NHD Flowlines..."
        # copy original tables and import copies
        print "importing tables..."
        arcpy.MakeTableView_management(indat.Elevslope, "Elevslope")
        arcpy.MakeTableView_management(indat.PlusflowVAA, "PlusflowVAA")

        # delete all unneeded fields
        fields2keep = ["comid",
                     "divergence",
                     "lengthkm",
                     "thinnercod",
                     "maxelevsmo",
                     "minelevsmo",
                     "hydroseq",
                     "uphydroseq",
                     "dnhydroseq",
                     "reachcode",
                     "streamorde",
                     "arbolatesu",
                     "fcode",
                     "levelpathI",
                     "uplevelpat",
                     "dnlevelpat"]
        fields2keep = [x.lower() for x in fields2keep]
        self.ofp.write('Joining {0:s} with {1:s}: fields kept:\n'.format(indat.Elevslope, indat.Flowlines))
        self.ofp.write('%s\n' % ('\n'.join(fields2keep)))
        print "\nkeeping: %s fields;\ndeleting the rest:" % (' '.join(fields2keep))
        Join = True  # whether or not PlusflowVAA and Elevslope need to be joined
        for table in ["PlusflowVAA", "Elevslope"]:
            fields2delete = []
            allfields = arcpy.ListFields(table)

            for field in allfields:
                name = field.name
                namel = name.lower()
                if namel in fields2keep:
                    if table == "PlusflowVAA" and namel == "maxelevsmo":
                        Join = False
                    continue
                elif namel == "oid":
                    continue
                else:
                    fields2delete.append(field.name)
                    print "\r{} ".format(field.name),
            if len(fields2delete) > 0:
                arcpy.DeleteField_management(table, fields2delete)
        # Join PlusflowVAA and Elevslope to Flowlines
        arcpy.MakeFeatureLayer_management(indat.Flowlines, "Flowlines")

        # If not already, permanently join PlusflowVAA and Elevslope, then to Flowlines
        if Join:
            print "\nJoining Elevslope to PlusflowVAA...\n"
            SFRoperations.getfield("PlusflowVAA", "comid", "comid1")
            SFRoperations.getfield("Elevslope", "comid", "comid2")
            arcpy.JoinField_management("PlusflowVAA",
                                       SFRoperations.joinnames['comid1'],
                                       "Elevslope",
                                       SFRoperations.joinnames['comid2'])
        else:
            print "PlusflowVAA and Elevslope already joined from previous run..."

        print "Joining PlusflowVAA to NHDFlowlines...\n"
        SFRoperations.getfield("Flowlines", "comid", 'comid1')

        # join to Flowlines, keeping only common
        SFR_arcpy.general_join(indat.Flowlines,
                               "Flowlines",
                               SFRoperations.joinnames['comid1'],
                               indat.PlusflowVAA,
                               "comid",
                               True)

        # reopen flowlines as "Flowlines" --> clunky a bit to save and reopen, but must do so
        arcpy.MakeFeatureLayer_management(indat.Flowlines, "Flowlines")

        print "\n"
        self.ofp.write('\n' + 25*'#' + '\nRemoving segments with no elevation information, and with ThinnerCod = -9..\n')
        print "Removing segments with no elevation information, and with ThinnerCod = -9..."

        SFRoperations.getfield("Flowlines", "thinnercod", "ThinnerCod")
        SFRoperations.getfield("Flowlines", "maxelevsmo", "MaxEl")
        SFRoperations.getfield("Flowlines", "comid", "comid")
        FLtable = arcpy.UpdateCursor("Flowlines")
        zerocount = 0
        tcount = 0
        for segments in FLtable:
            if segments.getValue(SFRoperations.joinnames['MaxEl']) == 0 or \
                            segments.getValue(SFRoperations.joinnames['MaxEl']) == -9998:
                print "%d no elevation data" % segments.getValue(SFRoperations.joinnames['comid'])
                self.ofp.write("%d no elevation data\n" % segments.getValue(SFRoperations.joinnames['comid']))
                FLtable.deleteRow(segments)
                zerocount += 1
            elif segments.getValue(SFRoperations.joinnames['ThinnerCod']) == -9:
                print "%d ThinnerCod=-9" % segments.getValue(SFRoperations.joinnames['comid'])
                self.ofp.write("%d ThinnerCod=-9\n" % segments.getValue(SFRoperations.joinnames['comid']))
                FLtable.deleteRow(segments)
                tcount += 1

        # update the "node" field in indat.MFgrid
        # if there is a field with unique values, assume it's ok
        # otherwise delete the column if it already exists
        # NB --> only evaluating the first 100 rows...
        if not indat.node_attribute:
            print 'verifying that there is a "cellnum" field in {0:s}'.format(indat.MFgrid)
            hascellnum = False
            cellnumunique = 0
            MFgridflds = arcpy.ListFields(indat.MFgrid)
            for cfield in MFgridflds:
                if cfield.name.lower() == indat.node_attribute:
                    hascellnum = True
            if hascellnum:
                # now check to see that there are unique values in cell
                cursor = arcpy.SearchCursor(indat.MFgrid)
                cellvals = []
                crows = 0
                for row in cursor:
                    if crows > 100:
                        break
                    else:
                        crows += 1
                        cellvals.append(row.getValue(indat.node_attribute))
                cellnumunique = len(set(cellvals))
                del cellvals
                del row
                del cursor

            if cellnumunique > 1:
                print '"cellnum" field in place with unique values in {0:s}'.format(indat.MFgrid)
            else:
                for fld in arcpy.ListFields(indat.MFgrid):
                    if fld == indat.node_attribute:
                        arcpy.DeleteField_management(indat.MFgrid, indat.node_attribute)
                indat.node_attribute = 'cellnum' # PFJ: specify field name as cellnum
                arcpy.AddField_management(indat.MFgrid, indat.node_attribute, 'LONG')
                calcexpression = '((!row!-1)*{0:d}) + !column!'.format(indat.NCOL)
                arcpy.CalculateField_management(indat.MFgrid, indat.node_attribute, calcexpression, 'PYTHON')
                print 'updated "cellnum" field in {0:s}'.format(indat.MFgrid)




        print "...removed %s segments with no elevation data" %(zerocount)
        self.ofp.write("...removed %s segments with no elevation data\n" %(zerocount))
        print "...removed %s segments with ThinnerCod = -9\n" %(tcount)
        self.ofp.write("...removed %s segments with ThinnerCod = -9\n" %(tcount))
        print ("Performing spatial join (one-to-many) of NHD flowlines to model grid to get river cells..." +
               "(this step may take several minutes or more)\n")
        arcpy.SpatialJoin_analysis(indat.MFgrid, "Flowlines",
                                   indat.CELLS,
                                   "JOIN_ONE_TO_MANY",
                                   "KEEP_COMMON")

        # add in cellnum field for backwards compatibility
 #       arcpy.AddField_management(indat.CELLS, "CELLNUM", "LONG")
 #       arcpy.CalculateField_management(indat.CELLS, "CELLNUM", "!node!", "PYTHON")

        print "Dissolving river cells on cell number to isolate unique cells...\n"
        #SFRoperations.getfield(indat.CELLS, "cellnum", "cellnum")
        SFRoperations.getfield(indat.CELLS, indat.node_attribute, indat.node_attribute)
        arcpy.Dissolve_management(indat.CELLS, indat.CELLS_DISS, SFRoperations.joinnames[indat.node_attribute])

        print "Exploding NHD segments to grid cells using Intersect and Multipart to Singlepart..."
        arcpy.Intersect_analysis([indat.CELLS_DISS, "Flowlines"], os.path.join(indat.working_dir, "tmp_intersect.shp"))
        arcpy.MultipartToSinglepart_management(os.path.join(indat.working_dir, "tmp_intersect.shp"), indat.intersect)
        print "\n"
        print "Adding in stream geometry"
        #set up list and dictionary for fields, types, and associated commands
        fields = ('X_start', 'Y_start', 'X_end', 'Y_end', 'LengthFt')
        types = {'X_start': 'DOUBLE',
                 'Y_start': 'DOUBLE',
                 'X_end': 'DOUBLE',
                 'Y_end': 'DOUBLE',
                 'LengthFt': 'DOUBLE'}
        commands = {'X_start': "float(!SHAPE.firstpoint!.split()[0])",
                    'Y_start': "float(!SHAPE.firstpoint!.split()[1])",
                    'X_end': "float(!SHAPE.lastpoint!.split()[0])",
                    'Y_end': "float(!SHAPE.lastpoint!.split()[1])",
                    'LengthFt': "float(!SHAPE.length!)"}

        #add fields for start, end, and length
        for fld in fields:
            arcpy.AddField_management(indat.intersect, fld, types[fld])

        #calculate the fields
        for fld in fields:
            print "\tcalculating %s(s)..." % (fld)
            arcpy.CalculateField_management(indat.intersect, fld, commands[fld], "PYTHON")
        self.ofp.write('\n' + 25*'#' + '\nRemoving reaches with lengths less than or equal to %s...\n' % indat.reach_cutoff)
        print "\nRemoving reaches with lengths less than or equal to %s..." % indat.reach_cutoff
        SFRoperations.getfield(indat.intersect, "comid", "comid")
        #SFRoperations.getfield(indat.intersect, "cellnum", "cellnum")
        SFRoperations.getfield(indat.intersect, indat.node_attribute, indat.node_attribute)
        SFRoperations.getfield(indat.intersect, "lengthft", "Length")
        table = arcpy.UpdateCursor(indat.intersect)
        count = 0
        for reaches in table:
            if reaches.getValue(SFRoperations.joinnames["Length"]) <= indat.reach_cutoff:
                print "segment: %d, cell: %s, length: %s" % (reaches.getValue(SFRoperations.joinnames["comid"]),
                                                            reaches.getValue(SFRoperations.joinnames[indat.node_attribute]),
                                                            reaches.getValue(SFRoperations.joinnames["Length"]))
                self.ofp.write("segment: %d, cell: %s, length: %s\n"
                          %(reaches.getValue(SFRoperations.joinnames["comid"]),
                            reaches.getValue(SFRoperations.joinnames[indat.node_attribute]),
                            reaches.getValue(SFRoperations.joinnames["Length"])))
                table.deleteRow(reaches)
                count += 1
        print "removed %s reaches with lengths <= %s\n" % (count, indat.reach_cutoff)
        self.ofp.write("removed %s reaches with lengths <= %s\n" % (count, indat.reach_cutoff))
        print "removing cells corresponding to those reaches..."

        # temporarily join river_cells_dissolve to river explode; record nodes with no elevation information
        arcpy.MakeFeatureLayer_management(indat.CELLS_DISS, "river_cells_dissolve")
        arcpy.MakeFeatureLayer_management(indat.intersect, "river_explode")
        arcpy.AddJoin_management("river_cells_dissolve",
                                 indat.node_attribute,
                                 "river_explode",
                                 indat.node_attribute,
                                 "KEEP_ALL")  # this might not work as-in in stand-alone mode
        SFRoperations.getfield("river_cells_dissolve", indat.node_attribute, indat.node_attribute)
        SFRoperations.getfield("river_cells_dissolve", "maxelevsmo", "maxelevsmo")
        table = arcpy.SearchCursor("river_cells_dissolve")
        nodes2delete = []
        for row in table:
            if row.isNull(SFRoperations.joinnames["maxelevsmo"]):
                nodes2delete.append(row.getValue(SFRoperations.joinnames[indat.node_attribute]))
        # this RemoveJoin operation was causing an arcpy error (000800: The value is not a member of | river_explode)
        # and may not be necessary
        #arcpy.RemoveJoin_management("river_cells_dissolve", "river_explode")

        # preserve the FragID code as FragID for all future reference
        arcpy.AddField_management(indat.intersect, "FragID", "LONG")
        arcpy.CalculateField_management(indat.intersect, "FragID", "!FID!", "PYTHON")

        # remove cellnums with no elevation information from river_explode
        self.ofp.write('\n' + 25*'#' + '\nRemoving nodes with no elevation information from river_explode\n')
        print '\nRemoving nodes with no elevation information from river_explode'
        SFRoperations.getfield(indat.CELLS_DISS, indat.node_attribute, indat.node_attribute)
        table = arcpy.UpdateCursor(indat.CELLS_DISS)
        count = 0
        for cells in table:
            if cells.getValue(SFRoperations.joinnames[indat.node_attribute]) in nodes2delete:
                print "%d" % (cells.getValue(SFRoperations.joinnames[indat.node_attribute]))
                self.ofp.write('%d\n' % (cells.getValue(SFRoperations.joinnames[indat.node_attribute])))
                table.deleteRow(cells)
                count += 1
        print "\nremoved %s cells\n" % (count)
        self.ofp.write("removed %s cells" % (count))

        print "removing any remaining disconnected reaches..."
        self.ofp.write('\n' + 25*'#' + '\nremoving any remaining disconnected reaches...\n')
        SFRoperations.getfield(indat.intersect, indat.node_attribute, indat.node_attribute)
        table = arcpy.UpdateCursor(indat.intersect)
        count = 0
        for reaches in table:
            if reaches.getValue(SFRoperations.joinnames[indat.node_attribute]) in nodes2delete:
                print "%d" % (reaches.getValue(SFRoperations.joinnames[indat.node_attribute]))
                self.ofp.write('%d\n' % (reaches.getValue(SFRoperations.joinnames[indat.node_attribute])))
                table.deleteRow(reaches)
                count += 1
        if count > 0:
            print "removed %s disconnected reaches" % count
            self.ofp.write("removed %s disconnected reaches\n" % count)
        else:
            print "no disconnected reaches found!"
            self.ofp.write("no disconnected reaches found!\n")
        print "\n"
        print "Done with pre-processing, ready to run intersect!"
        self.ofp.write('\n' + '#' * 25 + '\nDone with pre-processing, ready to run intersect!')
        ts = time.time()
        end_time = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')

        self.ofp.write('\n' + '#' * 25 + '\nEnd Time: %s\n' % end_time + '#' * 25)
        self.ofp.close()

    def intersect_contours(self, SFRdata):

        # set flag for plotting indicating that this was run for contours
        SFRdata.calculated_contour_elevs = True

        # GIS preprocessing to intersect topographic contours with river_explode.shp

        print "\nIntersecting NHD flowlines with elevation contours..."
        # create set of points where elevation contours intersect the streams
        arcpy.Intersect_analysis([SFRdata.ELEVcontours, SFRdata.intersect], SFRdata.Contours_intersect, "ALL", "", "POINT")

        # add field to river_explode that uniquely identifies each FID (but isn't named FID)
        arcpy.AddField_management(SFRdata.intersect, "FragID", "LONG")
        arcpy.CalculateField_management(SFRdata.intersect, "FragID", "!FID!", "PYTHON")

        # add field of 0s to river_explode to designate starting points for route measurements
        # (from_measure_field for routes)
        arcpy.AddField_management(SFRdata.intersect, "FromM", "SHORT")

        # convert stream linework to Arc routes (necessary for LocateFeaturesAlongRoutes_lr)
        # maintain m-data from original NHD, by specifying starting and ending distance fields for each FragID
        # all "FromM" values are 0 (measurements are from the start of each segment, as defined in the m-direction);
        # "LengthFt" is used for the to_measure_field field
        arcpy.CreateRoutes_lr(SFRdata.intersect, "FragID", SFRdata.Routes, "TWO_FIELDS", "FromM", "LengthFt")

        # returns location of contour intersections along each FragID of stream in river_explode,
        # as a distance from the line segment start identified from CreateRoutes_lr
        arcpy.LocateFeaturesAlongRoutes_lr(SFRdata.Contours_intersect, SFRdata.Routes, "FragID", "1 feet",
                                SFRdata.Contours_intersect_distances, "Route POINT fmp", "", "", "", "", "M_DIRECTON")

    def intersect_DEM(self, SFRdata):

        # set flag for plotting indicating that this was run for DEM
        SFRdata.calculated_DEM_elevs = True

        # GIS preprocessing to intersect DEM with river_explode.shp

        print "\n\nIntersecting NHD flowline fragments with DEM..."
        # convert end vertices of river_explode Fragments to points
        arcpy.FeatureVerticesToPoints_management(SFRdata.intersect, SFRdata.intersect_points)

        # add FragID field identifying unique river_explode segments if one doesn't exist
        arcpy.MakeFeatureLayer_management(SFRdata.intersect, "river_explode")
        Fields = [f.name.lower() for f in arcpy.ListFields("river_explode")]
        if "fragid" not in Fields:
            arcpy.AddField_management(SFRdata.intersect, "FragID", "LONG")
            arcpy.CalculateField_management(SFRdata.intersect, "FragID", "!FID!", "PYTHON")

        # extract DEM values at locations of points
        arcpy.sa.ExtractMultiValuesToPoints(SFRdata.intersect_points, [[SFRdata.DEM, "DEM_elev"]])


class SFROperations:
    """
    class to make operations on data
    """
    def __init__(self, SFRdata):
        self.SFRdata = SFRdata
        self.newCOMIDdata = dict()
        self.joinnames = dict()

    def getfield(self, table, joinname, returnedname):
    # get name of field (useful for case issues and appended field names in joined tables, etc)
        Fields = arcpy.ListFields(table)
        joinname = joinname.lower()
        for field in Fields:
            if joinname in field.name.lower():
                joinname = field.name
                break
        self.joinnames[returnedname] = joinname

    def intersect(self):
        """
        merge the NHDplus information with the model grid
        this was originally the intersect.py script
        """
        # check to make sure that preprocessing has been run
        if not os.path.isfile(self.SFRdata.intersect):
            raise NoPreProcessing


        #convert to layers
        arcpy.MakeFeatureLayer_management(self.SFRdata.Flowlines_unclipped, 'nhd_lyr')
        arcpy.MakeFeatureLayer_management(self.SFRdata.MFdomain, 'grid_lyr')
        print "selecting streams that cross model grid boundary..."
        arcpy.SelectLayerByLocation_management('nhd_lyr', 'CROSSED_BY_THE_OUTLINE_OF', 'grid_lyr')
        arcpy.CopyFeatures_management('nhd_lyr', os.path.join(self.SFRdata.working_dir, 'intersect.shp'))

        #copy the model NHD streams to a temp shapefile, find the ends that match
        #the streams that were clipped by the boundary and update the lengthKm,
        #minsmoothelev, maxsmoothelev in the MODLNHD
        if arcpy.Exists(os.path.join(self.SFRdata.working_dir, 'temp.shp')):
            print 'first removing old version of temp.shp'
            arcpy.Delete_management(os.path.join(self.SFRdata.working_dir, 'temp.shp'))
        arcpy.CopyFeatures_management(self.SFRdata.Flowlines, os.path.join(self.SFRdata.working_dir, 'temp.shp'))

        #add fields for start, end, and length to the temp and intersect
        #shapefiles  (use LENKM as field name because temp already has LENGHTKM)
        print "adding fields for start, end and length..."
        shapelist = (os.path.join(self.SFRdata.working_dir, 'temp.shp'), os.path.join(self.SFRdata.working_dir, 'intersect.shp'))
        for shp in shapelist:
            arcpy.AddField_management(shp, 'STARTX', 'DOUBLE')
            arcpy.AddField_management(shp, 'STARTY', 'DOUBLE')
            arcpy.AddField_management(shp, 'ENDX', 'DOUBLE')
            arcpy.AddField_management(shp, 'ENDY', 'DOUBLE')
            arcpy.AddField_management(shp, 'LENKM', 'DOUBLE')

        print "calculating new info for fields..."
        #calculate the fields, convert length to km - projection is in feet
        for shp in shapelist:
            arcpy.CalculateField_management(shp, 'STARTX', "float(!SHAPE.firstpoint!.split()[0])", "PYTHON")
            arcpy.CalculateField_management(shp, 'STARTY', "float(!SHAPE.firstpoint!.split()[1])", "PYTHON")
            arcpy.CalculateField_management(shp, 'ENDX', "float(!SHAPE.lastpoint!.split()[0])", "PYTHON")
            arcpy.CalculateField_management(shp, 'ENDY', "float(!SHAPE.lastpoint!.split()[1])", "PYTHON")
            arcpy.CalculateField_management(shp, 'LENKM', "float(!SHAPE.length@kilometers!)", "PYTHON")

        #go through intersect, identify which end matches COMID in temp
        #find new length in temp and use linear interpolation to get new elev
        #finally put COMID out to <outfile> (named above) to indicate
        #if the cut-off is flowing out of grid or into grid so the
        #routing tables in the final datasets do not route out of grid
        #and back in. Also identifies to user ends of streams that
        #could accept inflow conditions for SFR
        print "fixing routing for streams that cross the grid boundary..."

        eps = self.SFRdata.eps
        # make a list of COMIDs that are found for later reference
        comidlist = []
        intersects = arcpy.SearchCursor(os.path.join(self.SFRdata.working_dir, 'intersect.shp'))
        manual_intervention = 0
        cutoff_comids = []
        for stream in intersects:
            comid = int(stream.COMID)
            query = "COMID = {0:d}".format(comid)
            stx = stream.STARTX
            sty = stream.STARTY
            endx = stream.ENDX
            endy = stream.ENDY
            lenkm = stream.LENKM
            clippedstream = arcpy.SearchCursor('temp.shp', query)
            for clip in clippedstream:
                clstx = clip.STARTX
                clsty = clip.STARTY
                clendx = clip.ENDX
                clendy = clip.ENDY
                cllen = clip.LENKM
                clmaxel = clip.MAXELEVSMO
                clminel = clip.MINELEVSMO
                #find the end that matches
                stdiffx = stx-clstx
                stdiffy = sty-clsty
                enddiffx = endx-clendx
                enddiffy = endy-clendy
                if np.fabs(stdiffx) < eps and np.fabs(stdiffy) < eps:
                    comidlist.append(comid)
                    self.newCOMIDdata[comid] = COMIDPropsForIntersect(
                        comid, 'OUT', clstx, clsty, clendx, clendy, clmaxel,
                        clminel, cllen, lenkm
                    )
                elif np.fabs(enddiffx) < eps and np.fabs(enddiffy) < eps:
                    comidlist.append(comid)
                    self.newCOMIDdata[comid] = COMIDPropsForIntersect(
                        comid, 'IN', clstx, clsty, clendx, clendy, clmaxel,
                        clminel, cllen, lenkm
                    )
                else:
                    manual_intervention += 1
                    if manual_intervention == 1:
                        ofp = open(os.path.join(self.SFRdata.working_dir, 'boundary_manual_fix_issues.txt'), 'w')
                        ofp.write('The following COMIDs identify streams that need manual attention.\n')
                        ofp.write('Fix in the files {0:s} and river_explode.shp. '
                                  'Then rerun main script without preprocessing\n'.format(self.SFRdata.Flowlines))
                        ofp.write('#' * 16 + '\n')
                    cutoff_comids.append(comid)
                    ofp.write('both ends are cut off for comid {0:d}\n'.format(comid))
        if len(cutoff_comids) > 0:
            print 'The following COMIDs have both ends cutoff:'
            ofp.write('The following COMIDs have both ends cutoff:\n')
            for comid in cutoff_comids:
                print comid
                ofp.write('{}\n'.format(comid))

        del intersects, stream
        if manual_intervention:
            ofp.write('#' * 16 + '\n')
            ofp.close()

        # ATL: COMIDs with multiple ends and those with cut off ends
        # (i.e., listed in fix_com_IDs.txt and boundary_manual_fix_issues.txt)
        # appear to be already removed in SFRdata.intersect (i.e. river_explode).
        # Remove these features by clipping SFRdata.Flowlines to SFRdata.intersect
        # output to SFRdata.NHD (instead of simply copying SFRdata.Flowlines, as below)
        #arcpy.Clip_analysis(self.SFRdata.Flowlines, self.SFRdata.intersect, self.SFRdata.NHD)

        #create a new working NHD shapefile incorporating new values just found

        print "\nCreating new shapefile {0:s}".format(self.SFRdata.NHD)
        arcpy.CopyFeatures_management(self.SFRdata.Flowlines, self.SFRdata.NHD)
        intersects = arcpy.UpdateCursor(self.SFRdata.NHD)
        for stream in intersects:
            comid = stream.COMID
            if comid in comidlist:
                stream.LENGTHKM = self.newCOMIDdata[comid].newlength
                stream.MAXELEVSMO = self.newCOMIDdata[comid].newmaxel
                stream.MINELEVSMO = self.newCOMIDdata[comid].newminel
                intersects.updateRow(stream)
        del stream, intersects

        # finally, save out the routing information
        print "Saving routing information to {0:s}".format(os.path.join(self.SFRdata.working_dir, 'boundaryClipsRouting.txt'))
        ofp = open(os.path.join(self.SFRdata.working_dir, 'boundaryClipsRouting.txt'), 'w')
        ofp.write("FROMCOMID,TOCOMID\n")
        for currcomid in self.newCOMIDdata:
            if self.newCOMIDdata[currcomid].inout == 'IN':
                ofp.write("99999,{0:d}\n".format(self.newCOMIDdata[currcomid].comid))
            else:
                ofp.write("{0:d},99999\n".format(self.newCOMIDdata[currcomid].comid))
        ofp.close()
        if manual_intervention:
            print '\nSome manual intervention required:\n' \
                  'See "{}" for details'.format(os.path.join(self.SFRdata.working_dir, 'boundary_manual_fix_issues.txt'))

    def make_rivers_table(self, FragIDdata):
        """
        from assign_and_route -->
        """
        print "\nSetting up the elevation table --> {0:s}".format(self.SFRdata.rivers_table)
        if arcpy.Exists(self.SFRdata.rivers_table):
            arcpy.Delete_management(self.SFRdata.rivers_table)
        arcpy.CreateTable_management(self.SFRdata.working_dir, os.path.split(self.SFRdata.rivers_table)[-1])
        arcpy.AddField_management(self.SFRdata.rivers_table, "OLDFragID", "LONG")
        arcpy.AddField_management(self.SFRdata.rivers_table, self.SFRdata.node_attribute, "LONG")
        arcpy.AddField_management(self.SFRdata.rivers_table, "ELEVMAX", "DOUBLE")
        arcpy.AddField_management(self.SFRdata.rivers_table, "ELEVAVE", "DOUBLE")
        arcpy.AddField_management(self.SFRdata.rivers_table, "ELEVMIN", "DOUBLE")
        segcounter = 0
        rows = arcpy.InsertCursor(self.SFRdata.rivers_table)
        fix_comids_summary = []
        fixcomids_flag = False
        for comid in FragIDdata.allcomids:
            segcounter += 1
            FragIDlist = FragIDdata.comid_FragID[comid]
            start_has_end = dict()
            end_has_start = dict()

            for i in FragIDlist:
                haveend = False
                havestart = False
                for j in FragIDlist:
                    if j != i:
                        diffstartx = FragIDdata.allFragIDs[i].startx - FragIDdata.allFragIDs[j].endx
                        diffstarty = FragIDdata.allFragIDs[i].starty - FragIDdata.allFragIDs[j].endy
                        diffendx = FragIDdata.allFragIDs[i].endx - FragIDdata.allFragIDs[j].startx
                        diffendy = FragIDdata.allFragIDs[i].endy - FragIDdata.allFragIDs[j].starty
                        if np.fabs(diffstartx) < self.SFRdata.rfact and np.fabs(diffstarty) < self.SFRdata.rfact:
                            start_has_end[i] = j
                            haveend = True
                        if np.fabs(diffendx) < self.SFRdata.rfact and np.fabs(diffendy) < self.SFRdata.rfact:
                            end_has_start[i] = j
                            havestart = True
                        if haveend and havestart:
                            break

            #find key in start_has_end that didn't match an end and
            #key in end_has_start that didn't match a start
            numstart = 0
            numend = 0
            startingFragID = []
            endingFragID = []
            for test in FragIDlist:
                if test not in start_has_end:
                    startingFragID.append(test)
                    numstart += 1
                if test not in end_has_start:
                    endingFragID.append(test)
                    numend += 1
            if numstart != 1 or numend != 1:
                if not fixcomids_flag:
                    outfile = open(os.path.join(self.SFRdata.working_dir, "fix_com_IDs.txt"), 'w')
                    fixcomids_flag = True
                outfile.write("numstart =" + str(numstart) + " \n")
                outfile.write("numend = " + str(numend) + " \n")
                outfile.write("starting FragIDs: " + ",".join(map(str, startingFragID)) + "\n")
                outfile.write("ending FragIDs: " + ",".join(map(str, endingFragID)) + "\n")
                outfile.write("manually fix COMID = %d\n" %comid)
                fix_comids_summary.append('%d\n' %comid)
                FragIDdata.noelev[comid] = 1  #set flag
                continue
            orderedFragID = []
            orderedFragID.append(startingFragID[0])
            for i in range(1, len(end_has_start)):
                try:
                    orderedFragID.append(end_has_start[orderedFragID[i-1]])
                except KeyError:
                    raise(MissingFragID('end_has_start', orderedFragID[i-1], FragIDdata.allFragIDs[orderedFragID[i-1]].comid))

            if orderedFragID[-1] != endingFragID[0]:       #  don't repeat the last entry FragID...
                orderedFragID.append(endingFragID[0])
            #total length read through lengthkm didn't always match up
            #to the sum of the lengths of the segments (exactly), sumup total length
            totallength = 0
            for i in range(0, len(orderedFragID)):
                totallength += FragIDdata.allFragIDs[orderedFragID[i]].lengthft
            if totallength == 0:
                exit('check length ft for FragIDs in COMID= %d' % comid)
            FragIDdata.return_smoothelev_comid(comid)
            slope = (FragIDdata.maxelev - FragIDdata.minelev)/totallength
            distance = 0.
            FragIDdata.COMID_orderedFragID[comid] = orderedFragID
            for i in range(0, len(orderedFragID)):
                maxcellrivelev = FragIDdata.maxelev - slope*distance
                distance += FragIDdata.allFragIDs[orderedFragID[i]].lengthft
                mincellrivelev = FragIDdata.maxelev - slope*distance
                avecellrivelev = 0.5*(maxcellrivelev + mincellrivelev)

                # assign NHD elevations to FragID properties
                FragIDdata.allFragIDs[orderedFragID[i]].NHDPlus_elev_min = mincellrivelev
                FragIDdata.allFragIDs[orderedFragID[i]].NHDPlus_elev_max = maxcellrivelev
                FragIDdata.allFragIDs[orderedFragID[i]].NHDPlus_elev_mean = avecellrivelev
                FragIDdata.allFragIDs[orderedFragID[i]].NHDPlus_slope = slope

                row = rows.newRow()
                row.OLDFragID = orderedFragID[i]
                row.CELLNUM = FragIDdata.allFragIDs[orderedFragID[i]].cellnum
                row.ELEVMAX = maxcellrivelev
                row.ELEVAVE = avecellrivelev
                row.ELEVMIN = mincellrivelev
                rows.insertRow(row)
        # write out the summary of comids to fix, then close the outfile
        if len(fix_comids_summary) > 0:
            print '\nSome cells have multiple COMIDs entering and/or leaving.\nSee file "{}"'.format(
                os.path.join(self.SFRdata.working_dir, 'fix_com_IDs.txt'))
            outfile.write('#' * 30 + '\nSummary of COMIDS to fix:\n' +
                          'Delete these COMIDs from river_explode.shp, \n' +
                          'then run CleanupRiverCells.py Then rerun main script without preprocessing\n')
            [outfile.write(line) for line in fix_comids_summary]
            outfile.close()
        del row
        del rows

    def reach_ordering(self, COMIDdata, FragIDdata, LevelPathdata):
        """
        crawl through the hydrosequence values and
        set a preliminary reach ordering file
        """
        indat = self.SFRdata
        SFRseq = 0
        # write out the candidate routing information by hydrosequence ordered
        ofp = open(indat.RCH, 'w')
        ofp.write('CELLNUM, COMID, hydroseq, uphydroseq, dnhydroseq, '
                  'levelpathID, dnlevelpath, SFRseq, localseq\n')
        for currhydroseq in COMIDdata.hydrosequence_sorted:
            ccomid = COMIDdata.hydrosequence_comids[currhydroseq]
            if ccomid not in FragIDdata.noelev.keys():
                # update the levelpathID hydrosequence data
                LevelPathdata.allids[COMIDdata.allcomids[ccomid].levelpathID].ordered_hydrosequence.append(currhydroseq)
                SFRseq += 1
                localseq = 0
                for cFragID in FragIDdata.COMID_orderedFragID[ccomid]:
                    localseq += 1
                    ofp.write('{0:d},{1:d},{2:d},{3:d},{4:d},{5:d},{6:d},{7:d},{8:d}\n'.format(
                    FragIDdata.allFragIDs[cFragID].cellnum,
                    ccomid,
                    COMIDdata.allcomids[ccomid].hydrosequence,
                    COMIDdata.allcomids[ccomid].uphydrosequence,
                    COMIDdata.allcomids[ccomid].downhydrosequence,
                    COMIDdata.allcomids[ccomid].levelpathID,
                    LevelPathdata.allids[COMIDdata.allcomids[ccomid].levelpathID].down_levelpathID,
                    SFRseq,
                    localseq
                    ))
        ofp.close()

        # calculate and list unique cellnums and FragIDs in downstream order by levelpathID
        for clevelpathid in LevelPathdata.level_ordered:
            LevelPathdata.allids[clevelpathid].ordered_hydrosequence = \
                sorted(list(set(LevelPathdata.allids[clevelpathid].ordered_hydrosequence)), reverse=True)
            for currhydroseq in LevelPathdata.allids[clevelpathid].ordered_hydrosequence:
                ccomid = COMIDdata.hydrosequence_comids[currhydroseq]
                if ccomid not in FragIDdata.noelev.keys():
                    for cFragID in FragIDdata.COMID_orderedFragID[ccomid]:
                        if FragIDdata.allFragIDs[cFragID].cellnum not in LevelPathdata.allids[clevelpathid].ordered_cellnums:
                            LevelPathdata.allids[clevelpathid].ordered_cellnums.append(
                            FragIDdata.allFragIDs[cFragID].cellnum
                            )
                        LevelPathdata.allids[clevelpathid].ordered_FragIDs.append(cFragID)

    def assign_layers(self, SFRdata):
        """
        from Assign_Layers.py --> uses model layer information to assign SFR cells to proper layers
        """

        BotcorPDF = os.path.join(SFRdata.working_dir, "Corrected_Bottom_Elevations.pdf")  # PDF file showing original and corrected bottom elevations
        Layerinfo = os.path.join(SFRdata.working_dir, "SFR_layer_assignments.txt")  # text file documenting how many reaches are in each layer as assigned
        DX, DY, NLAY, NROW, self.NCOL, i = disutil.read_meta_data(self.SFRdata.MFdis)
        print "\nRead in model grid top elevations from {0:s}".format(SFRdata.MFdis)
        topdata, i = disutil.read_nrow_ncol_vals(SFRdata.MFdis, SFRdata.NROW, SFRdata.NCOL, np.float, i)
        print "Read in model grid bottom layer elevations from {0:s}".format(SFRdata.MFdis)
        bots = np.zeros([SFRdata.NLAY, SFRdata.NROW, SFRdata.NCOL])
        for clay in np.arange(SFRdata.NLAY):
            print 'reading layer {0:d}'.format(clay+1)
            bots[clay, :, :], i = disutil.read_nrow_ncol_vals(SFRdata.MFdis, SFRdata.NROW, SFRdata.NCOL, np.float, i)
        SFRinfo = np.genfromtxt(SFRdata.MAT1, delimiter=',', names=True, dtype=None)

        print '\nNow assiging stream cells to appropriate layers...'
        below_bottom = open(os.path.join(SFRdata.working_dir, 'below_bot.csv'), 'w')
        below_bottom.write('SFRbot,ModelBot,Land_surf,cellnum,segment\n')
        below_bot_adjust = defaultdict()  # list row,column locations where SFR goes below model bottom
        nbelow = 0
        New_Layers = []
        ofp_bottom_warnings = open(os.path.join(SFRdata.working_dir, 'bottom_warnings.dat'), 'w')
        for i in range(len(SFRinfo)):
            r = SFRinfo['row'][i]
            c = SFRinfo['column'][i]
            STOP = SFRinfo['top_streambed'][i]
            cellbottoms = list(bots[:, r-1, c-1])
            SFRbot = STOP - SFRdata.bedthick - SFRdata.buff
            for b in range(SFRdata.NLAY):
                if SFRbot < cellbottoms[b]:
                    if b+1 >= SFRdata.NLAY:
                        warnstring1 = 'Streambottom elevation={0:f}, Model bottom={1:f} at ' \
                              'row {2:d}, column {3:d}, cellnum {4:d}'.format(
                              SFRbot, cellbottoms[-1], r, c, (r-1)*SFRdata.NCOL + c)
                        warnstring2 = 'Land surface is {0:f}'.format(topdata[r-1, c-1])
                        ofp_bottom_warnings.write(warnstring1  + warnstring2 + '\n')
                        print warnstring1
                        print warnstring2
                        below_bottom.write('{0:f},{1:f},{2:f},{3:d},{4:d}\n'.format(
                            SFRbot, cellbottoms[-1], topdata[r-1, c-1], (r-1)*SFRdata.NCOL+c, SFRinfo['segment'][i]))
                        below_bot_adjust[(r-1, c-1)] = cellbottoms[-1] - SFRbot  # diff between SFR bottom and model bot
                        nbelow += 1
                        New_Layers.append(b+1)
                else:
                    New_Layers.append(b+1)
                    break
        ofp_bottom_warnings.close()
        below_bottom.close()
        New_Layers = np.array(New_Layers)
        bots_orig = bots[-1, :, :].copy() # keep a copy of non-changed bottom elevations for plotting

        # create a new array of bottom elevations with dimensions like topdata
        if SFRdata.Lowerbot:
            print "\n\nAdjusting model bottom to accommodate SFR cells that were below bottom"
            print "see {0:s}\n".format(BotcorPDF)
            below_bot_adjust_tuples = below_bot_adjust.keys()
            #for (r, c) in below_bot_adjust_tuples:
            #    bots[-1, r, c] -= below_bot_adjust[(r, c)]

            # PFJ:  This for loop was generated by MNF to improve indexing to py and MF arrays.
            for r in range(1, SFRdata.NROW + 1):
                for c in range(1, SFRdata.NCOL + 1):
                    if (r, c) in below_bot_adjust_tuples:
                        bots[-1, r-1, c-1] -= below_bot_adjust[(r, c)]

            '''
            for r in range(SFRdata.NROW):
                for c in range(SFRdata.NCOL):
                    if r == 142 and c == 63:
                        j=2
                    if (r, c) in below_bot_adjust_tuples:
                        bots[-1, r, c] -= below_bot_adjust[(r, c)]
            '''
            outarray = os.path.join(SFRdata.working_dir, 'SFR_Adjusted_bottom_layer.dat')

            with file(outarray, 'w') as outfile:
                for slice_2d in bots:
                    np.savetxt(outfile, slice_2d, fmt='%.5e')
            outfile.close()

            # show bottom before and after corrections
            outpdf = PdfPages(BotcorPDF)
            befaft = ['before adjustment', 'after adjustment', 'adjustments']
            diff = bots_orig - bots[-1, :, :]
            for i, mat in enumerate([bots_orig, bots[-1, :, :], diff]):
                plt.figure()
                plt.imshow(mat)
                plt.colorbar()
                plt.title('Bottom of model {0:s}'.format(befaft[i]))
                outpdf.savefig()
            outpdf.close()

        # histogram of layer assignments
        freq = np.histogram(list(New_Layers), range(SFRdata.NLAY + 2)[1:])

        # writeout info on Layer assignments
        ofp = open(Layerinfo, 'w')
        ofp.write('Layer\t\tNumber of assigned reaches\n')
        print '\nLayer assignments:'
        for i in range(SFRdata.NLAY):
            ofp.write('{0:d}\t\t{1:d}\n'.format(freq[1][i], freq[0][i]))
            print '{0:d}\t\t{1:d}\n'.format(freq[1][i], freq[0][i])
        ofp.close()

        # write new SFRmat1 file
        print "Writing output to %s..." %(SFRdata.MAT1)
        ofp = open(SFRdata.MAT1,'w')
        ofp.write(','.join(SFRinfo.dtype.names)+'\n')

        SFRinfo['layer'] = New_Layers
        for i in range(len(SFRinfo)):
            line = list(SFRinfo[i])
            line = ','.join(map(str, line))
            ofp.write('{0:s}\n'.format(line))
        ofp.close()

        if not SFRdata.Lowerbot:
            if nbelow > 0:
                print "Warning {0:d} SFR streambed bottoms were below the model bottom. See below_bots.csv".format(
                    nbelow)


class ElevsFromContours:
    def __init__(self, SFRdata):
        self.intersect_dist_table = SFRdata.Contours_intersect_distances
        self.elevs_edited = []
        self.in_between_contours = []
        self.comids_with_no_contour_intersects = []
        self.verbose = False # prints interpolation information to screen for debugging/troubleshooting

    def interpolate_elevs(self, FragIDdata, start_fid, end_fid):
        # interpolate max/min elevation values for each FragID between start and end fids

        # first setup indexing
        # start_fid is between first reach in current comid and end_fid
        # (for cases where up/down contour is outide of current comid, it was already reset to the first or last reach)
        dist = 0
        end_reach = FragIDdata.COMID_orderedFragID[self.current_comid].index(end_fid)
        start_reach = FragIDdata.COMID_orderedFragID[self.current_comid].index(start_fid)

        # compute max elevation in first FragID
        if end_fid == FragIDdata.COMID_orderedFragID[self.current_comid][-1]:
            # downstream contour was below current comid
            dist += self.downstream_dist
            FragIDdata.allFragIDs[end_fid].interpolated_contour_elev_min = self.end_elev + self.slope * dist
            dist += FragIDdata.allFragIDs[end_fid].lengthft
            FragIDdata.allFragIDs[end_fid].interpolated_contour_elev_max = self.end_elev + self.slope * dist

        else:
            dist = np.min(FragIDdata.allFragIDs[end_fid].contour_elev_distance)
            FragIDdata.allFragIDs[end_fid].interpolated_contour_elev_max = self.end_elev + self.slope * dist
        if self.verbose:
            print 'fid:%s min/max elevs %s %s' % (end_fid, FragIDdata.allFragIDs[end_fid].interpolated_contour_elev_min, FragIDdata.allFragIDs[end_fid].interpolated_contour_elev_max)

        # compute max and min elevations in subsequent FragIDs
        ind = end_reach
        for FragID in FragIDdata.COMID_orderedFragID[self.current_comid][start_reach:end_reach][::-1]:
            ind -= 1
            previous_fid = FragIDdata.COMID_orderedFragID[self.current_comid][ind + 1]
            FragIDdata.allFragIDs[FragID].interpolated_contour_elev_min = FragIDdata.allFragIDs[previous_fid].interpolated_contour_elev_max
            dist += FragIDdata.allFragIDs[FragID].lengthft
            FragIDdata.allFragIDs[FragID].interpolated_contour_elev_max = self.end_elev + self.slope * dist
            if self.verbose:
                print 'FragID:%s min/max elevs %s %s' % (FragID, FragIDdata.allFragIDs[FragID].interpolated_contour_elev_min, FragIDdata.allFragIDs[FragID].interpolated_contour_elev_max)


    def get_dist_slope(self, comid, FragIDdata, COMIDdata):

            for FragID in self.reachlist:
                if self.verbose:
                    print "comid %s, FragID %s" %(comid, FragID)

                # looking downstream for first contour below current comid
                if self.downstream:

                    # at end of comid, check for outlet
                    if FragID == self.reachlist[-1]:
                        # another downcomid exists, continue
                        try:
                            COMIDdata.allcomids[COMIDdata.allcomids[comid].to_comid[0]]
                            # additional test to see if COMID listed in FragIDdata
                            # if not, is in the list of COMIDs that were flagged earlier as having multiple starts/ends
                            # or the stream leaves the grid
                            FragIDdata.COMID_orderedFragID[COMIDdata.allcomids[comid].to_comid[0]]
                        # at outlet; assign downstream elev using NHD
                        except KeyError:
                            self.downstream_dist += FragIDdata.allFragIDs[FragID].lengthft
                            self.downstream_contour_comid = comid
                            self.end_elev = FragIDdata.allFragIDs[FragID].minsmoothelev
                            self.downstream = False
                            break
                    # found contour; assign downstream elev
                    elif FragID in self.elevs_edited:
                        self.downstream_dist += np.min(FragIDdata.allFragIDs[FragID].contour_elev_distance)
                        self.downstream_contour_comid = comid
                        self.end_elev = np.max(FragIDdata.allFragIDs[FragID].contour_elev)
                        self.downstream = False
                        break
                    else:
                        self.downstream_dist += FragIDdata.allFragIDs[FragID].lengthft

                # moving upstream
                if self.upstream:

                    # most upstream reach in comid, check for headwaters or edge of grid
                    # if from_comid is 999999, headwater; or, if no from_comid, stream originates outside of grid
                    # in both cases, use NHD to set end elevation, unless slope is negative, then set slope to 0
                    if FragID == FragIDdata.COMID_orderedFragID[comid][0]:

                        # in all cases, add distance first
                        self.upstream_dist += FragIDdata.allFragIDs[FragID].lengthft

                        # for multiple from_comids (confluence)
                        for upid in COMIDdata.allcomids[comid].from_comid:
                            try:
                                COMIDdata.allcomids[upid]
                                if upid == 999999:
                                    headwater = True
                                else:
                                    # if an upstream segment is encountered that isn't a headwater,
                                    # continue interpolating into that segment
                                    headwater = False
                                    break
                            except KeyError:
                                headwater = True

                        if headwater:

                            # set start_elev and calculate slope
                            self.start_elev = FragIDdata.allFragIDs[FragID].NHDPlus_elev_max
                            self.slope = (self.start_elev - self.end_elev) / self.upstream_dist
                            if self.slope < 0:
                                self.slope = 0

                            # set end_fid (default is last FragID in current_comid)
                            # if downstream contour found in current comid:
                            #if self.downstream_contour_comid == self.current_comid:
                                #self.end_FragID = self.reachlist[0]

                            # set start_FragID
                            # if upstream contour found in current comid:
                            if comid == self.current_comid:
                                self.start_FragID = FragID
                            # or if upstream contour was found above current comid
                            else:
                                self.start_FragID = FragIDdata.COMID_orderedFragID[self.current_comid][0]

                            # go to interpolation
                            break

                        else:
                            continue

                    # found a contour
                    elif FragID in self.elevs_edited:

                        # calculate distance; reset downstream distance
                        self.upstream_dist += FragIDdata.allFragIDs[FragID].lengthft - np.max(FragIDdata.allFragIDs[FragID].contour_elev_distance)

                        # set start_elev and calculate slope
                        self.start_elev = np.min(FragIDdata.allFragIDs[FragID].contour_elev)
                        self.slope = (self.start_elev - self.end_elev) / self.upstream_dist

                        # set end_FragID (default is last FragID in current_comid)
                        # if downstream contour found in current comid:
                        #if self.downstream_contour_comid == self.current_comid:
                            #self.end_FragID = self.reachlist[0]

                        # set start_FragID
                        # if upstream contour found in current comid:
                        if comid == self.current_comid:
                            self.start_FragID = FragID
                        # or if upstream contour was found above current comid
                        else:
                            self.start_FragID = FragIDdata.COMID_orderedFragID[self.current_comid][0]

                        # break out of get_dist_slope and go to interpolation
                        break

                    # no contour found, keep moving upstream
                    else:
                        self.upstream_dist += FragIDdata.allFragIDs[FragID].lengthft


    def get_contour_intersections(self, FragIDdata, COMIDdata):


        #The slow way with Arc
        # loop through rows in contours intersect table
        distances = arcpy.SearchCursor(self.intersect_dist_table)

        print "getting elevations at topographic contour intersections..."
        for row in distances:
            comid = int(row.COMID)
            cellnum = int(row.cellnum)

            # loop through FragIDs in comid, if cellnum matches, get contour elevation and distance
            try:
                for FragID in FragIDdata.COMID_orderedFragID[comid]:
                    if FragIDdata.allFragIDs[FragID].cellnum == cellnum:

                        # multiple contours in a FragID are appended to a list
                        # up and downstream slopes from an intersected FragID are calculated using
                        # respective max and min intersected contour values
                        FragIDdata.allFragIDs[FragID].contour_elev.append(float(row.ContourEle))
                        FragIDdata.allFragIDs[FragID].contour_elev_distance.append(float(row.fmp))

                        self.elevs_edited.append(FragID)
            # a key error means that the COMID was removed earlier in the SFR setup process
            except KeyError:
                continue

        # build list of FragIDs that do not intersect any contours
        for FragID in FragIDdata.allFragIDs.keys():
            if FragID not in self.elevs_edited:
                self.in_between_contours.append(FragID)
        '''
        # faster way with python
        distances = np.genfromtxt(self.intersect_dist_table, names=True, delimiter=',', dtype=None)

        for row in distances:
            comid = int(row['COMID'])
            cellnum = int(row['node'])

            # loop through FragIDs in comid, if cellnum matches, get contour elevation and distance
            for FragID in FragIDdata.COMID_orderedFragID[comid]:
                if FragIDdata.allFragIDs[FragID].cellnum == cellnum:
                    # multiple contours in a FragID are appended to a list
                    # up and downstream slopes from an intersected FragID are calculated using
                    # respective max and min intersected contour values
                    FragIDdata.allFragIDs[FragID].contour_elev.append(float(row['ContourEle']))
                    FragIDdata.allFragIDs[FragID].contour_elev_distance.append(float(row['fmp']))

                    self.elevs_edited.append(FragID)
        '''


    def assign_elevations_to_FragID(self, FragIDdata, COMIDdata):

        print "interpolating elevations to FragIDs..."
        knt = 0
        for comid in FragIDdata.COMID_orderedFragID.keys():
            knt += 1
            print "\rCOMID: {0} ({1} of {2})".format(comid, knt, len(FragIDdata.COMID_orderedFragID.keys())),
            if self.verbose:
                print "comid number %s" % (knt)
            from_comid = COMIDdata.allcomids[comid].from_comid
            to_comid = COMIDdata.allcomids[comid].to_comid
            self.current_comid = comid

            # within each COMID, iterate going upstream (reversed)
            self.start_elev = None
            self.end_elev = None
            self.start_reach = 0
            self.start_FragID = FragIDdata.COMID_orderedFragID[comid][-1]
            self.end_FragID = FragIDdata.COMID_orderedFragID[comid][-1]
            self.end_reach = len(FragIDdata.COMID_orderedFragID[comid])
            self.dist = 0
            self.downstream_dist = 0
            self.upstream_dist = 0
            self.interp = False
            self.slope = -9999
            self.upstream = False
            self.downstream = False
            outlet = False

            # If comid is outlet, place "contour" elevation from NHD at end
            # (and upstream interpolation will start from there)
            if COMIDdata.allcomids[comid].to_comid[0] == 0:
                if self.verbose:
                    print "outlet!"
                outlet = True
            try:
                COMIDdata.allcomids[COMIDdata.allcomids[comid].to_comid[0]]
                # additional test to see if COMID listed in FragIDdata
                # if not, is in the list of COMIDs that were flagged earlier as having multiple starts/ends
                # or the stream leaves the grid (to_comid ==[999999])
                FragIDdata.COMID_orderedFragID[to_comid[0]]
            except KeyError:
                if self.verbose:
                    print "outlet!"
                outlet = True
            if outlet:
                self.end_elev = FragIDdata.allFragIDs[self.end_FragID].NHDPlus_elev_min
                self.downstream_contour_comid = comid

            # Go downstream to find downstream contour in pair
            if not outlet:
                if self.verbose:
                    print "looking downstream for first contour"
                self.downstream = True
                self.in_current_comid = False
                self.reachlist = FragIDdata.COMID_orderedFragID[to_comid[0]]
                downcomid = to_comid[0]
                while self.downstream:
                    # reachlist is iterated over;
                    # self.start_reach is always 0 (termination is based on finding a contour)
                    # self.end_reach is updated each time a pair of contours is found
                    self.get_dist_slope(downcomid, FragIDdata, COMIDdata)
                    # if self.downstream = False, a contour or an outlet was found
                    if self.downstream:
                        downcomid = COMIDdata.allcomids[downcomid].to_comid[0]
                        self.reachlist = FragIDdata.COMID_orderedFragID[downcomid]

            # Go upstream to find upstream contour in pair; add distance from downstream
            if self.verbose:
                print "looking upstream for second contour"
            self.upstream_dist += self.downstream_dist
            self.upstream = True
            self.interp = True
             # upstream from end of comid
            while self.interp:
                # reachlist is iterated over;
                # self.start_reach is always 0 (termination is based on finding a contour)
                # self.end_reach is updated each time a pair of contours is found
                self.reachlist = FragIDdata.COMID_orderedFragID[comid][self.start_reach:self.end_reach][::-1]
                self.get_dist_slope(comid, FragIDdata, COMIDdata)

                # if no slope yet, get_dist_slope ran to upstream end of comid without finding contour
                # break out of while loop to move to next comid
                if self.slope == -9999:
                    if self.verbose:
                        print "upstream contour not found in this comid, moving to next"
                    break

                # otherwise, currently an even number of contours in current comid (an upstream contour was found);
                # interpolate between before moving to next
                if self.verbose:
                    print "interpolating within comid..."
                # self.start_FragID was updated to current upstream contour by get_dist_slope
                # self.end_FragID was updated to next contour downstream by get_dist_slope
                self.interpolate_elevs(FragIDdata, self.start_FragID, self.end_FragID)

                # reset everything for next contour pair
                # upstream_dist begins with remainder of end_FragID (upstream contour)
                self.end_FragID = self.start_FragID
                self.end_elev = self.start_elev
                try:
                    self.upstream_dist = np.min(FragIDdata.allFragIDs[self.end_FragID].contour_elev_distance)
                except ValueError:
                    self.upstream_dist = 0
                self.end_reach = FragIDdata.COMID_orderedFragID[comid].index(self.start_FragID)
                self.slope = -9999

                # if updated end_reach is 0, a contour intersects the most upstream FragID in current comid
                # (i.e. entire comid has been interpolated); break out of while loop
                if self.end_reach == 0:
                    self.interp = False
                    break

            # Next upstream contour not in current comid, go to (upstream) from_comid
            if self.interp and self.slope == -9999:
                if self.verbose:
                    print "looking in upstream comids..."
                slopes = []

                # loop through comids in from_comid
                for upcomid in from_comid:
                    if self.verbose:
                        print upcomid
                    self.interp = True  # self.interp might be false from another upcomid in from_comid list

                    # look in subsequent from_comids until upstream contour is found
                    # get_dist_slope will terminate
                    # Note: this does not guarantee that the minimum slope in the tree will be returned
                    # because the loop will terminate upon the first instance of an upstream contour in that branch
                    upcomids = [upcomid]
                    while self.interp:
                        for upid in upcomids:

                            # if multiple upids, and only first is headwater, get_dist_slope will leave interp = True.
                            # Skip past headwater.
                            try:
                                self.reachlist = FragIDdata.COMID_orderedFragID[upid][::-1]
                            except KeyError:
                                # if on the last upcomid, set interp to False and slope to 0
                                if upid == upcomids[-1]:
                                    self.interp = False
                                    self.slope = 0
                                continue

                            self.get_dist_slope(upid, FragIDdata, COMIDdata)
                            if not self.interp:
                                break
                            # need something here that handles missing keys for
                            # streams that originate outside of grid
                            upcomids = COMIDdata.allcomids[upid].from_comid
                            if self.verbose:
                                print "looking further upstream in %s" %(upcomids)

                    # append slope to list of slopes
                    slopes.append(self.slope)

                # slope below confluence is minimum of all slopes going through confluence
                # (avoids possible downstream elevation increase at confluence)
                self.slope = np.min(slopes)

                # to limit interpolation to current comid, reset start_FragID to first reach in current comid
                # end_FragID is preserved from downstream elevation contour in pair
                self.start_FragID = FragIDdata.COMID_orderedFragID[comid][0]

                self.interpolate_elevs(FragIDdata, self.start_FragID, self.end_FragID)

                # add "contour" elevation at first FragID downstream of confluence,
                # so that any subsequent interpolation from upstream will stop here
                # (avoids problems with multiple interpolations through confluence)

                FragIDdata.allFragIDs[self.start_FragID].contour_elev = [FragIDdata.allFragIDs[self.start_FragID].interpolated_contour_elev_max]
                FragIDdata.allFragIDs[self.start_FragID].contour_elev_distance = [0]
                self.elevs_edited.append(self.start_FragID)

            # assign mean elevation and slope for each FragID
            for FragID in FragIDdata.COMID_orderedFragID[comid]:
                FragIDdata.allFragIDs[FragID].interpolated_contour_elev_mean = 0.5 * \
                (FragIDdata.allFragIDs[FragID].interpolated_contour_elev_max + FragIDdata.allFragIDs[FragID].interpolated_contour_elev_min)
                FragIDdata.allFragIDs[FragID].interpolated_contour_slope = \
                (FragIDdata.allFragIDs[FragID].interpolated_contour_elev_max - FragIDdata.allFragIDs[FragID].interpolated_contour_elev_min) / \
                FragIDdata.allFragIDs[FragID].lengthft
        print "\n"


class ElevsFromDEM:

    def __init__(self):

        self.MFgrid_elev_field = 'MEAN'
        self.DEM_elevs_by_cellnum = dict()
        self.ElevsbyFragID = defaultdict(list)
        self.up_increment = 1.0
        self.end_interp_report = "end_interp_report.txt"
        self.adjusted_elev = -9999
        self.verbose = False # prints interpolation information to screen for debugging/troubleshooting


    def DEM_elevs_by_cellnum(self, SFRData, SFROperations):

        # makes dictionary of DEM elevations by cellnum, using table for grid shapefile (column: cellnum)
        # i.e. DEM elevations are at the grid scale
        # (one elevation per cell; multiple FragIDs in that cell would have a uniform elevation)
        arcpy.MakeFeatureLayer_management(SFRData.MFgrid, "MFgrid")

        # if DEM elevations not in grid shapefile,
        # sample DEM with grid shapefile using ZonalStatisticsAsTable function in Arc
        # Note: as of 1/22/2014, this function requires a uniform grid
        Fields = [f.name.lower() for f in arcpy.ListFields("MFgrid")]
        if self.MFgrid_elev_field.lower() not in Fields:
            SFR_arcpy.compute_zonal(self.NROW, self.NCOL, self.DX, SFRData.DEM_z_conversion, SFRData.MFgrid, SFRData.DEM)

        # reference fields for cellnum and elevation, case-insensitive
        SFROperations.getfield("MFgrid", SFRData.node_attribute, SFRData.node_attribute)
        SFROperations.getfield("MFgrid", self.MFgrid_elev_field, self.MFgrid_elev_field.lower())

        # get elevations in each stream cell, from grid shapefile table
        MFgridtable = arcpy.SearchCursor(self.MFgrid)
        self.DEM_elevs_by_cellnum = dict()
        for row in MFgridtable:
            cellnum = row.getValue(SFROperations.joinnames[SFRData.node_attribute])
            self.DEM_elevs_by_cellnum[cellnum] = row.getValue(SFROperations.joinnames[self.MFgrid_elev_field.lower()])


    def DEM_elevs_by_FragID(self, SFRdata, SFROperations):

        # assigns max and min DEM elevations to each FragID following intersect_DEM preprocessing
        try:
            Frag_vertices = arcpy.SearchCursor(SFRdata.intersect_points)
            arcpy.MakeFeatureLayer_management(SFRdata.intersect_points, "intersect_points")
            SFROperations.getfield("intersect_points", "fragid", "fragid")
            SFROperations.getfield("intersect_points", "dem_elev", "dem_elev")
        except:
            raise IOError("Cannot open Fragment end vertices shapefile {0}. "
                          "intersect_DEM in SFRpreproc must be run first.".format(SFRdata.intersect_points))

        # append all DEM elevations sampled along each fragment to a dictionary
        for point in Frag_vertices:
            FragID = point.getValue(SFROperations.joinnames["fragid"])
            elev = point.getValue(SFROperations.joinnames["dem_elev"]) * SFRdata.DEM_z_conversion
            self.ElevsbyFragID[FragID].append(elev)


    def interpolate_between_minima(self, FragIDdata, reachlist, start_reach, end_reach, dS):

        # calculate distance
        dist_to_end = 0
        for FragID in reachlist[start_reach: end_reach]:
            dist_to_end += FragIDdata.allFragIDs[FragID].lengthft
        slope = dS / dist_to_end

        # calculate elevs
        for i in range(len(reachlist))[start_reach: end_reach]:
            #if i > start_reach:
            if reachlist[i] != reachlist[0]:
                FragIDdata.allFragIDs[reachlist[i]].smoothed_DEM_elev_max = \
                FragIDdata.allFragIDs[reachlist[i-1]].smoothed_DEM_elev_min

            #if i+1 < len(reachlist): # don't calculate minimum for end; already assigned
            FragIDdata.allFragIDs[reachlist[i]].smoothed_DEM_elev_min = \
            FragIDdata.allFragIDs[reachlist[i]].smoothed_DEM_elev_max - \
            FragIDdata.allFragIDs[reachlist[i]].lengthft * slope

            FragIDdata.allFragIDs[reachlist[i]].slope = slope

            self.ofp.write('{0}\n'.format(FragIDdata.allFragIDs[reachlist[i]].smoothed_DEM_elev_min))
            if self.verbose:
                print "Reach {0}, max_elev: {1}, min_elev: {2}".format(i,
            FragIDdata.allFragIDs[reachlist[i]].smoothed_DEM_elev_max, FragIDdata.allFragIDs[reachlist[i]].smoothed_DEM_elev_min)


    def assign_elevation_and_slope(self, FragIDs, reach, FragIDdata):
        FragIDdata.allFragIDs[FragIDs[reach-1]].smoothed_DEM_elev_min = self.adjusted_elev
        FragIDdata.allFragIDs[FragIDs[reach-1]].slope = \
            (FragIDdata.allFragIDs[FragIDs[reach-1]].smoothed_DEM_elev_max - self.adjusted_elev) / \
            FragIDdata.allFragIDs[FragIDs[reach-1]].lengthft

    def assign_mean_elevation(self, comid, FragIDdata):
        # assign mean elevation and slope for each FragID
        for FragID in FragIDdata.COMID_orderedFragID[comid]:
            FragIDdata.allFragIDs[FragID].smoothed_DEM_elev_mean = 0.5 * \
            (FragIDdata.allFragIDs[FragID].smoothed_DEM_elev_max + FragIDdata.allFragIDs[FragID].smoothed_DEM_elev_min)
            FragIDdata.allFragIDs[FragID].smoothed_DEM_slope = \
            (FragIDdata.allFragIDs[FragID].smoothed_DEM_elev_max - FragIDdata.allFragIDs[FragID].smoothed_DEM_elev_min) /\
            FragIDdata.allFragIDs[FragID].lengthft

    def connect_downhill(self, FragIDdata):

        print "smoothing DEM elevations along streams... "
        self.ofp = open(self.end_interp_report, 'w')

        knt = 0
        for comid in FragIDdata.COMID_orderedFragID.keys():
            knt += 1
            print "\rCOMID: {0} ({1} of {2})".format(comid, knt, len(FragIDdata.COMID_orderedFragID.keys())),

            FragIDs = FragIDdata.COMID_orderedFragID[comid]
            self.ind_current_minimum = 0

            # set segment end elevations to NHD:
            FragIDdata.allFragIDs[FragIDs[0]].smoothed_DEM_elev_max = FragIDdata.allFragIDs[FragIDs[0]].NHDPlus_elev_max
            FragIDdata.allFragIDs[FragIDs[-1]].smoothed_DEM_elev_min = FragIDdata.allFragIDs[FragIDs[-1]].NHDPlus_elev_min

            if self.verbose:
                print "start_elev: {0}, end_elev {1}".format(FragIDdata.allFragIDs[FragIDs[0]].smoothed_DEM_elev_max,
                                                         FragIDdata.allFragIDs[FragIDs[-1]].smoothed_DEM_elev_min)
            # for segments with only end reaches
            if len(FragIDs) < 2: # no interior points; keep ending elevations and assign mean
                self.assign_mean_elevation(comid, FragIDdata)
                continue
            end = FragIDdata.allFragIDs[FragIDs[-1]].smoothed_DEM_elev_min

            for i in range(len(FragIDs))[1:]:

                if FragIDs[i] == 1926:
                    j=2

                # if COMID start and end elevations are equal, assign all interior elevations to same value
                if FragIDdata.allFragIDs[FragIDs[0]].smoothed_DEM_elev_max == end:
                    self.adjusted_elev = end
                    self.assign_elevation_and_slope(FragIDs, i, FragIDdata)
                    dS = 0
                    self.interpolate_between_minima(FragIDdata, FragIDs, i, len(FragIDs), dS)
                    break

                # assign DEM elevations at Fragment ends using intersections with adjacent fragment pairs
                #FragIDdata.allFragIDs[FragIDs[i-1]].DEM_elev_min = \
                #list(set(self.ElevsbyFragID[FragIDs[i-1]]).intersection(set(self.ElevsbyFragID[FragIDs[i]])))[0]

                # assign DEM elevation at previous Fragment end (same as current Fragment start),
                # using minimum of sampled DEM elevations
                FragIDdata.allFragIDs[FragIDs[i-1]].DEM_elev_min = np.min(self.ElevsbyFragID[FragIDs[i-1]])

                FragIDdata.allFragIDs[FragIDs[i]].DEM_elev_max = FragIDdata.allFragIDs[FragIDs[i-1]].DEM_elev_min

                # adjusted elevation initially set at DEM elevation for beginning of current reach
                self.adjusted_elev = FragIDdata.allFragIDs[FragIDs[i-1]].DEM_elev_min
                if i == 37:
                    j=2

                # in case DEM elevation for current reach is below segment end, set elev. to increment above, interpolate to end
                if self.adjusted_elev <= end:

                    # make sure that initially adjusted elevation is still below previous minimum
                    if end + self.up_increment < FragIDdata.allFragIDs[FragIDs[self.ind_current_minimum]].smoothed_DEM_elev_max:
                        self.adjusted_elev = end + self.up_increment

                    # otherwise set it to the average of the upstream elevation and the end
                    else:
                        self.adjusted_elev = 0.5 * (FragIDdata.allFragIDs[FragIDs[self.ind_current_minimum]].smoothed_DEM_elev_max + end)

                    self.ofp.write("Streambed top at COMID {0} reach {1} is {2}, segment end elevation is {3}\n"
                                   "Readjusting to {4}\n".format(comid, i, FragIDdata.allFragIDs[FragIDs[i-1]].DEM_elev_min,
                                                                 end, self.adjusted_elev))
                    self.ofp.write("interpolating to segment end elevation...\n")

                    # if current reach is past previous minimum, interpolate to current before interpolating to end
                    if i > self.ind_current_minimum:
                        dS = FragIDdata.allFragIDs[FragIDs[self.ind_current_minimum]].smoothed_DEM_elev_max - self.adjusted_elev
                        self.interpolate_between_minima(FragIDdata, FragIDs, self.ind_current_minimum, i, dS)

                    # first assign adjusted elevation and slope, then interpolate to end
                    # this might be redundant
                    self.assign_elevation_and_slope(FragIDs, i, FragIDdata)

                    dS = self.adjusted_elev - end
                    self.interpolate_between_minima(FragIDdata, FragIDs, i, len(FragIDs), dS)
                    break

                # DEM elevation for current reach is above end and below all previous
                elif self.adjusted_elev < FragIDdata.allFragIDs[FragIDs[self.ind_current_minimum]].smoothed_DEM_elev_max:

                    # if previous reach was minimum, append the current DEM elevation
                    if self.ind_current_minimum == i-1:

                        # set minimum elevation in previous reach and assign slope
                        self.assign_elevation_and_slope(FragIDs, i, FragIDdata)

                    # if previous reach was above or equal to minimum, interpolate back to current minimum
                    else:
                        dS = FragIDdata.allFragIDs[FragIDs[self.ind_current_minimum]].smoothed_DEM_elev_max - self.adjusted_elev
                        self.interpolate_between_minima(FragIDdata, FragIDs, self.ind_current_minimum, i, dS)

                    # update max elev to minimum in adjacent upstream reach (interp routine stops after upstream reach)
                    FragIDdata.allFragIDs[FragIDs[i]].smoothed_DEM_elev_max = \
                    FragIDdata.allFragIDs[FragIDs[i-1]].smoothed_DEM_elev_min

                    # update minimum to current reach
                    self.ind_current_minimum = i

                # on the last reach, and DEM elevation for current reach (max DEM elevation) is above end, but >= previous minimum
                # interpolate back to current minimum
                elif i == len(FragIDs)-1:
                    dS = FragIDdata.allFragIDs[FragIDs[self.ind_current_minimum]].smoothed_DEM_elev_max - end
                    self.interpolate_between_minima(FragIDdata, FragIDs, self.ind_current_minimum, i, dS)

                    # update max elev to minimum in adjacent upstream reach (interp routine stops after upstream reach)
                    FragIDdata.allFragIDs[FragIDs[i]].smoothed_DEM_elev_max = \
                    FragIDdata.allFragIDs[FragIDs[i-1]].smoothed_DEM_elev_min

                # DEM elevation for current reach is not below all previous
                else:
                    continue

                if self.verbose and FragIDdata.allFragIDs[FragIDs[i]].smoothed_DEM_elev_max:
                    print "Reach {0}, max_elev: {1}, min_elev: {2}".format(i-1, FragIDdata.allFragIDs[FragIDs[i-1]].smoothed_DEM_elev_max,
                                                                       FragIDdata.allFragIDs[FragIDs[i-1]].smoothed_DEM_elev_min)
            # use min and max smoothed elevations to assign means
            self.assign_mean_elevation(comid, FragIDdata)

        print "\n"
        self.ofp.close()


class SFRoutput:
    """
    Class with methods to write the output both for external processors (e.g. GWV)
    and also to write the SFR file
    """
    def __init__(self, SFRdata):
        self.indat = SFRdata

    def write_SFR_tables(self, SFRSegsAll):
        #open GWV files
        indat = self.indat
        mat1out = open(indat.MAT1, 'w')
        mat1out.write('row,column,layer,stage,top_streambed,reach,segment,'
                      'width_in_cell,length_in_cell,bed_K,bed_thickness,bed_slope,bed_roughness\n')

        mat2out = open(indat.MAT2, 'w')
        mat2out.write('segment,icalc,outseg,iupseg,iprior,nstrpts,flow,runoff,'
                      'etsw,pptsw,roughch,roughbk,cdepth,fdepth,awdth,bwdth\n')

        seglist = sorted(list(SFRSegsAll.allSegs.keys()))
        nss = len(seglist)  # calcualte the total number of segments

        for cseg in seglist:
            reachlist = sorted(SFRSegsAll.allSegs[cseg].seg_reaches)
            curr_reaches = SFRSegsAll.allSegs[cseg].seg_reaches
            #
            # write to the mat1out file with information on each reach
            #
            for creach in reachlist:
                printstring = (curr_reaches[creach].row,
                curr_reaches[creach].column, 1)
                mat1out.write(','.join(map(str, printstring)))
                mixedlist = (curr_reaches[creach].elevreach,
                    curr_reaches[creach].elevreach - indat.stream_depth,
                    creach,
                    cseg,
                    curr_reaches[creach].eff_width,
                    curr_reaches[creach].eff_length,
                    indat.bedK,
                    indat.bedthick,
                    curr_reaches[creach].eff_slope,
                    indat.roughch)
                printstring = ',{0:.4e},{1:.4e},{2:d},{3:d},{4:.4e},{5:.4e},{6:.4e},{7:.4e},{8:.4e},{9:.4e}\n'.format(
                    *mixedlist)
                mat1out.write(printstring)
            #
            # now write to the mat2out file with information on each segment
            #
            curr_seg = SFRSegsAll.allSegs[cseg]
            iupseg = 0
            iprior = 0
            flow = 0
            #write ouput to GWV matrix 2 file
            mlist1 = (cseg, curr_seg.icalc, curr_seg.outseg, iupseg, iprior, indat.nstrpts)
            mlist2 = (flow, curr_seg.runoff, curr_seg.etsw, curr_seg.pptsw, curr_seg.roughch, indat.roughbk)
            mlist3 = (indat.cdepth, indat.fdepth, indat.awdth, indat.bwdth)
            printstring = '{0:d},{1:d},{2:d},{3:d},{4:d},{5:d}'.format(*mlist1)
            printstring = printstring+',{0:.5e},{1:.5e},{2:.5e},{3:.5e},{4:.5e},{5:.5e}'.format(*mlist2)
            printstring = printstring+',{0:.5e},{1:.5e},{2:.5e},{3:.5e}'.format(*mlist3)
            mat2out.write('{0:s}\n' .format(printstring))

        mat1out.close()
        mat2out.close()

    def build_SFR_package(self):

        print "Building new SFR package file {0}...".format(self.indat.OUT)
        Mat1 = np.genfromtxt(self.indat.MAT1, delimiter=',', names=True, dtype=None)
        Mat2 = np.genfromtxt(self.indat.MAT2, names=True, delimiter=',', dtype=None)

        SFRoutfile = self.indat.OUT

        if self.indat.tpl:
            SFRoutfile = '_'.join(self.indat.OUT.split('.'))+'.tpl'
        else:
            #  MNF note ---> commenting this out for now to allow for variable conductance values
            '''
            # if building an SFR package (after PEST has written it)
            # read in PEST-adjusted value for streambed conductance, SFRoutfile must exist
            if os.path.isfile(SFRoutfile):
                SFRinputdata = open(SFRoutfile, 'r').readlines()
                if len(SFRinputdata) > 0:  #in case the file is empty...
                    Cond = float(SFRinputdata[3].strip().split()[-1])
                    Mat1['bed_K'] = np.ones(len(Mat1))*Cond
            '''
        nreaches = len(Mat1)
        nseg = int(np.max(Mat1['segment']))
        ofp = open(SFRoutfile, 'w')
        if self.indat.tpl:
            ofp.write("ptf ~\n")
        ofp.write("#SFRpackage file generated by Python Script using NHDPlus v2 data\n")
        ofp.write('{0:d} {1:d} {2:d} {3:d} {4:e} {5:e} {6:d} {7:d} {8:d} {9:d} {10:d} {11:d}\n'.format(
            -1*nreaches,
            nseg,
            self.indat.nsfrpar,
            self.indat.nparseg,
            self.indat.const,
            self.indat.dleak,
            self.indat.istcb1,
            self.indat.istcb2,
            self.indat.isfropt,
            self.indat.nstrail,
            self.indat.isuzn,
            self.indat.nsfrsets
        ))

        for i in range(len(Mat1)):
            slope = Mat1['bed_slope'][i]
            if slope <= self.indat.minimum_slope:  # one last check for negative or zero slopes
                slope = self.indat.minimum_slope
            if self.indat.tpl:
                bedK = '~SFRc~'
            else:
                bedK = '{0:e}'.format(Mat1['bed_K'][i])
            ofp.write('{0:d} {1:d} {2:d} {3:d} {4:d} {5:e} {6:e} {7:e} {8:e} {9:s}\n'.format(
                Mat1['layer'][i],
                Mat1['row'][i],
                Mat1['column'][i],
                int(Mat1['segment'][i]),
                int(Mat1['reach'][i]),
                Mat1['length_in_cell'][i],
                Mat1['top_streambed'][i],
                slope,
                self.indat.bedthick,
                bedK))

        ofp.write('{0:d} 0 0 0\n'.format(nseg))
        for i in range(len(Mat2)):
            seg = Mat2['segment'][i]
            seg_Mat1_inds = np.where(Mat1['segment'] == seg)
            seginfo = Mat1[seg_Mat1_inds]
            ofp.write('{0:d} {1:d} {2:d} 0 {3:e} 0.0 0.0 0.0 {4:e}\n'.format(
                seg,
                self.indat.icalc,
                Mat2['outseg'][i],
                Mat2['flow'][i],
                self.indat.roughch
                ))

            #need to add the input of modify_segment_widths
            #if self.indat.modify_segment_widths and seg in segments2modify:
            if seginfo['width_in_cell'][0] > 0.:
                ofp.write('{0:e}\n'.format(seginfo['width_in_cell'][0]))
                ofp.write('{0:e}\n'.format(seginfo['width_in_cell'][-1]))
            else:
                if self.indat.icalc == 0:
                    ofp.write('{0:e} {1:e}\n'.format(seginfo['width_in_cell'][0], self.indat.stream_depth))
                    ofp.write('{0:e} {1:e}\n'.format(seginfo['width_in_cell'][-1], self.indat.stream_depth))
                else:
                    ofp.write('{0:e}\n'.format(seginfo['width_in_cell'][0]))
                    ofp.write('{0:e}\n'.format(seginfo['width_in_cell'][-1]))
        ofp.close()

    def build_SFR_shapefile(self, SFRSegsAll):
        print 'building a shapefile of final SFR segments and reaches'
        Mat1 = np.genfromtxt(self.indat.MAT1, delimiter=',', names=True, dtype=None)
        #path = os.getcwd()
        #arcpy.env.workspace = path
        arcpy.env.workspace = self.working_dir
        outshapefile = self.indat.GISSHP
        rivcells = self.indat.CELLS_DISS
        if arcpy.Exists(outshapefile):
            arcpy.Delete_management(outshapefile)

        #add fields including outseg for checking, note could have multiple
        #lines for the same cellnum.
        arcpy.CreateFeatureclass_management(self.working_dir, outshapefile, "POLYGON", rivcells)
        arcpy.AddField_management(outshapefile, self.indat.node_attribute, "LONG")
        arcpy.AddField_management(outshapefile, "ROW", "LONG")
        arcpy.AddField_management(outshapefile, "COLUMN", "LONG")
        arcpy.AddField_management(outshapefile, "LAYER", "LONG")
        arcpy.AddField_management(outshapefile, "SEGMENT", "LONG")
        arcpy.AddField_management(outshapefile, "REACH", "LONG")
        arcpy.AddField_management(outshapefile, "OUTSEG", "LONG")
        arcpy.AddField_management(outshapefile, "LEVELPATH", "LONG")
        newrows = arcpy.InsertCursor(outshapefile)
        shapeName = arcpy.Describe(rivcells).shapeFieldName

        for seg in SFRSegsAll.allSegs.iterkeys():
            reachdict = SFRSegsAll.allSegs[seg].seg_reaches
            for rch in reachdict.iterkeys():
                localcell = reachdict[rch].cellnum

                # get layer from Mat 1
                Mat1_ind = np.where((Mat1['segment'] == seg) & (Mat1['reach'] == rch))[0][0]
                layer = Mat1['layer'][Mat1_ind]

                query = "CELLNUM={0}".format(localcell)
                poly = arcpy.SearchCursor(rivcells, query)
                # skips if poly is empty, should only go through it once
                for entry in poly:
                    newvals = newrows.newRow()
                    feature = entry.getValue(shapeName)
                    newvals.Shape = feature
                    newvals.CELLNUM = localcell
                    newvals.ROW = reachdict[rch].row
                    newvals.COLUMN = reachdict[rch].column
                    newvals.LAYER = layer
                    newvals.SEGMENT = seg
                    newvals.REACH = rch
                    newvals.OUTSEG = SFRSegsAll.allSegs[seg].outseg
                    newvals.LEVELPATH = SFRSegsAll.allSegs[seg].levelpathID
                    newrows.insertRow(newvals)

        del newvals, newrows

    def buildSFRshapefile2(self):

        try:
            from fiona import collection
            from shapely.geometry import shape, mapping
            import pandas as pd
        except:
            print "Need these modules for class buildSFRshapefile2: pandas\nfiona\nshapely\n"
            quit()

        print "reading {0} and {1} into pandas DataFrames...".format(self.indat.MAT1, self.indat.MAT2)
        Mat1 = pd.read_csv(self.indat.MAT1)

        # calculate column of cellnums for Mat1; index by cellnum
        Mat1[self.indat.node_attribute] = Mat1.apply(lambda x: (x['row'] - 1) * self.indat.NCOL + x['column'], axis=1).astype('int32')
        #Mat1[self.indat.node_attribute] = Mat1[self.indat.node_attribute]
        Mat1 = Mat1.set_index([Mat1[self.indat.node_attribute], Mat1.index])
        Mat2 = pd.read_csv(self.indat.MAT2)
        Mat2.index = Mat2['segment']
        Mat2upsegs = pd.Series(Mat2['segment'].values, index=Mat2['outseg']).to_dict()
        print "making shapefile of SFR network using {}...".format(self.indat.CELLS_DISS)
        knt = 0
        nSFRcells = len(Mat1)
        with collection(self.indat.CELLS_DISS, "r") as input:
            # schema = input.schema.copy()
            schema = {'geometry': 'Polygon',
                      'properties': {self.indat.node_attribute: 'int',
                                     'row': 'int',
                                     'column': 'int',
                                     'layer': 'int',
                                     'segment': 'int',
                                     'reach': 'int',
                                     'outseg': 'int',
                                     'upseg': 'int',
                                     'sb_elev': 'float',
                                     'modeltop': 'float'}}

            with collection(self.indat.GISSHP, "w", "ESRI Shapefile", schema) as output:
                for node in input:
                    cellnum = node['properties'][self.indat.node_attribute]
                    if cellnum == 381937:
                        j=2
                    # handle cells that were subsequently deleted after creation of input shapefile
                    try:
                        Mat1.ix[cellnum]

                    except KeyError:
                        continue

                    print "\r{:d}%".format(100 * knt / nSFRcells),
                    # iterate over multi-index (allows for multiple reaches in one cell)
                    for uniquereach in Mat1.ix[cellnum].index:
                        knt += 1
                        segment = int(Mat1.ix[(cellnum, uniquereach), 'segment'])
                        outseg = int(Mat2.ix[segment, 'outseg'])
                        # handle headwaters (no upseg)
                        try:
                            Mat2upsegs[segment]
                        except KeyError:
                            Mat2upsegs[segment] = 0

                        # shapefiles are incompatible with int64.
                        # but apparently they are compatible with float64 ARGH!
                        output.write({'properties': {
                                         self.indat.node_attribute: Mat1.ix[(cellnum, uniquereach), self.indat.node_attribute],
                                         'row': Mat1.ix[(cellnum, uniquereach), 'row'].astype('int32'),
                                         'column': Mat1.ix[(cellnum, uniquereach), 'column'].astype('int32'),
                                         'layer': Mat1.ix[(cellnum, uniquereach), 'layer'].astype('int32'),
                                         'segment': int(segment),
                                         'reach': Mat1.ix[(cellnum, uniquereach), 'reach'].astype('int32'),
                                         'outseg': outseg,
                                         'upseg': np.int32(Mat2upsegs[segment]),
                                         'sb_elev': Mat1.ix[(cellnum, uniquereach), 'top_streambed'],
                                         'modeltop': 99.99},
                                      'geometry': mapping(shape(node['geometry']))})

        # copy over prj file
        shutil.copyfile("{}.prj".format(self.indat.CELLS_DISS[:-4]), "{}.prj".format(self.indat.GISSHP[:-4]))



def widthcorrelation(arbolate):
    #estimate widths, equation from Feinstein and others (Lake
    #Michigan Basin model) width=0.1193*(x^0.5032)
    # x=arbolate sum of stream upstream of the COMID in meters
    #NHDPlus has arbolate sum in kilometers.
    #print a table with reachcode, order, estimated width, Fcode
    estwidth = 0.1193*math.pow(1000*arbolate, 0.5032)
    return estwidth

"""
###################
ERROR EXCEPTION CLASSES
###################
"""
class InputFileMissing(Exception):
    def __init__(self, infile):
        self.infile = infile
    def __str__(self):
        return('\n\nCould not open or parse input file {0}.\nCheck for errors in XML formatting.'.format(self.infile))

class NoPreProcessing(Exception):
    def __init__(self, infile):
        self.infile = infile
    def __str__(self):
        return('\n\nCannot find intermediate (output) files from pre-processing routine. '
               'Make sure that <preproc> is set to True in input XML file and try again.')

class BadElevChoice(Exception):
    def __init__(self, choice):
        self.choice = choice
    def __str__(self):
        return('\n\nElevation Choice in XML is not allowable {0}\n'.format(self.choice))

class MissingFragID(Exception):
    def __init__(self, container, FragID, COMID):
        self.container = container
        self.FragID = FragID
        self.COMID = COMID

    def __str__(self):
        return('\n\nFragID {0} not in {1}. This FragID belongs to COMID {2}. '
               'Check NHD input; this COMID may need to be deleted or edited'.format(self.FragID, self.container, self.COMID))
