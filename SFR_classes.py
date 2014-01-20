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

        self.compute_zonal = self.tf2flag(inpars.findall('.//compute_zonal')[0].text)
        self.preproc = self.tf2flag(inpars.findall('.//preproc')[0].text)
        self.reach_cutoff = float(inpars.findall('.//reach_cutoff')[0].text)
        self.rfact = float(inpars.findall('.//rfact')[0].text)
        self.MFgrid = inpars.findall('.//MFgrid')[0].text
        self.MFdomain = inpars.findall('.//MFdomain')[0].text
        self.MFdis = inpars.findall('.//MFdis')[0].text
        self.DEM = inpars.findall('.//DEM')[0].text
        self.intersect = inpars.findall('.//intersect')[0].text
        self.rivers_table = inpars.findall('.//rivers_table')[0].text
        self.PlusflowVAA = inpars.findall('.//PlusflowVAA')[0].text
        self.Elevslope = inpars.findall('.//Elevslope')[0].text
        self.Flowlines_unclipped = inpars.findall('.//Flowlines_unclipped')[0].text
        self.arcpy_path = inpars.findall('.//arcpy_path')[0].text
        self.FLOW = inpars.findall('.//FLOW')[0].text
        self.FTab = inpars.findall('.//FTab')[0].text
        self.Flowlines = inpars.findall('.//Flowlines')[0].text
        self.ELEV = inpars.findall('.//ELEV')[0].text
        self.CELLS = inpars.findall('.//CELLS')[0].text
        self.CELLS_DISS = inpars.findall('.//CELLS_DISS')[0].text
        self.NHD = inpars.findall('.//NHD')[0].text
        self.OUT = inpars.findall('.//OUT')[0].text
        self.MAT1 = inpars.findall('.//MAT1')[0].text
        self.MAT2 = inpars.findall('.//MAT2')[0].text
        self.WIDTH = inpars.findall('.//WIDTH')[0].text
        self.MULT = inpars.findall('.//MULT')[0].text
        self.ELEVcontours = inpars.findall('.//ELEVcontours')[0].text
        self.Routes = inpars.findall('.//Routes')[0].text
        self.Contours_intersect = inpars.findall('.//Contours_intersect')[0].text
        self.Contours_intersect_distances = inpars.findall('.//Contours_intersect_distances')[0].text
        self.RCH = inpars.findall('.//RCH')[0].text

        try:
            self.eps = float(inpars.findall('.//eps')[0].text)
        except:
            self.eps = 1.0000001e-02  # default value used if not in the input file

        # conversion for vertical length units between NHDPlus v2. (usually cm) and model
        try:
            self.z_conversion = float(inpars.findall('.//z_conversion')[0].text)
        except:
            self.z_conversion = 1.0/(2.54 *12)  # default value used if not in the input file

        #cutoff to check stream length in cell against fraction of cell dimension
        #if the stream length is less than cutoff*side length, the piece of stream is dropped
        try:
            self.cutoff = float(inpars.finall('.//cutoff')[0].text)
        except:
            self.cutoff = 0.0

        #read the Fcode-Fstring table and save it into a dictionary, Fstring
        descrips = arcpy.SearchCursor(self.FTab)
        self.Fstring = dict()
        for description in descrips:
            Fcodevalue = int(description.FCode)
            if not Fcodevalue in self.Fstring:
                self.Fstring[Fcodevalue]=description.Descriptio
        del descrips

        # initialize the arcpy environment
        arcpy.env.workspace = os.getcwd()
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


class FragIDprops(object):
    """
    Properties for each COMID
    """
    '''
    __slots__ = ['comid', 'startx', 'starty', 'endx', 'endy', 'FragID',
                 'maxsmoothelev', 'minsmoothelev', 'lengthft',
                 'cellnum', 'contour_elev','elev', 'sidelength',
                 'segelevinfo', 'start_has_end', 'end_has_start', 'elev_distance', 'segelevinfo',
                 'elev_min','elev_max','elev_mean',]
    '''
    #  using __slots__ makes it required to declare properties of the object here in place
    #  and saves significant memory
    def __init__(self, comid, startx, starty, endx, endy, FragID,
                 maxsmoothelev, minsmoothelev, lengthft, cellnum,contour_elev,segelevinfo,elev_min,elev_max,elev_mean,elev_distance):
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
        self.segelevinfo = segelevinfo
        self.elev_min = elev_min
        self.elev_max = elev_max
        self.elev_mean = elev_mean
        self.elev_distance = elev_distance


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
                if reachlength < CELLdata.allcells[cellnum].sidelength*SFRdata.cutoff:
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
    __slots__ = ['delx', 'dely', 'sidelength', 'row', 'column']
    '''
    def __init__(self, delx, dely, sidelength, row, column):
        self.delx = delx
        self.dely = dely
        self.sidelength = sidelength
        self.row = row
        self.column = column

class CellPropsAll:
    def __init__(self):
        self.allcells = dict()

    def populate_cells(self, SFRdata):
        #use the CELLS shapefile to get the length of the cell sides, used
        #to weed out short river reaches.  If the length in the cell is
        #less than sidelength * the input 'cutoff', also get row and column
        #for each cellnumber

        cells = arcpy.SearchCursor(SFRdata.CELLS)

        for cell in cells:
            cellnum = int(cell.CELLNUM)
            dx = float(cell.delx)
            dy = float(cell.dely)
            minside = float(cell.delx)
            if float(cell.dely) < minside:
                minside = float(cell.dely)
            row = int(cell.row)
            column = int(cell.column)
            self.allcells[cellnum] = CellProps(dx, dy, minside, row, column)


class COMIDPropsAll:
    def __init__(self):
        self.allcomids = dict()
        self.hydrosequence_sorted = list()
        self.hydrosequence_comids = dict()

    def populate_routing(self, SFRdata, FragIDdata, LevelPathdata):
        """
        Read the COMID routing information from the SFRdata.FLOW file
        """
        print ('Reading in routing information from {0:s}'.format(SFRdata.FLOW))
        # open the SFRdata.FLOW file as read-only (using SearchCursor)

        CLIP = np.loadtxt('boundaryClipsRouting.txt', skiprows=1, delimiter=',', dtype=int)

        for ccomid in FragIDdata.allcomids:
            self.allcomids[ccomid] = COMIDprops()

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
        comidseen = list()
        with arcpy.da.SearchCursor(SFRdata.PlusflowVAA,
                                   ("ComID", "Hydroseq", "uphydroseq", "dnhydroseq",
                                   "ReachCode","StreamOrde","ArbolateSu","Fcode", "levelpathI")) as cursor:
            for crow in cursor:
                comid = int(crow[0])
                hydrosequence = int(crow[1])
                uphydrosequence = int(crow[2])
                downhydrosequence = int(crow[3])
                reachcode = crow[4]
                stream_order = int(crow[5])
                arbolate_sum = float(crow[6])
                Fcode = int(crow[7])
                levelpathid = int(crow[8])
                if int(comid) in FragIDdata.allcomids:
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

            # find unique levelpathIDs
            LevelPathdata.level_ordered = sorted(list(set(LevelPathdata.level_ordered)), reverse=True)

            for clevelpathid in LevelPathdata.level_ordered:
                LevelPathdata.allids[clevelpathid] = LevelPathIDprops()
                LevelPathdata.levelpath_FragID[clevelpathid] = []


            # assign levelpathID routing
        del crow, cursor

        with arcpy.da.SearchCursor(SFRdata.PlusflowVAA,
                                   ("ComID", "LevelPathI", "DnLevelPat")) as cursor:
            for crow in cursor:
                comid = int(crow[0])
                levelpathid = int(crow[1])
                downlevelpathid = int(crow[2])
                if levelpathid in LevelPathdata.level_ordered:
                    if downlevelpathid != levelpathid:
                        LevelPathdata.allids[levelpathid].down_levelpathID = downlevelpathid
                    if comid in FragIDdata.comid_FragID:
                        LevelPathdata.levelpath_FragID[levelpathid].extend(FragIDdata.comid_FragID[comid])


        comid_missing = list(set(FragIDdata.allcomids).difference(comidseen))
        if len(comid_missing) > 0:
            print "WARNING! the following COMIDs are missing from \n{0:s}".format('\n'.join(map(str(comid_missing))))

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
                None,None,None,None,None,None)

            self.allcomids.append(int(seg.COMID))
        self.allcomids = list(set(self.allcomids))




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
    __slots__ = ['cellnum','rch_number','eff_length','eff_slope','elevreach','bedthick','bedK','roughch']
    '''
    def __init__(self):
        self.cellnum = None
        self.rch_number = None
        self.eff_length = None
        self.eff_slope = None
        self.elevreach = None
        self.bedthick = None
        self.bedK = None
        self.roughch = None


class SFRSegmentProps(object):
    """
    class object to hold final SFR segment information
    """
    '''
    __slots__ = ['seg_cells','icalc','iupseg','outseg','runoff','etsw','pptsw']
    '''
    def __init__(self):
        self.seg_cells = list()
        self.outseg = None
        self.icalc = None
        self.iupseg = 0    #iupseg default zero right now, diversions could be added
        self.outseg = None
        self.runoff = None
        self.estw = None
        self.pptsw = None
        self.seg_reaches = dict()   #reach object with properties from SFRReachProps keyed by rch_number


class SFRSegmentsAll:
    """
    class that makes up a dictionary of SFR objects
    and contains methods to accumulate properties of the
    same levelpathID within a call, find confluences, subdivide
    levelpathIDs and establish the final SFR segments
    """
    def __init__(self):
        self.allSegs = dict()
        self.confluences = defaultdict(list)    # dictionary of lists keyed by levelpathID

    def divide_at_confluences(self,LevelPathdata, FragIDdata, COMIDdata):
        #establish provisional segment numbers from downstream (decending)
        #list of levelpathIDs
        provSFRseg = dict()
        for i,lpID in enumerate(LevelPathdata.level_ordered):
            provSFRseg[lpID] = i

        for clevelpathid in LevelPathdata.allids.iterkeys():
            iseg = provSFRseg[clevelpathid]
            nextlevelpath = LevelPathdata.allids[clevelpathid].down_levelpathID
            cellist = LevelPathdata.allids[clevelpathid].ordered_cellnums
            if len(cellist) == 0:
                print 'check levelpathID {0}, no cells are assigned'.format(clevelpathid)
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
                    outsegbeg = LevelPathdata.allids[outid].ordered_cellnums[0]
                    confl = -1    # flag in case an end is not found
                    for cell in LevelPathdata.allids[outid].ordered_cellnums:
                        if cell == isegend:
                            confl = cell
                else:
                    #no downstream levelpath
                    outsegbeg = 0
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
                #put it in downstream order with no repeats
                self.confluences[clevelpathid]=[cell for cell in LevelPathdata.allids[clevelpathid].ordered_cellnums \
                    if cell in set(self.confluences[clevelpathid])]


    def accumulate_same_levelpathID(self, LevelPathdata, COMIDdata, FragIDdata):
        """
        method to add lengths and weight widths and slopes
        for parts of a stream that have the same levelpathID
        within a cell
        """
        self.totlength=dict()   #these should be pointed to reach properties....
        self.weightwidth=dict()
        self.elevcell=dict()
        self.weightedslope=dict()
        for lpID in LevelPathdata.level_ordered:
            segment=SFRprovseg[lpID]
            for localcell in yyyyyy:
                tt=0.
                ww=0.
                ws=0.
                el=0.
                for cFragID in xxxxxx:
                    knt=0
                    if lpcell==localcell:
                        knt=knt+1
                        tt = tt + FragIDdata.allFragIDs[cFragID].lengthft
                        ww = ww + COMIDdata[ccomid].est_width* FragIDdata.allFragIDs[cFragID].lengthft
                       # el = el +
                       # ws = ws +
                if knt==0:
                    totlength[cFragID] = FragIDdata.allFragIDs[cFragID].lengthft
                    #elevcell[lpID][localcell] = riverelev[localcell][0]
                    #weightedslope[lpID][localcell] = cellslope[localcell][0]
                    weightwidth[cFragID]= COMIDdata[ccomid].est_width
                else:
                    totlength[cFragID] = tt
                    #elevcell[cFragID]=el/float(knt)
                    if tt>0:
                        weightwidth[cFragID] = ww/tt
                        #weightedslope[cFragID] = ws/tt
                    else:
                        weightwidth[cFragID] = 99999.
                        #weightedslope[cFragID] = 99999.

class SFRpreproc:
    def __init__(self, SFRdata):
        self.joinnames = dict()
        self.indata = SFRdata
        ts = time.time()
        st_time = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S')
        self.ofp = open('SFR_preproc.log', 'w')
        self.ofp.write('SFR_preproc log.')
        self.ofp.write('\n' + '#' * 25 + '\nStart Time: {0:s}\n'.format(st_time) + '#' * 25 + '\n')


    def getfield(self, table, joinname, returnedname):
    # get name of field (useful for case issues and appended field names in joined tables, etc)
        Fields = arcpy.ListFields(table)
        joinname = joinname.lower()
        for field in Fields:
            if joinname in field.name.lower():
                joinname = field.name
                break
        self.joinnames[returnedname] = joinname

    def clip_and_join_attributes(self):
        indat = self.indata

        print "Clip original NHD flowlines to model domain..."
        # this creates a new file so original dataset untouched
        arcpy.Clip_analysis(indat.Flowlines_unclipped,
                            indat.MFdomain,
                            indat.Flowlines)

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
        print "\nkeeping: %s fields; deleting the rest" % ('\n'.join(fields2keep))
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
                    print field.name
            if len(fields2delete) > 0:
                arcpy.DeleteField_management(table, fields2delete)
        # Join PlusflowVAA and Elevslope to Flowlines
        arcpy.MakeFeatureLayer_management(indat.Flowlines, "Flowlines")

        # If not already, permanently join PlusflowVAA and Elevslope, then to Flowlines
        if Join:
            print "\nJoining Elevslope to PlusflowVAA...\n"
            self.getfield("PlusflowVAA", "comid", "comid1")
            self.getfield("Elevslope", "comid", "comid2")
            arcpy.JoinField_management("PlusflowVAA",
                                       self.joinnames['comid1'],
                                       "Elevslope",
                                       self.joinnames['comid2'])
        else:
            print "PlusflowVAA and Elevslope already joined from previous run..."

        print "Joining PlusflowVAA to NHDFlowlines...\n"
        self.getfield("Flowlines", "comid", 'comid1')

        # join to Flowlines, keeping only common
        SFR_arcpy.general_join(indat.Flowlines,
                               "Flowlines",
                               self.joinnames['comid1'],
                               indat.PlusflowVAA,
                               "comid",
                               True)

        # reopen flowlines as "Flowlines" --> clunky a bit to save and reopen, but must do so
        arcpy.MakeFeatureLayer_management(indat.Flowlines, "Flowlines")

        print "\n"
        self.ofp.write('\n' + 25*'#' + '\nRemoving segments with no elevation information, and with ThinnerCod = -9..\n')
        print "Removing segments with no elevation information, and with ThinnerCod = -9..."

        self.getfield("Flowlines", "thinnercod", "ThinnerCod")
        self.getfield("Flowlines", "maxelevsmo", "MaxEl")
        self.getfield("Flowlines", "comid", "comid")
        FLtable = arcpy.UpdateCursor("Flowlines")
        zerocount = 0
        tcount = 0
        for segments in FLtable:
            if segments.getValue(self.joinnames['MaxEl']) == 0:
                print "%d no elevation data" % segments.getValue(self.joinnames['comid'])
                self.ofp.write("%d no elevation data\n" % segments.getValue(self.joinnames['comid']))
                FLtable.deleteRow(segments)
                zerocount += 1
            elif segments.getValue(self.joinnames['ThinnerCod']) == -9:
                print "%d ThinnerCod=-9" % segments.getValue(self.joinnames['comid'])
                self.ofp.write("%d ThinnerCod=-9\n" % segments.getValue(self.joinnames['comid']))
                FLtable.deleteRow(segments)
                tcount += 1
        #  read in discretization information
        DX, DY, NLAY, NROW, NCOL, i = disutil.read_meta_data(indat.MFdis)

        # update the "node" field in indat.MFgrid
        # if there is a field with unique values, assume it's ok
        # otherwise delete the column if it already exists
        # NB --> only evaluating the first 100 rows...
        print 'verifying that there is a "cellnum" field in {0:s}'.format(indat.MFgrid)
        hascellnum = False
        cellnumunique = 0
        MFgridflds = arcpy.ListFields(indat.MFgrid)
        for cfield in MFgridflds:
            if cfield.name.lower() == 'cellnum':
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
                    cellvals.append(row.getValue('cellnum'))
            cellnumunique = len(set(cellvals))
            del cellvals
            del row
            del cursor

        if cellnumunique > 1:
            print '"cellnum" field in place with unique values in {0:s}'.format(indat.MFgrid)
        else:
            for fld in arcpy.ListFields(indat.MFgrid):
                if fld == 'cellnum':
                    arcpy.DeleteField_management(indat.MFgrid, 'cellnum')
            arcpy.AddField_management(indat.MFgrid, 'cellnum', 'LONG')
            calcexpression = '((!row!-1)*{0:d}) + !column!'.format(NCOL)
            arcpy.CalculateField_management(indat.MFgrid, 'cellnum', calcexpression, 'PYTHON')
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
        self.getfield(indat.CELLS, "cellnum", "cellnum")
        arcpy.Dissolve_management(indat.CELLS, indat.CELLS_DISS, self.joinnames['cellnum'])

        print "Exploding NHD segments to grid cells using Intersect and Multipart to Singlepart..."
        arcpy.Intersect_analysis([indat.CELLS_DISS, "Flowlines"], "tmp_intersect.shp")
        arcpy.MultipartToSinglepart_management("tmp_intersect.shp", indat.intersect)
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
        self.getfield(indat.intersect, "comid", "comid")
        self.getfield(indat.intersect, "cellnum", "cellnum")
        self.getfield(indat.intersect, "lengthft", "Length")
        table = arcpy.UpdateCursor(indat.intersect)
        count = 0
        for reaches in table:
            if reaches.getValue(self.joinnames["Length"]) <= indat.reach_cutoff:
                print "segment: %d, cell: %s, length: %s" % (reaches.getValue(self.joinnames["comid"]),
                                                            reaches.getValue(self.joinnames["cellnum"]),
                                                            reaches.getValue(self.joinnames["Length"]))
                self.ofp.write("segment: %d, cell: %s, length: %s\n"
                          %(reaches.getValue(self.joinnames["comid"]),
                            reaches.getValue(self.joinnames["cellnum"]),
                            reaches.getValue(self.joinnames["Length"])))
                table.deleteRow(reaches)
                count += 1
        print "removed %s reaches with lengths <= %s\n" % (count, indat.reach_cutoff)
        self.ofp.write("removed %s reaches with lengths <= %s\n" % (count, indat.reach_cutoff))
        print "removing cells corresponding to those reaches..."

        # temporarily join river_cells_dissolve to river explode; record nodes with no elevation information
        arcpy.MakeFeatureLayer_management(indat.CELLS_DISS, "river_cells_dissolve")
        arcpy.MakeFeatureLayer_management(indat.intersect, "river_explode")
        arcpy.AddJoin_management("river_cells_dissolve",
                                 "cellnum",
                                 "river_explode",
                                 "cellnum",
                                 "KEEP_ALL")  # this might not work as-in in stand-alone mode
        self.getfield("river_cells_dissolve", "cellnum", "cellnum")
        self.getfield("river_cells_dissolve", "maxelevsmo", "maxelevsmo")
        table = arcpy.SearchCursor("river_cells_dissolve")
        nodes2delete = []
        for row in table:
            if row.isNull(self.joinnames["maxelevsmo"]):
                nodes2delete.append(row.getValue(self.joinnames["cellnum"]))
        arcpy.RemoveJoin_management("river_cells_dissolve", "river_explode")

        # preserve the FragID code as FragID for all future reference
        arcpy.AddField_management(indat.intersect, "FragID", "LONG")
        arcpy.CalculateField_management(indat.intersect, "FragID", "!FID!", "PYTHON")

        # remove cellnums with no elevation information from river_explode
        self.ofp.write('\n' + 25*'#' + '\nRemoving nodes with no elevation information from river_explode\n')
        print 'Removing nodes with no elevation information from river_explode'
        self.getfield(indat.CELLS_DISS, "cellnum", "cellnum")
        table = arcpy.UpdateCursor(indat.CELLS_DISS)
        count = 0
        for cells in table:
            if cells.getValue(self.joinnames["cellnum"]) in nodes2delete:
                print "%d" % (cells.getValue(self.joinnames["cellnum"]))
                self.ofp.write('%d\n' % (cells.getValue(self.joinnames["cellnum"])))
                table.deleteRow(cells)
                count += 1
        print "removed %s cells\n" % (count)
        self.ofp.write("removed %s cells" % (count))

        print "removing any remaining disconnected reaches..."
        self.ofp.write('\n' + 25*'#' + '\nremoving any remaining disconnected reaches...\n')
        self.getfield(indat.intersect, "cellnum", "cellnum")
        table = arcpy.UpdateCursor(indat.intersect)
        count = 0
        for reaches in table:
            if reaches.getValue(self.joinnames["cellnum"]) in nodes2delete:
                print "%d" % (reaches.getValue(self.joinnames["cellnum"]))
                self.ofp.write('%d\n' % (reaches.getValue(self.joinnames["cellnum"])))
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

class SFROperations:
    """
    class to make operations on data
    """
    def __init__(self, SFRdata):
        self.SFRdata = SFRdata
        self.newCOMIDdata = dict()

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
        # bring data in as layers
        #set workspace


        #convert to layers
        arcpy.MakeFeatureLayer_management(self.SFRdata.Flowlines_unclipped, 'nhd_lyr')
        arcpy.MakeFeatureLayer_management(self.SFRdata.MFdomain, 'grid_lyr')
        print "selecting streams that cross model grid boundary..."
        arcpy.SelectLayerByLocation_management('nhd_lyr', 'CROSSED_BY_THE_OUTLINE_OF', 'grid_lyr')
        arcpy.CopyFeatures_management('nhd_lyr', 'intersect.shp')

        #copy the model NHD streams to a temp shapefile, find the ends that match
        #the streams that were clipped by the boundary and update the lengthKm,
        #minsmoothelev, maxsmoothelev in the MODLNHD
        if arcpy.Exists('temp.shp'):
            print 'first removing old version of temp.shp'
            arcpy.Delete_management('temp.shp')
        arcpy.CopyFeatures_management(self.SFRdata.Flowlines, 'temp.shp')

        #add fields for start, end, and length to the temp and intersect
        #shapefiles  (use LENKM as field name because temp already has LENGHTKM)
        print "adding fields for start, end and length..."
        shapelist = ('temp.shp', 'intersect.shp')
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
        intersects = arcpy.SearchCursor('intersect.shp')
        manual_intervention = 0
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
                        ofp = open('boundary_manual_fix_issues.txt', 'w')
                        ofp.write('The following COMIDs identify streams that need manual attention.\n')
                        ofp.write('Fix in the files {0:s} and river_explode.shp. Then rerun intersect.py\n'.format(self.SFRdata.Flowlines))
                        ofp.write('#' * 16 + '\n')
                    print 'both ends are cut off for comid {0:d}\n'.format(comid)
                    ofp.write('both ends are cut off for comid {0:d}\n'.format(comid))
                    print 'need to fix it manually'
        del intersects, stream, clip, clippedstream
        if manual_intervention:
            ofp.write('#' * 16 + '\n')
            ofp.close()

        #create a new working NHD shapefile incorporating new values just found

        print "Creating new shapefile {0:s}".format(self.SFRdata.NHD)
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
        print "Saving routing information to {0:s}".format('boundaryClipsRouting.txt')
        ofp = open('boundaryClipsRouting.txt', 'w')
        ofp.write("FROMCOMID,TOCOMID\n")
        for currcomid in self.newCOMIDdata:
            if self.newCOMIDdata[currcomid].inout == 'IN':
                ofp.write("99999,{0:d}\n".format(self.newCOMIDdata[currcomid].comid))
            else:
                ofp.write("{0:d},99999\n".format(self.newCOMIDdata[currcomid].comid))
        ofp.close()
        if manual_intervention:
            print 'Some manual intervention required:\n' \
                  'See boundary_manual_fix_issues.txt for details'

    def make_rivers_table(self, FragIDdata):
        """
        from assign_and_route -->
        """
        print "Setting up the elevation table --> {0:s}".format(self.SFRdata.rivers_table)
        if arcpy.Exists(self.SFRdata.rivers_table):
            arcpy.Delete_management(self.SFRdata.rivers_table)
        arcpy.CreateTable_management(os.getcwd(), self.SFRdata.rivers_table)
        arcpy.AddField_management(self.SFRdata.rivers_table, "OLDFragID", "LONG")
        arcpy.AddField_management(self.SFRdata.rivers_table, "CELLNUM", "LONG")
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
                    outfile = open("fix_com_IDs.txt", 'w')
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
                orderedFragID.append(end_has_start[orderedFragID[i-1]])
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
                FragIDdata.allFragIDs[orderedFragID[i]].segelevinfo = avecellrivelev
                row = rows.newRow()
                row.OLDFragID = orderedFragID[i]
                row.CELLNUM = FragIDdata.allFragIDs[orderedFragID[i]].cellnum
                row.ELEVMAX = maxcellrivelev
                row.ELEVAVE = avecellrivelev
                row.ELEVMIN = mincellrivelev
                rows.insertRow(row)
        # write out the summary of comids to fix, then close the outfile
        if len(fix_comids_summary) > 0:
            print 'Some cells have multiple COMIDs entering and/or leaving.\n See file "fix_comids.txt"'
            outfile.write('#' * 30 + '\nSummary of COMIDS to fix:\n' +
                          'Delete these COMIDs from river_explode.shp, \n' +
                          'then run CleanupRiverCells.py and rerun Assign_and_Route.py\n')
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

def widthcorrelation(arbolate):
    #estimate widths, equation from Feinstein and others (Lake
    #Michigan Basin model) width=0.1193*(x^0.5032)
    # x=arbolate sum of stream upstream of the COMID in meters
    #NHDPlus has arbolate sum in kilometers.
    #print a table with reachcode, order, estimated width, Fcode
    estwidth = 0.1193*math.pow(1000*arbolate,0.5032)


"""
###################
ERROR EXCEPTION CLASSES
###################
"""
class InputFileMissing(Exception):
    def __init__(self, infile):
        self.infile = infile
    def __str__(self):
        return('\n\nCould not open or parse input file {0}\n'.format(self.infile))