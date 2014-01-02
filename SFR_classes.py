__author__ = 'Fienen, Reeves, Leaf - USGS'

import xml.etree.ElementTree as ET
import arcpy
import os


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
        self.reach_cutoff = float(inpars.findall('.//reach_cutoff')[0].text)
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
        self.NHD = inpars.findall('.//NHD')[0].text
        self.OUT = inpars.findall('.//OUT')[0].text
        self.MAT1 = inpars.findall('.//MAT1')[0].text
        self.MAT2 = inpars.findall('.//MAT2')[0].text
        self.WIDTH = inpars.findall('.//WIDTH')[0].text
        self.MULT = inpars.findall('.//MULT')[0].text
        try:
            self.eps = float(inpars.findall('.//eps')[0].text)
        except:
            self.eps = 1.0000001e-02 # default value used if not in the input file

    def tf2flag(self, intxt):
        # converts text written in XML file to True or False flag
        if intxt.lower() == 'true':
            return True
        else:
            return False


class COMIDProps:
    """
    Properties for each COMID
    """
    def __init__(self, comid, startx, starty, endx, endy, FID,
                 maxsmoothelev, minsmoothelev, lengthft, cellnum):
        self.comid = comid
        self.startx = startx
        self.starty = starty
        self.endx = endx
        self.endy = endy
        self.FID = FID
        self.maxsmoothelev = maxsmoothelev
        self.minsmoothelev = minsmoothelev
        self.lengthft = lengthft
        self.cellnum = cellnum
        self.elev = -9898989898
        self.from_comid = []
        self.to_comid = []


class COMIDPropsAll:
    def __init__(self):
        self.allcomids = dict()

    def populate(self, SFRdata):
        """
        read in the main COMID-related properties from the intersect file
        """
        segments = arcpy.SearchCursor(SFRdata.intersect)

        for seg in segments:
            comid = int(seg.COMID)
            self.allcomids[comid] = COMIDProps(
                int(comid),
                float(seg.X_start),
                float(seg.Y_start),
                float(seg.X_end),
                float(seg.Y_end),
                int(seg.FID),
                float(seg.MAXELEVSMO)*3.2808,
                float(seg.MINELEVSMO)*3.2808,
                float(seg.LengthFt),
                seg.node,
                )

    def populate_elevations(self, SFRdata):
        """
        Read elevation information, per COMID, from the SFRdata.ELEV file
        """
        arcpy.RefreshCatalog(os.getcwd())
        arcpy.env.qualifiedFieldNames = False
        arcpy.env.workspace = os.getcwd()
        if arcpy.Exists(SFRdata.ELEV):
            arcpy.Delete_management(SFRdata.ELEV)
        if arcpy.Exists("grid_temp"):
            arcpy.Delete_management("grid_temp")
        print "\nJoining {0:s} to {1:s}; \nsaving as {2:s} (might take a minute)".format(
              SFRdata.rivers_table, SFRdata.intersect, SFRdata.ELEV)
        arcpy.MakeFeatureLayer_management(SFRdata.intersect, "grid_temp")
        arcpy.AddJoin_management("grid_temp", "FID", SFRdata.rivers_table, "OLDFID")
        currelev = SFRdata.ELEV

        arcpy.CopyFeatures_management("grid_temp", currelev)

        with arcpy.da.SearchCursor(currelev, ("COMID", "ELEVAVE")) as cursor:
            for crow in cursor:
                self.allcomids[int(crow[0])].elev = float(crow[1])

    def populate_routing(self, SFRdata):
        """
        Read the COMID routing information from the SFRdata.FLOW file
        """
        allCOMIDs = self.allcomids.keys()
        print ('Reading in routing information from {0:s}'.format(SFRdata.FLOW))
        # open the SFRdata.FLOW file as read-only (using SearchCursor)
        with arcpy.da.SearchCursor(SFRdata.FLOW, ("FROMCOMID", "TOCOMID")) as cursor:
            for crow in cursor:
                if int(crow[0]) in allCOMIDs:
                    self.allcomids[int(crow[0])].to_comid.append(int(crow[1]))
                if int(crow[1]) in allCOMIDs:
                    self.allcomids[int(crow[1])].from_comid.append(int(crow[0]))

class SFRReachProps:
    """
    class containing just the data for each reach
    """
    def __init__(self):
        self.poo = poo


class SFRReachesAll:
    """
    class that makes up a list of SFRReachProps objects
    and also contains methods to access them
    """
    def __init__(self):
        j = 1


class InputFileMissing(Exception):
    def __init__(self, infile):
        self.infile = infile
    def __str__(self):
        return('\n\nCould not open or parse input file {0}\n'.format(self.infile))