__author__ = 'Fienen, Reeves, Leaf - USGS'

import xml.etree.ElementTree as ET
import arcpy
import os
import numpy as np
import SFR_arcpy
import time
import datetime
import discomb_utilities as disutil #  utility to read in a dis file


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

        #cutoff to check stream length in cell against fraction of cell dimension
        #if the stream length is less than cutoff*side length, the piece of stream is dropped
        try:
            self.cutoff = float(inpars.finall('.//cutoff')[0].text)
        except:
            self.cutoff = 0.0

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


class FIDprops(object):
    """
    Properties for each COMID
    """
    __slots__ = ['comid', 'startx', 'starty', 'endx', 'endy', 'FID',
                 'maxsmoothelev', 'minsmoothelev', 'lengthft',
                 'cellnum', 'elev', 'sidelength'
                 'segelevinfo', 'start_has_end', 'end_has_start','elev_distance']
    #  using __slots__ makes it required to declare properties of the object here in place
    #  and saves significant memory
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
        self.elev = None
        self.segelevinfo = None
        self.elev_distance = None


class COMIDprops(object):
    """
    routing information by COMIDs
    """
    __slots__ = ['from_comid', 'to_comid', 'hydrosequence', 'levelpathID']

    def __init__(self):
        self.from_comid = list()
        self.to_comid = list()
        self.hydrosequence = None
        self.levelpathID = None

class LevelPathIDprops(object):
    """
    routing of LevelPathIDs
    """
    __slots__ = ["down_levelpathID"]
    def __init__(self):
        self.down_levelpathID = None

class LevelPathIDpropsAll:
    def __init__(self):
        self.allids = dict()
        self.level_ordered = list()
        self.levelpath_fid = dict()

    def return_cutoffs(self, FIDdata, SFRdata):
        self.rmlist = []
        self.indata = SFRdata
        for lpID in self.level_ordered:
            for fid in self.levelpath_fid[lpID]:
                reachlength=FIDdata.allfids[fid].lengthft
                if reachlength < FIDdata.allfids[fid].sidelength*self.indata.cutoff:
                    print '{0}, {1}\n'.format(fid, self.indata.cutoff)
                    self.rmlist.append(fid)



class COMIDPropsAll:
    def __init__(self):
        self.allcomids = dict()
        self.hydrosequence_sorted = list()
        self.hydrosequence_comids = dict()

    def populate_routing(self, SFRdata, FIDdata, LevelPathdata):
        """
        Read the COMID routing information from the SFRdata.FLOW file
        """
        print ('Reading in routing information from {0:s}'.format(SFRdata.FLOW))
        # open the SFRdata.FLOW file as read-only (using SearchCursor)

        CLIP = np.loadtxt('boundaryClipsRouting.txt', skiprows=1, delimiter=',', dtype=int)

        for ccomid in FIDdata.allcomids:
            self.allcomids[ccomid] = COMIDprops()

        with arcpy.da.SearchCursor(SFRdata.FLOW, ("FROMCOMID", "TOCOMID")) as cursor:
            for crow in cursor:
                if int(crow[0]) in FIDdata.allcomids:
                    if (crow[0]) in CLIP[:, 0]:
                        self.allcomids[crow[0]].to_comid.append(999999)
                    else:
                        self.allcomids[crow[0]].to_comid.append(int(crow[1]))
                if int(crow[1]) in FIDdata.allcomids:
                    if crow[1] in CLIP[:, 1]:
                        self.allcomids[crow[1]].from_comid.append(999999)
                    else:
                        self.allcomids[crow[1]].from_comid.append(int(crow[0]))
        del crow, cursor
        comidseen = list()
        with arcpy.da.SearchCursor(SFRdata.PlusflowVAA,
                                   ("ComID", "Hydroseq", "LevelPathI", "UpLevelPat", "DnLevelPat")) as cursor:
            for crow in cursor:
                comid = int(crow[0])
                hydrosequence = int(crow[1])
                levelpathid = int(crow[2])
                if int(comid) in FIDdata.allcomids:
                    self.allcomids[crow[0]].hydrosequence = hydrosequence
                    self.allcomids[crow[0]].levelpathID = levelpathid
                    LevelPathdata.level_ordered.append(levelpathid)
                    comidseen.append(comid)

            # find unique levelpathIDs
            LevelPathdata.level_ordered = sorted(list(set(LevelPathdata.level_ordered)), reverse=True)

            for clevelpathid in LevelPathdata.level_ordered:
                LevelPathdata.allids[clevelpathid] = LevelPathIDprops()
                LevelPathdata.levelpath_fid[clevelpathid]=[]


            # assign levelpathID routing
        del crow, cursor

        with arcpy.da.SearchCursor(SFRdata.PlusflowVAA,
                                   ("ComID", "Hydroseq", "LevelPathI", "UpLevelPat", "DnLevelPat")) as cursor:
            for crow in cursor:
                comid = int(crow[0])
                levelpathid = int(crow[2])
                uplevelpathid = int(crow[3])
                downlevelpathid = int(crow[4])
                if levelpathid in LevelPathdata.level_ordered:
                    if downlevelpathid != levelpathid:
                        LevelPathdata.allids[levelpathid].down_levelpathID = downlevelpathid
                    if comid in FIDdata.comid_fid:
                        LevelPathdata.levelpath_fid[levelpathid].extend(FIDdata.comid_fid[comid])


        comid_missing = list(set(FIDdata.allcomids).difference(comidseen))
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


class FIDPropsAll:
    def __init__(self):
        self.allfids = dict()
        self.allcomids = list()  # comprehensive list of all comids
        self.unique_cells = list()  # list of unique cellnum values in the grid/streams intersection
        self.maxelev = None
        self.minelev = None
        self.noelev = dict()
        self.COMID_orderedFID = dict()
        self.comid_fid = None

    def return_fid_comid_list(self):
        """
        Return a dict of lists of comids linked with each FID
        """
        self.comid_fid = {ccomid: [] for ccomid in self.allcomids}
        allfids = self.allfids.keys()
        for cfid in allfids:
            self.comid_fid[self.allfids[cfid].comid].append(cfid)

    def return_smoothelev_comid(self, comid):
        self.maxelev = -np.inf
        self.minelev = np.inf
        fids = self.comid_fid[comid]
        for cid in fids:
            if self.allfids[cid].maxsmoothelev > self.maxelev:
                self.maxelev = self.allfids[cid].maxsmoothelev
            if self.allfids[cid].minsmoothelev < self.minelev:
                self.minelev = self.allfids[cid].minsmoothelev

    def return_unique_cells(self):
        for cid in self.allfids.keys():
            self.unique_cells.append(self.allfids[cid].cellnum)
        self.unique_cells = set(self.unique_cells)

    def populate(self, SFRdata):
        """
        read in the main COMID-related properties from the intersect file
        """
        segments = arcpy.SearchCursor(SFRdata.intersect)

        for seg in segments:
            fid = int(seg.FID)
            self.allfids[fid] = FIDprops(
                int(seg.COMID),
                float(seg.X_start),
                float(seg.Y_start),
                float(seg.X_end),
                float(seg.Y_end),
                int(seg.FID),
                float(seg.MAXELEVSMO)*3.2808,  # UNIT CONVERSION
                float(seg.MINELEVSMO)*3.2808,  # UNIT CONVERSION
                float(seg.LengthFt),
                seg.node)

            self.allcomids.append(int(seg.COMID))
        self.allcomids = list(set(self.allcomids))

    def populate_elevations(self, SFRdata):
        """
        Read elevation information, per COMID, from the SFRdata.ELEV file
        """


        arcpy.MakeFeatureLayer_management(SFRdata.intersect, "grid_temp")
        SFR_arcpy.general_join(SFRdata.ELEV, "grid_temp", "FID", SFRdata.rivers_table, "OLDFID", keep_common=False)

        with arcpy.da.SearchCursor(SFRdata.ELEV, ("FID", "ELEVAVE")) as cursor:
            for crow in cursor:
                self.allfids[int(crow[0])].elev = float(crow[1])


class SFRReachProps(object):
    """
    class containing just the data for each reach
    """
    __slots__ = ['junk']
    def __init__(self):
        self.junk = junk


class SFRReachesAll:
    """
    class that makes up a list of SFRReachProps objects
    and also contains methods to access them
    """
    def __init__(self):
        j = 1


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
        print 'verifying that there is a "node" field in {0:s}'.format(indat.MFgrid)
        hasnode = False
        nodeunique = 0
        MFgridflds = arcpy.ListFields(indat.MFgrid)
        for cfield in MFgridflds:
            if cfield.name.lower() == 'node':
                hasnode = True
        if hasnode:
            # now check to see that there are unique values in node
            cursor = arcpy.SearchCursor(indat.MFgrid)
            nodevals = []
            crows = 0
            for row in cursor:
                if crows > 100:
                    break
                else:
                    crows += 1
                    nodevals.append(row.getValue('node'))
            nodeunique = len(set(nodevals))
            del nodevals
            del row
            del cursor
        else:
            arcpy.AddField_management(indat.MFgrid, 'node', 'LONG')

        if nodeunique > 1:
            print '"node" field in place with unique values in {0:s}'.format(indat.MFgrid)
        else:
            arcpy.DeleteField_management(indat.MFgrid, 'node')
            arcpy.AddField_management(indat.MFgrid, 'node', 'LONG')
            print 'Updating "node" field in {0:s}'.format(indat.MFgrid)
            # now loop through and set "node" which is equivalent to "cellnum"
            cursor = arcpy.UpdateCursor(indat.MFgrid)
            for row in cursor:
                crow = row.getValue('row')
                ccol = row.getValue('column')
                row.setValue('node', ((int(crow)-1) * NCOL) + int(ccol))
                cursor.updateRow(row)
            del row
            del cursor





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
        arcpy.AddField_management(indat.CELLS, "CELLNUM", "LONG")
        arcpy.CalculateField_management(indat.CELLS, "CELLNUM", "!node!", "PYTHON")

        print "Dissolving river cells on cell number to isolate unique cells...\n"
        self.getfield(indat.CELLS, "node", "node")
        arcpy.Dissolve_management(indat.CELLS, indat.CELLS_DISS, self.joinnames['node'])

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
        self.getfield(indat.intersect, "node", "node")
        self.getfield(indat.intersect, "lengthft", "Length")
        table = arcpy.UpdateCursor(indat.intersect)
        count = 0
        for reaches in table:
            if reaches.getValue(self.joinnames["Length"]) <= indat.reach_cutoff:
                print "segment: %d, cell: %s, length: %s" % (reaches.getValue(self.joinnames["comid"]),
                                                            reaches.getValue(self.joinnames["node"]),
                                                            reaches.getValue(self.joinnames["Length"]))
                self.ofp.write("segment: %d, cell: %s, length: %s\n"
                          % (reaches.getValue(self.joinnames["comid"]),
                            reaches.getValue(self.joinnames["node"]),
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
                                 "node",
                                 "river_explode",
                                 "node",
                                 "KEEP_ALL")  # this might not work as-in in stand-alone mode
        self.getfield("river_cells_dissolve", "node", "node")
        self.getfield("river_cells_dissolve", "maxelevsmo", "maxelevsmo")
        table = arcpy.SearchCursor("river_cells_dissolve")
        nodes2delete = []
        for row in table:
            if row.isNull(self.joinnames["maxelevsmo"]):
                nodes2delete.append(row.getValue(self.joinnames["node"]))
        arcpy.RemoveJoin_management("river_cells_dissolve", "river_explode")

        # remove nodes with no elevation information from river_explode
        self.ofp.write('\n' + 25*'#' + '\nRemoving nodes with no elevation information from river_explode\n')
        print 'Removing nodes with no elevation information from river_explode'
        self.getfield(indat.CELLS_DISS, "node", "node")
        table = arcpy.UpdateCursor(indat.CELLS_DISS)
        count = 0
        for cells in table:
            if cells.getValue(self.joinnames["node"]) in nodes2delete:
                print "%d" % (cells.getValue(self.joinnames["node"]))
                self.ofp.write('%d\n' % (cells.getValue(self.joinnames["node"])))
                table.deleteRow(cells)
                count += 1
        print "removed %s cells\n" % (count)
        self.ofp.write("removed %s cells" % (count))

        print "removing any remaining disconnected reaches..."
        self.ofp.write('\n' + 25*'#' + '\nremoving any remaining disconnected reaches...\n')
        self.getfield(indat.intersect, "node", "node")
        table = arcpy.UpdateCursor(indat.intersect)
        count = 0
        for reaches in table:
            if reaches.getValue(self.joinnames["node"]) in nodes2delete:
                print "%d" % (reaches.getValue(self.joinnames["node"]))
                self.ofp.write('%d\n' % (reaches.getValue(self.joinnames["node"])))
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

    def make_rivers_table(self, FIDdata):
        """
        from assign_and_route -->
        """
        print "Setting up the elevation table --> {0:s}".format(self.SFRdata.rivers_table)
        if arcpy.Exists(self.SFRdata.rivers_table):
            arcpy.Delete_management(self.SFRdata.rivers_table)
        arcpy.CreateTable_management(os.getcwd(), self.SFRdata.rivers_table)
        arcpy.AddField_management(self.SFRdata.rivers_table, "OLDFID", "LONG")
        arcpy.AddField_management(self.SFRdata.rivers_table, "CELLNUM", "LONG")
        arcpy.AddField_management(self.SFRdata.rivers_table, "ELEVMAX", "DOUBLE")
        arcpy.AddField_management(self.SFRdata.rivers_table, "ELEVAVE", "DOUBLE")
        arcpy.AddField_management(self.SFRdata.rivers_table, "ELEVMIN", "DOUBLE")
        segcounter = 0
        rows = arcpy.InsertCursor(self.SFRdata.rivers_table)
        fix_comids_summary = []
        fixcomids_flag = False
        for comid in FIDdata.allcomids:
            segcounter += 1
            fidlist = FIDdata.comid_fid[comid]
            start_has_end = dict()
            end_has_start = dict()

            for i in range(0, len(fidlist)):
                haveend = False
                havestart = False
                for j in range(0, len(fidlist)):
                    if j == i:
                        continue
                    diffstartx = FIDdata.allfids[i].startx - FIDdata.allfids[j].endx
                    diffstarty = FIDdata.allfids[i].starty - FIDdata.allfids[j].endy
                    diffendx = FIDdata.allfids[i].endx - FIDdata.allfids[j].startx
                    diffendy = FIDdata.allfids[i].endy - FIDdata.allfids[j].starty
                    if np.fabs(diffstartx) < self.SFRdata.rfact and np.fabs(diffstarty) < self.SFRdata.rfact:
                        start_has_end[fidlist[i]] = fidlist[j]
                        haveend = True
                    if np.fabs(diffendx) < self.SFRdata.rfact and np.fabs(diffendy) < self.SFRdata.rfact:
                        end_has_start[fidlist[i]] = fidlist[j]
                        havestart = True
                    if haveend and havestart:
                        break

            #find key in start_has_end that didn't match an end and
            #key in end_has_start that didn't match a start
            numstart = 0
            numend = 0
            startingFID = []
            endingFID = []
            for test in fidlist:
                if test not in start_has_end:
                    startingFID.append(test)
                    numstart += 1
                if test not in end_has_start:
                    endingFID.append(test)
                    numend += 1
            if numstart != 1 or numend != 1:
                if not fixcomids_flag:
                    outfile = open("fix_com_IDs.txt", 'w')
                    fixcomids_flag = True
                outfile.write("numstart =" + str(numstart) + " \n")
                outfile.write("numend = " + str(numend) + " \n")
                outfile.write("starting FIDs: " + ",".join(map(str, startingFID)) + "\n")
                outfile.write("ending FIDs: " + ",".join(map(str, endingFID)) + "\n")
                outfile.write("manually fix COMID = %d\n" %comid)
                fix_comids_summary.append('%d\n' %comid)
                FIDdata.noelev[comid] = 1  #set flag
                continue

            orderedFID = []
            orderedFID.append(startingFID[0])
            for i in range(1, len(end_has_start)):
                orderedFID.append(end_has_start[orderedFID[i-1]])
            if orderedFID[-1] != endingFID[0]:       #  don't repeat the last entry FID...
                orderedFID.append(endingFID[0])
            #total length read through lengthkm didn't always match up
            #to the sum of the lengths of the segments (exactly), sumup total length
            totallength = 0
            for i in range(0, len(orderedFID)):
                totallength += FIDdata.allfids[orderedFID[i]].lengthft
            if totallength == 0:
                exit('check length ft for FIDs in COMID= %d' % comid)
            FIDdata.return_smoothelev_comid(comid)
            slope = (FIDdata.maxelev - FIDdata.minelev)/totallength
            distance = 0.
            FIDdata.COMID_orderedFID[comid] = orderedFID
            for i in range(0, len(orderedFID)):
                maxcellrivelev = FIDdata.maxelev - slope*distance
                distance += FIDdata.allfids[orderedFID[i]].lengthft
                mincellrivelev = FIDdata.maxelev - slope*distance
                avecellrivelev = 0.5*(maxcellrivelev + mincellrivelev)
                FIDdata.allfids[orderedFID[i]].segelevinfo = avecellrivelev
                row = rows.newRow()
                row.OLDFID = orderedFID[i]
                row.CELLNUM = FIDdata.allfids[orderedFID[i]].cellnum
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

    def reach_ordering(self):
        """
        crawl through the hydrosequence values and
        set a preliminary reach ordering file
        """
        indat = self.SFRdata
        ofp = open(indat.RCH, 'w')
        ofp.write('CELLNUM,COMID, hydroseq, uphydroseq, dnhydroseq, '
                  'levelpathID, uplevelpath, dnlevelpath, SFRseq, localseq\n')

class Elevs_from_contours:
    def __init__(self,SFRdata,FIDdata):
        self.contours_intersect_distances = SFRdata.Contours_intersect_distances
        self.FIDdata = FIDdata

    def assign_elevations_to_FID(self,COMIDdata):

        # loop through rows in contours intersect table
        distances = arcpy.SearchCursor(self.contours_intersect_distances)
        elevs_edited = []
        for row in distances:
            comid = int(row.COMID)
            cellnum = int(row.node)

            # loop through FIDdata: if comid and cellnum match, update elevation from contours intersect table
            for fid in self.FIDdata.allfids.keys():
                if self.FIDdata.allfids[fid].cellnum == cellnum and self.FIDdata.allfids[fid].comid == comid:
                    self.FIDdata.allfids[fid].elev == float(row.ContourEle)
                    self.FIDdata.allfids[fid].elev_distance = float(row.fmp)
                    elevs_edited.append(fid)

        # reset all elevations not updated from the contours to 0
        for fid in self.FIDdata.allfids.keys():
            if fid not in elevs_edited:
                self.FIDdata.allfids[fid].elev = 0

        # update elevations by COMID
        for comid in self.FIDdata.COMID_orderedFID.keys():
            from_comid = COMIDdata.allcomids[comid].from_comid
            to_comid = COMIDdata.allcomids[comid].to_comid

            # within each COMID, iterate going upstream (reversed)
            start_elev,end_elev = None,None
            dist = 0
            interp = False

            def get_dist_slope(FIDdata,comid,start_elev,end_elev,dist,interp):

                for fid in self.FIDdata.COMID_orderedFID[comid].reverse():

                    # if the FID intersects a contour, record elevation and position
                    if not interp and fid in elevs_edited:
                        end_elev = self.FIDdata.allfids[fid].elev
                        dist +=self.FIDdata.allfids[fid].elev_distance
                        interp = True
                        continue

                    # between contours, add distance
                    if interp and fid not in elevs_edited:
                        dist += self.FIDdata.allfids[fid].lengthft
                        continue

                    # at upstream contour, record elevation and calculate slope
                    if interp and fid in elevs_edited:
                        start_elev = self.FIDdata.allfids[fid].elev
                        dist += self.FIDdata.allfids[fid].lengthft - self.FIDdata.allfids[fid].elev_distance
                        slope = (start_elev - end_elev)/dist
                        interp = False
                        return dist,slope,interp

                    else: #upstream contour not in current comid
                        slope = None
                        return dist,slope,interp


            dist,slope,end_elev,interp = get_dist_slope(self.FIDdata,comid,start_elev,end_elev,dist,interp)

            # if upstream contour not in current comid, go to from_comid
            if interp and slope == None:
                while interp:
                    upcomid = COMIDdata.allcomids[comid].from_comid
                    dist,slope,end_elev,interp = get_dist_slope(self.FIDdata,upcomid,start_elev,end_elev,dist,interp)











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