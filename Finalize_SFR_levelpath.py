# FinalizeRiverInput_SFR2.py
# Description: Takes river_w_elevation shapefile after RouteStreamNetwork.py
# and assigns final cell elevations.  Also uses routed_cells.txt and reach_ordering.txt
# that are output from Assign_and_Route.py
#
# New version using levelpathID information, now provided by Assign_and_Route,
# to build SFR segments and reaches
#
# Output file for the SFR2 package
#
# Requirements: os, re, arcpy, defaultdict, itemgetter, math, numpy, pdb
#
# Author: Howard Reeves; USGS Michigan Water Science Center
#         Mike Fienen and Andy Leaf, USGS Wisconsin Water Science Center
# Date: 1/2/2014
#
import os
import re
import arcpy
import numpy as np
from collections import defaultdict
from operator import itemgetter
import math
import pdb

# Global Input file for SFR utilities
infile="SFR_setup.in"

# Get input files locations
infile=open(infile,'r').readlines()
inputs=defaultdict()
for line in infile:
    if "#" not in line.split()[0]:
        varname=line.split("=")[0]
        var=line.split("=")[1].split()[0]
        inputs[varname]=var.replace("'","") # strip extra quotes

ELEV=inputs["ELEV"]
ROUTE="routed_cells.txt"
ORDER="reach_ordering.txt"
VAA=inputs["PlusflowVAA"]
FTab=inputs["FTab"]
MFgrid=inputs["MFgrid"]
#nrows=inputs["ROWS"]
#numcol=inputs["COLUMNS"]

# Output file
OUT=inputs["OUT"]
GISSHP=inputs["GISSHP"]
MAT1=inputs["MAT1"]
MAT2=inputs["MAT2"]
WIDTH=inputs["WIDTH"]
MULT=inputs["MULT"]
CELLS=inputs["CELLS"]      #used as template for GISSHP

# Step 14 in Howard's SFR notes
#arcpy.JoinField_management(ELEV,"node","river_cells.shp","node")

print "getting grid information.."
numrow=arcpy.SearchCursor(MFgrid,"","","","row D").next().getValue("row")
numcol=arcpy.SearchCursor(MFgrid,"","","","column D").next().getValue("column")

print '%s rows, %s columns' %(numrow,numcol)
'''
def getgridinfo(MFgrid):
    print "getting grid information.."
    
    table=arcpy.SearchCursor(MFgrid)
    numrow=[]
    numcols=[]
    delx=[]
    dely=[]
    for row in table:
        numrow.append(row.getValue("row"))
        numcols.append(row.getValue("column"))
        delx.append(row.getValue("delx"))
        dely.append(row.getValue("dely"))

    nrows=np.max(numrow)
    numcol=np.max(numcols)

    print "\t%s rows x %s columns" %(nrows,numcol)
    if len(np.unique(delx))==1:
        spacing="uniform"
        print "\tuniform spacing of %s " %(delx[1])
    else:
        spacing="non-uniform"
    return(nrows,numcol)'''

   

#set cut off for total stream length to include cell (summed over
#all the segments (fraction of side edge as a decimal value)
cutoff=0.0

#choose how to assign river elevation, max, ave or min

elev_type='min'

#set a minimum slope for a cell if the computed slope is very small or zero
#such as for cells in a meander

minslope=0.00001

#set elevation difference to assign same value to all cells
#so that check does not get messed up by rounding floats
eps=1.0e-02

#set arcpy path and environment
path=os.getcwd()
arcpy.env.workspace=path

#delete any working layers if found
if arcpy.Exists("temp_lyr"):
    arcpy.Delete_management("temp_lyr")

#go through river cells.  Use cell_num index as dictionary key and
#use defaultdict for river information in case there are multiple rivers
#in the cell.

reachlength=defaultdict(list)
riverelev=defaultdict(list)
old_fid=defaultdict(list)
old_fid_comid=dict()
row=dict()
column=dict()
sidelength=dict()
comid=defaultdict(list)
cellslope=defaultdict(list)
comidseen=dict()
newreachlength=dict()     #try setting up dicts keyed by cellnum+comid (as strings)
newcellslope=dict()
newcellelev=dict()

cellrows=arcpy.SearchCursor(ELEV)
for cell in cellrows:
    cellnum=int(cell.CELLNUM)
    comidin=int(cell.comid)
    comid[cellnum].append(comidin)
    comidseen[comidin]=1
    row[cellnum]=cell.row
    column[cellnum]=cell.column
    sidelength[cellnum]=float(cell.delx)
    if float(cell.dely) < sidelength[cellnum]:
        sidelength[cellnum]=float(cell.dely)
    reachlength[cellnum].append(cell.LengthFT)
    if re.match("min",elev_type,re.I):
        elevchoice=cell.ELEVMIN
    elif re.match("max",elev_type,re.I):
        elevchoice=cell.ELEVMAX
    elif re.match("ave",elev_type,re.I):
        elevchoice=cell.ELEVAVE
    else:
        print "check elevchoice\n"
        exit()
    if cell.LengthFT>0.000001:
        localslope=(cell.ELEVMAX - cell.ELEVMIN)/cell.LengthFT
    else:
        localslope=0.
    if localslope < 0.:
        print 'check elevations, min>max, cellnum= %d' % cellnum
        #exit('check elevations, min > max')
    cellslope[cellnum].append(localslope)
    riverelev[cellnum].append(elevchoice)
    old_fid[cellnum].append(cell.ORIG_FID)
    old_fid_comid[cell.ORIG_FID]=comidin
    #dictionaries storing cellnum/comid pair information using
    #cellnum+comid (as strings) for the key; makes assembly easier
    #than going through defaultdict lists that are keyed by cellnum
    #and have comid and other information that has to be retrieved by
    #order in the list
    newkey=str(int(cell.CELLNUM))+str(int(cell.comid))
    newreachlength[newkey]=float(cell.LengthFT)
    newcellslope[newkey]=localslope
    newcellelev[newkey]=elevchoice
    
del cellrows

#read in NHDPlus tables for arbolate sum and stream order
#and assign to cells by COMID, then use to estimate
#stream width for each COMID, stored by CELLNUM for printing

#make a dictionary of arbolate sum, stream order, Fcode, and reach code by COMID
#also get hydroseq, uphydroseq, and dnhydroseq
reachcode=dict()
streamorder=dict()
arbolate=dict()
Fcode=dict()
Fcodeseen=dict()

segments=arcpy.SearchCursor(VAA)
for segment in segments:
    comidseg=int(segment.comid)
    if comidseg in comidseen:
        reachcode[comidseg]=long(segment.ReachCode)
        streamorder[comidseg]=int(segment.StreamOrde)
        arbolate[comidseg]=float(segment.ArbolateSu)
        Fcode[comidseg]=int(segment.Fcode)
        Fcodeseen[Fcode[comidseg]]=1
del segments

#get text that corresponds to the Fcode
Fstring=dict()
descrips=arcpy.SearchCursor(FTab)
for description in descrips:
    Fcodevalue=int(description.FCode)
    if Fcodevalue in Fcodeseen:
        Fstring[Fcodevalue]=description.Descriptio
del descrips

#estimate widths, equation from Feinstein and others (Lake
#Michigan Basin model) width=0.1193*(x^0.5032)
# x=arbolate sum of stream upstream of the COMID in meters
#NHDPlus has arbolate sum in kilometers.
#print a table with reachcode, order, estimated width, Fcode

widthout=open(WIDTH,'w')
widthout.write('cellnum,row,column,comid,stream_order,arbolate_sum,est_width,reachlength,REACHCODE,FCode,Description\n')

estwidth=defaultdict(list)
newestwidth=dict()
for cellnum in row:
    for i in range(0,len(comid[cellnum])):
        # added 1000 to convert from km to m for arbolate sum (widths were too small otherwise)
        widthcorrelation=0.1193*math.pow(1000*arbolate[comid[cellnum][i]],0.5032)
        estwidth[cellnum].append(widthcorrelation)
        #dictionary keyed by cellnum+comid as strings...
        newkey=str(cellnum)+str(comid[cellnum][i])
        newestwidth[newkey]=widthcorrelation
    for i in range(0,len(comid[cellnum])):
        comidcell=comid[cellnum][i]
        printstring=(cellnum,
                     row[cellnum],
                     column[cellnum],
                     comidcell,
                     streamorder[comidcell],
                     arbolate[comidcell],
                     estwidth[cellnum][i],
                     reachlength[cellnum][i],
                     reachcode[comidcell],
                     Fcode[comidcell],
                     Fstring[Fcode[comidcell]])
        widthout.write(",".join(map(str,printstring))+'\n')
widthout.close()

#get the cell-to-cell routing that was generated by RouteRiverCells.py

routetable=open(ROUTE,'r')
routetable.readline()   # get the header
vals=[]
fromCell=defaultdict(list)
for line in routetable:
    vals=re.split(",",line,2)
    fromCell[int(vals[0])].append(int(vals[1]))

            
#use levelpathID and hydrosequence numbering to generate lists of segments and
#reaches in hydrosequence order; get the information for SFR2 that was
#generated by Assign_and_Route.py, make dictionaries linking
#cellnumbers, hydrosequence numbers, and levelpathIDs; after unique
#levelpathIDs are identified, assign each levelpathID to a SFR segment
#and build the reaches.  The levelpathID is the downstream hydrosequence
#number for a stretch of the stream, so working in decending order will
#build the SFR segments in downstream order.

ordering=open(ORDER,'r')
ordering.readline()  # skip the header
hydrocomid=dict()        #keyed by comid
inv_hydrocomid=dict()    #keyed by hydrosequence
orderedcell=defaultdict(list)  #list of cells by comid already in downstream order
dnhydroseq=dict()          #keyed by comid -> gives downstream hydrosequence, used in confluence determination
uplevelpath=dict()
dnlevelpath=dict()
levelpathpair=defaultdict(list)  #levelpathpair will be a dictionary of lists, each entry will be [comid,cellnum]
levelpathlist=[]
for line in ordering:
    vals=re.split(",",line)
    cellnum=int(vals[0])
    localcomid=int(vals[1])
    localhydroseq=int(vals[2])
    dnhydroseq[localcomid]=int(vals[4])
    hydrocomid[localcomid]=localhydroseq
    inv_hydrocomid[localhydroseq]=localcomid
    orderedcell[localcomid].append(cellnum)
    lpID=int(vals[5])
    levelpathlist.append(lpID)
    uplevelpathin=int(vals[6])   #individual entries point to itself, only one points to the next up or down levelpathID
    if uplevelpathin!=lpID:
        uplevelpath[lpID]=uplevelpathin
    dnlevelpathin=int(vals[7])
    if dnlevelpathin!=lpID:
        dnlevelpath[lpID]=dnlevelpathin
    levelpathpair[lpID].append([cellnum,localcomid])
ordering.close()

#sort the levelpathIDs in downstream order and get the unique list of levelpathIDs
level_ordered=sorted(set(list(levelpathlist)), reverse=True)

#loop over levelpathIDs and check reachlengths in cells; if any reachlength
#is less than the input cutoff, delete that piece from the list.
rmpair=defaultdict(list)
for lpID in level_ordered:
    for entry in range(0,len(levelpathpair[lpID])):
        lpcell=levelpathpair[lpID][entry][0]
        lpcomid=levelpathpair[lpID][entry][1]
        for i in range(0,len(comid[lpcell])):
            locallength=reachlength[lpcell][i]
            localcomid=comid[lpcell][i]
            if localcomid==lpcomid:
                if (locallength < sidelength[lpcell]*cutoff):
                    #length of comid,cellnum in the levelpath is short, remove that
                    #comid and cellnum from the levelpath lists
                    rmpair[lpID].append([lpcell,lpcomid])

for lpID in rmpair:
    #find indices of matches
    indices=set()
    for i in range(0,len(rmpair[lpID])):
        rmcell=rmpair[lpID][i][0]
        rmcomid=rmpair[lpID][i][1]
        for j in range(0,len(levelpathpair[lpID])):
            lpcell=levelpathpair[lpID][j][0]
            lpcomid=levelpathpair[lpID][j][1]
            if rmcell==lpcell and rmcomid==lpcomid:
                indices.add(j)
    #now use enumerate and list comprehension to build a new list without
    #the pairs that matched
    newlist=[pair for i, pair in enumerate(levelpathpair[lpID]) if i not in indices]
    levelpathpair[lpID]=newlist

#If any levelpath ends up empty, remove it from the level_ordered list
indices=set()
for i in range(0,len(level_ordered)):
    if len(levelpathpair[level_ordered[i]])==0:
        indices.add(i)
newlist=[lp for i, lp in enumerate(level_ordered) if i not in indices]
level_ordered=newlist

#loop over level_ordered and assign each as a SFR segment,
#build the reaches from the hydrosequences and comids within each
#levelpathID
SFRprovreachlist=dict()                # dictionary keyed by levelpathID
                                        # each entry is list of cellnumbers
SFRprovseg=dict()                      # dictionary keyed by levelpathID, value is provisional segment number
SFRprovcomid=dict()                    # final comid in an SFR segment, keyed by levelpathID
totalreach=0
for i in range(0,len(level_ordered)):
    iseg=i+1   # segment ordering starts at 1
    irch=0
    lpID=level_ordered[i]
    SFRprovseg[lpID]=iseg
    #make a dictionary of cell/hydrosequences in the levelpath by getting the comids from levelpair
    hydrodict=dict()
    for j in range(0,len(levelpathpair[lpID])):
        lpcell=levelpathpair[lpID][j][0]
        lpcomid=levelpathpair[lpID][j][1]
        hydrodict[lpcomid]=hydrocomid[lpcomid]   #key=comid, value=hydrosequence
        totalreach+=1
    #now sort hydrodict(which is just for the current levelpath) in downstream order
    hydro_ordered=hydrodict.items()
    hydro_ordered.sort(key=itemgetter(1), reverse=True)
    hydro_orderedcomids=map(itemgetter(0), hydro_ordered)
    levelpathcells=[]
    for localcomid in hydro_orderedcomids:
        levelpathcells.extend(orderedcell[localcomid])
    #get the unique cells in order, smooths away the meandering.
    uniq=[]
    for cell in levelpathcells:
        if cell not in  uniq:
            uniq.append(cell)
    SFRprovreachlist[lpID]=uniq
    SFRprovcomid[lpID]=hydro_orderedcomids[-1]

#go through level_ordered, for each levelpathID, there may be more than one
#comid in a cell.  Sum the reachlengths and get a weighted width and slope, elevation
#put into defaultdict(dict) keyed by lpID and cellnumber
totlength=defaultdict(dict)
weightwidth=defaultdict(dict)
elevcell=defaultdict(dict)
weightedslope=defaultdict(dict)
for lpID in level_ordered:
    segment=SFRprovseg[lpID]
    for localcell in SFRprovreachlist[lpID]:
        tt=0.
        ww=0.
        ws=0.
        el=0.
        for entry in range(0,len(levelpathpair[lpID])):
            lpcell=levelpathpair[lpID][entry][0]
            lpcomid=levelpathpair[lpID][entry][1]
            knt=0
            if lpcell==localcell:
                knt=knt+1
                newkey=str(lpcell)+str(lpcomid)
                tt=tt+newreachlength[newkey]
                ww=ww+newestwidth[newkey]*newreachlength[newkey]
                el=el+newcellelev[newkey]
                ws=ws+newcellslope[newkey]*newreachlength[newkey]
        if knt==0:
            totlength[lpID][localcell]=reachlength[localcell][0]
            elevcell[lpID][localcell]=riverelev[localcell][0]
            weightedslope[lpID][localcell]=cellslope[localcell][0]
            weightwidth[lpID][localcell]=estwidth[localcell][0]
        else:
            totlength[lpID][localcell]=tt
            elevcell[lpID][localcell]=el/float(knt)
            if tt>0:
                weightwidth[lpID][localcell]=ww/tt
                weightedslope[lpID][localcell]=ws/tt
            else:
                weightwidth[lpID][localcell]=99999.
                weightedslope[lpID][localcell]=99999.

nss=iseg  #store the number of segments
print "number of segments = %d and number of reaches = %d" % (nss,totalreach)
print 'now going to generate SFR files'

#go through levelpaths, find the cells where there are confluences and
#make a list of cells of confluences for each segment (levelpathID) that will
#be used to subdivide and renumber to get the final SFR segments

CHK2=open('segment_confluences.txt','w')
CHK2.write('segment,levelpathID,num_confluences,cells\n')

confluences=defaultdict(list)
for i in range(0,nss):
    iseg=i+1
    lpID=level_ordered[i]
    if not lpID in dnlevelpath:
        outseg=int(0)
    else:
        nextlevelpath=dnlevelpath[lpID]
        if nextlevelpath in SFRprovseg:
            outseg=SFRprovseg[nextlevelpath]
        else:
            outseg=int(0)
    isegend=SFRprovreachlist[lpID][-1]
    if outseg>0:
        outID=level_ordered[outseg-1]
        outsegbeg=SFRprovreachlist[outID][0]
        confl=-1                                        #flag in case an end is not found
        for j in range(0,len(SFRprovreachlist[outID])):
            if SFRprovreachlist[outID][j]==isegend:
                confl=SFRprovreachlist[outID][j]
    else:
        #no downstream levelpathID
        outsegbeg=0
        outID=0
        confl=0
    #no end found
    if confl==-1:
        lastcomid=SFRprovcomid[lpID]
        #check to see if dnhydroseq =0, if so the last cell of iseg is a downstream end within the model
        if lastcomid in dnhydroseq:
            if dnhydroseq[lastcomid]==0:
                confl=-2
    #put confluence cells into lists, iseg end is automatically a confluence for iseg
    confluences[lpID].append(isegend)
    #now check the downstream (receiving) segment, only put one on a list if confl is
    #greater than zero; otherwise confl is a downsteam outlet (or an error indicated by confl=-1)
    if confl>0:
        confluences[outID].append(confl)
    if confl==-1:
        print 'check downstream connection for levelpathID %d' % lpID
    #finalreachlist is in hydrosequence order, put confluences in hydrosequence order with no repeats from confluences
    confluences[lpID]=[cell for cell in SFRprovreachlist[lpID] if cell in set(confluences[lpID])]

subseg=dict()
subconfl=defaultdict(list)
subcell=dict()
subreaches=dict()
conf_count=1

for i in range(0,nss):
    iseg=i+1
    lpID=level_ordered[i]
    numconfls=len(confluences[lpID])
    #break up segments with multiple confluences
    strt=0
    plist1=(iseg,SFRprovreachlist[lpID])
    #CHK2.write(",".join(map(str,plist1))+'\n')
    #CHK2.write("**"+","+str(len(confluences[lpID]))+","+",".join(map(str,confluences[lpID]))+"\n")
    for confl in range(0,numconfls):
        sublabel=str(iseg)+'-'+str(confl)
        subseg[sublabel]=conf_count
        subconfl[iseg].append(sublabel)
        subcell[sublabel]=confluences[lpID][confl]
        conf_count+=1
        #build reach lists for each subsection of a provisional segment
        #find the index that matches the cell given by confluences[lpID][confl]
        endindx=SFRprovreachlist[lpID].index(confluences[lpID][confl])
        subreaches[sublabel]=SFRprovreachlist[lpID][strt:endindx+1]
        strt=endindx+1
        plist2=(sublabel,subseg[sublabel],subreaches[sublabel])
        #CHK2.write(",".join(map(str,plist2))+'\n')

# use confluence information and build final SFR segments and reaches

SFRfinalreachlist=dict()          #list of cells (reaches) keyed by final segment number
SFRfinaloutseg=dict()             #dictionary of outsegment keyed by final segment number
SFRfinalseg=dict()                #dictionary of segment labels from confluence step keyed by final segment number
SFRfinalID=dict()                 #dictionary of levelpathID keyed by final segment number
for i in range(0,nss):
    iseg=i+1
    lpID=level_ordered[i]
    for upstlabl in subconfl[iseg]:
        SFRfinalreachlist[subseg[upstlabl]]=subreaches[upstlabl]
        SFRfinalseg[subseg[upstlabl]]=upstlabl
        SFRfinalID[subseg[upstlabl]]=lpID
        lastcell=subreaches[upstlabl][-1]
        #find the provisional segment that iseg is connected
        if not lpID in dnlevelpath:
            outseg=int(0)
        else:
            nextlevelpath=dnlevelpath[lpID]
            if nextlevelpath in SFRprovseg:
                outseg=SFRprovseg[nextlevelpath]
            else:
                outseg=int(0)
        if outseg>0:
            outID=level_ordered[outseg-1]
            if outseg in subconfl:
                #find confluence
                printed=False
                for dnlabl in subconfl[outseg]:
                    #see if the beginning of an outsegment matches the end of the current segment
                    #or if the cell connected to the end of the current segment is the beginning
                    #of a segment (no overlap, end of current segment is connected to beginning of
                    #next segment in an adjacent cell)
                    overlapcell=subreaches[dnlabl][0]
                    if overlapcell==lastcell:
                        SFRfinaloutseg[subseg[upstlabl]]=subseg[dnlabl]
                        plist=('overlap',iseg, upstlabl, subseg[upstlabl], lastcell, outseg, dnlabl, subseg[dnlabl], subreaches[dnlabl])
                        CHK2.write(",".join(map(str,plist))+'\n')
                        printed=True
                    if not printed:
                        for nxtdwnstream in set(fromCell[lastcell]):
                            if nxtdwnstream==overlapcell:
                                SFRfinaloutseg[subseg[upstlabl]]=subseg[dnlabl]
                                plist=('offset',iseg, upstlabl, subseg[upstlabl], lastcell, outseg, dnlabl, subseg[dnlabl], subreaches[dnlabl])
                                CHK2.write(",".join(map(str,plist))+'\n')
                                printed=True

                if not printed:
                    SFRfinaloutseg[subseg[upstlabl]]=int(99999)
                    plist=('no connection',iseg,upstlabl,subseg[upstlabl],lastcell,99999,99999)
                    CHK2.write(",".join(map(str,plist))+'\n')
            else:
                SFRfinaloutseg[subseg[upstlabl]]=int(0)
                plist=('downstream end',iseg,upstlabl,subseg[upstlabl],lastcell,0,0)
                CHK2.write(",".join(map(str,plist))+'\n')
        else:
            SFRfinaloutseg[subseg[upstlabl]]=int(0)
            plist=('downstream end',iseg,upstlabl,subseg[upstlabl],lastcell,0,0)
            CHK2.write(",".join(map(str,plist))+'\n')


CHK2.write('segment, confluence-based label, outsegment, cells(reaches) in segment\n')
for seg in SFRfinalseg.iterkeys():
    CHK2.write("%d, %s, %d, " % (seg, SFRfinalseg[seg], SFRfinaloutseg[seg]))
    CHK2.write("; ".join(map(str,SFRfinalreachlist[seg]))+"\n")

CHK2.close()

# have it all now by cellnumber... loop over lists by hydrosequence and generate SFR2 package
try:
    outfile=open(OUT,'w')
except:
    print "could not open output file"
    exit()
outfile.write(r'# SFR2 Package generated by python scrips from grid shapefile and NHDPlus data'+'\n')

#make a shapefile of river cells and save result--
#more than one entry will be assigned to each cell (potentially)
#so a simple join cannot be used.
if arcpy.Exists(GISSHP):
    arcpy.Delete_management(GISSHP)
#CELLS shapefile used to create GISSHP then add necessary fields, 
print 'making shapefile %s' % GISSHP
arcpy.CreateFeatureclass_management(path,GISSHP,"POLYGON",CELLS)
arcpy.AddField_management(GISSHP,"CELLNUM","LONG")
arcpy.AddField_management(GISSHP,"ROW","LONG")
arcpy.AddField_management(GISSHP,"COLUMN","LONG")
arcpy.AddField_management(GISSHP,"LAYER","LONG")
arcpy.AddField_management(GISSHP,"SEGMENT","LONG")
arcpy.AddField_management(GISSHP,"REACH","LONG")
arcpy.DeleteField_management(GISSHP,"node")
newrows=arcpy.InsertCursor(GISSHP)
shapeName=arcpy.Describe(ELEV).shapeFieldName

#open GWV files
mat1out=open(MAT1,'w')
mat1out.write('row,column,layer,stage,top_streambed,reach,segment,width_in_cell,length_in_cell,')
mat1out.write('bed_K,bed_thickness,bed_slope,bed_roughness\n')

mat2out=open(MAT2,'w')
mat2out.write('segment,icalc,outseg,iupseg,iprior,nstrpts,flow,runoff,etsw,pptsw,')
mat2out.write('roughch,roughbk,cdepth,fdepth,awdth,bwdth\n')

#SFR2 values - defaults from Daniel Feinstein
nsfrpar=0
nparseg=0
const=1.486*86400.
dleak=0.0001
nstrail=10
isuzn=1
nsfrsets=30
istcb1=50
istcb2=66
isfropt=1
bedK=5.0
bedKmin=1.0e-08
bedthick=1.0
icalc=1
nstrpts=0
iprior=0
flow=0.
runoff=0.
etsw=0.
pptsw=0.
roughch=0.037
roughbk=0.
cdepth=0.
fdepth=0.
awdth=0.
bwdth=0.
thickm1=1.0
thickm2=1.0
Hc1fact=1.0
Hc2fact=1.0

#reset nss to final segment number
nss=len(SFRfinalseg)
printstring=(-totalreach,nss,nsfrpar,nparseg,const,dleak,istcb1,istcb2,isfropt,nstrail,isuzn,nsfrsets)
outfile.write(' '.join(map(str,printstring))+'\n')

for iseg in SFRfinalseg.iterkeys():
    levelpathcells=SFRfinalreachlist[iseg]
    lpID=SFRfinalID[iseg]
    progprint=True
    for j in range(0,len(levelpathcells)):
        irch=j+1
        localcell=levelpathcells[j]
        printstring=(int(1),  
                        row[localcell],
                        column[localcell],
                        iseg,
                        irch)
        outfile.write(' '.join(map(str,printstring)))
        #assume top of streambed is 1 ft below elevation from NHDPlus
        floatlist=[totlength[lpID][localcell],
                        elevcell[lpID][localcell]-1.0,
                        weightedslope[lpID][localcell],
                        bedthick,
                        bedK]
        mixedlist=(elevcell[lpID][localcell],
               elevcell[lpID][localcell]-1.0,
               irch,
               iseg,
               weightwidth[lpID][localcell],
               totlength[lpID][localcell],
               bedK,
               bedthick,
               weightedslope[lpID][localcell],
               roughch)
        printstring=' {0:.2f} {1:.2f} {2:.3e} {3:.2f} {4:.2e}'.format(*floatlist)
        outfile.write(printstring + '\n')

        #GWV matrix 1 file
        printstring=printstring=(  
                        row[localcell],
                        column[localcell],
                        int(1),)
        mat1out.write(','.join(map(str,printstring)))
        
        printstring=',{0:.2f},{1:.2f},{2:d},{3:d},{4:.2f},{5:.2f},{6:.2e},{7:.2f},{8:.3e},{9:.4f}'.format(*mixedlist)
        mat1out.write(printstring+'\n')
        
        #GISSHP information
        #get polygon from CELLS
        query="CELLNUM=%d"%localcell
        if (iseg % int(nss/5)==0) and progprint:
            print "%d percent done with shapefile %s"%(int(math.ceil(float(i)/nss*100.)),GISSHP)
            progprint=False
        poly=arcpy.SearchCursor(CELLS,query)
        newvals=newrows.newRow()
        for entry in poly:
            feature=entry.getValue(shapeName)
            newvals.Shape=feature
            newvals.CELLNUM=localcell
            newvals.ROW=row[localcell]
            newvals.COLUMN=column[localcell]
            newvals.LAYER=int(1)
            newvals.SEGMENT=iseg
            newvals.REACH=irch
            newrows.insertRow(newvals)
        
del newvals, newrows
arcpy.RefreshCatalog(path)

#continue with output files
outfile.write('1 0 0 0\n')   # item 5...

maxcells=numrow*numcol

for iseg in SFRfinalseg.iterkeys():
    levelpathcells=SFRfinalreachlist[iseg]
    lpID=SFRfinalID[iseg]
    outseg=SFRfinaloutseg[iseg]
    iupseg=0
    printstring=[iseg,icalc,outseg,iupseg]       #iupseg set to zero right now, no diversions, could be added....
##    if iupseg > 0:
##        printstring.append(iprior)
    if icalc==4:
        printstring.append(nstrpts)
    printstring.extend([runoff,etsw,pptsw])
    if icalc==1 or icalc==2:
        printstring.append(roughch)
    outfile.write(' '.join(map(str,printstring))+'\n')
    #for SFR2, K and thick are in previous line (reach-by-reach) only width is needed here...
    startcell=SFRfinalreachlist[iseg][0]
    endcell=SFRfinalreachlist[iseg][-1]
    printstring='{0:.1f}'.format(weightwidth[lpID][startcell])  #width1...
    outfile.write(printstring+'\n')
    printstring='{0:.1f}'.format(weightwidth[lpID][endcell])  #width2...
    outfile.write(printstring +'\n')
    #write ouput to GWV matrix 2 file
    mlist1=(iseg,icalc,outseg,iupseg,iprior,nstrpts)
    mlist2=(flow,runoff,etsw,pptsw,roughch,roughbk)
    mlist3=(cdepth,fdepth,awdth,bwdth)
    printstring='{0:d},{1:d},{2:d},{3:d},{4:d},{5:d}'.format(*mlist1)
    printstring=printstring+',{0:.2f},{1:.2f},{2:.2f},{3:.2f},{4:.4f},{5:.2f}'.format(*mlist2)
    printstring=printstring+',{0:.2f},{1:.2f},{2:.2f},{3:.2f}'.format(*mlist3)
    mat2out.write(printstring+'\n')

outfile.close()
mat1out.close()
mat2out.close()