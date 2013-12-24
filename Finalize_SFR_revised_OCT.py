# FinalizeRiverInput_SFR2.py
# Description: Takes river_w_elevation shapefile after checking RouteStreamNetwork.py
# and assigns final cell elevations.  Also uses routed_cells.txt and reach_ordering.txt
# that are output from
# RouteRiverCells.py
#
# Revised to work with the output from 'Assign and Route' script
#
# Output file for the SFR2 package
#
# Requirements: os, sys, re, arcpy, defaultdict, itemgetter, math
#
# Author: H.W. Reeves; USGS Michigan Water Science Center
# Date 10/28/13
#
import os
import sys
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
cutoff=0.1

#choose how to assign river elevation, max, ave or min

elev_type='min'

#set a minimum slope for a cell if the computed slope is very small or zero
#such as for cells in a meander

minslope=0.00001

#set elevation difference to assign same value to all cells
#so that check does not get messed up by rounding floats
eps=1.0e-02

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
#if a cell has more than one comid, assign dominantcomid as the
#one with the largest width*reachlength
widthout=open(WIDTH,'w')
widthout.write('cellnum,row,column,comid,stream_order,arbolate_sum,est_width,REACHCODE,FCode,Description\n')

dominantcomid=dict()
estwidth=defaultdict(list)
for cellnum in row:
    for i in range(0,len(comid[cellnum])):
         estwidth[cellnum].append(0.1193*math.pow(1000*arbolate[comid[cellnum][i]],0.5032)) # added 1000 to convert from km to m (widths were too small otherwise)

    biggest=0.
    for i in range(0,len(comid[cellnum])):
        comidcell=comid[cellnum][i]
        if (estwidth[cellnum][i]*reachlength[cellnum][i])>biggest:
            biggest=(estwidth[cellnum][i]*reachlength[cellnum][i])
            dominantcomid[cellnum]=comidcell
        printstring=(cellnum,
                     row[cellnum],
                     column[cellnum],
                     dominantcomid[cellnum],
                     streamorder[comidcell],
                     arbolate[comidcell],
                     estwidth[cellnum][i],
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

###make a defaultdict list of meander cells.
###the file ORDER has provisional SFR segment and reaches
###can use this file to make lists of celles within a meander.
##    
##ordering=open(ORDER,'r')
##ordering.readline()  # skip the header
##meandercell=defaultdict(list)
##segmentcells=defaultdict(list)
##cell_segment=defaultdict(list)
##seg_comid=dict()
##
##for line in ordering:
##    vals=re.split(',',line.strip())
##    cellnum=int(vals[0])
##    provsegment=int(vals[5])
##    provreach=int(vals[6])
##    segmentcells[provsegment].append(cellnum)
##    cell_segment[cellnum].append(provsegment)
##    seg_comid[provsegment]=int(vals[1])
##ordering.close()
##
##for cellnum in old_fid.iterkeys():
##    if cellnum in meandercell:
##        # check if already got it with another FID
##        continue
##    #see if there are repeated FIDS in the same cell using the set() data type
##    reptFID=set([x for x in old_fid[cellnum] if old_fid[cellnum].count(x)>1])
##    if len(reptFID)>0:
##        #assemble a list cells between the two entries
##        #that have the same FID in the cell
##        tempsegment=cell_segment[cellnum]
##        #see if any cell segments in tempsegment are repeated 
##        rept=set([x for x in tempsegment if tempsegment.count(x)>1])
##        #find the indices of the segmentcells list that match the cellnumber with a meander
##        for seg in list(rept):
##            j=[indx for indx in range(len(segmentcells[seg])) if segmentcells[seg][indx]==cellnum]
##            #get a slice of the vector from first to last index of the identified list
##            tmpslice=segmentcells[seg][j[0]:j[-1]]
##            #append the comid of the segment to the end of the list of cells
##            tmpslice.append(seg_comid[seg])
##            if len(tmpslice)>2:   #sometimes divergences had one-cell long meanders in the test, skip 
##                meandercell[cellnum].append(tmpslice)
##            
##del tempsegment,seg_comid
##   
### go through meander cells and make sure there are no
### repeats within the same meandercell list: set() takes the unique values from the original list
### and then list() turns the unique set back into a list
##for cellnum in meandercell:
##    for i in range(0,len(meandercell[cellnum])):
##        meanderlist=meandercell[cellnum][i][0:-1]
##        mlistcomid=meandercell[cellnum][i][-1]
##        uniq=set(meanderlist)
##        meandercell[cellnum][i]=sorted(list(uniq))
##        meandercell[cellnum][i].append(mlistcomid)
##
##print 'done with meander step'
### go through cells and assign elevation and conductance, also
### use information from shapefiles of cells tops and bottoms
### to determine layer of the SFR2 cell
##cellseen=defaultdict(dict)
##condfactor=defaultdict(dict)
##totlength=defaultdict(dict)
##weightwidth=defaultdict(dict)
##elevcell=defaultdict(dict)
##weightedslope=defaultdict(dict)
##meandernum=0
##for cellnum in meandercell.iterkeys():
##    meandernum=meandernum+1
##    for mlist in range(0,len(meandercell[cellnum])):
##        #comid has been appended to the list of cells in each meandercell entry, pull it out
##        meanderlist=meandercell[cellnum][mlist][0:-1]
##        mlistcomid=meandercell[cellnum][mlist][-1]
##        minelev=9999999
##        for cell in meanderlist:
##            if not cell in riverelev:
##                continue
##            elif min(riverelev[cell]) < minelev:
##                minelev=min(riverelev[cell])
##        ttlen=0.
##        sumcondfact=0.
##        ww=0.
##        i=0
##        for i in range(0,len(reachlength[cell])):
##            sumcondfact=sumcondfact+reachlength[cell][i]*estwidth[cell][i]
##            ttlen=ttlen+reachlength[cell][i]
##            ww=ww+estwidth[cell][i]*reachlength[cell][i]
##        reachout=0
##        layer=1
##        for cell in meanderlist:
##            if ttlen > (cutoff*sidelength[cell]):
##                totlength[cell][mlistcomid]=ttlen
##                condfactor[cell][mlistcomid]=sumcondfact
##                elevcell[cell][mlistcomid]=minelev
##                weightedslope[cell][mlistcomid]=minslope   # put a minimum slope in for the meanders
##                ww=ww/ttlen
##                weightwidth[cell][mlistcomid]=ww
##                cellseen[cell][mlistcomid]=1
##            #done with meanders
##
##multiple_segs=dict()
##cols=[]
##for cellnum in row.iterkeys():
##    multsegs=0
##    for i in range(0,len(riverelev[cellnum])):
##        if reachlength[cellnum][i]>(cutoff*sidelength[cellnum]):
##            comidcell=comid[cellnum][i]
##            if cellnum in cellseen:
##                if comidcell in cellseen[cellnum]: 
##                    continue  #skip if already set in the meander step....
##            multsegs+=1
##            condfactor[cellnum][comidcell]=reachlength[cellnum][i]*estwidth[cellnum][i]
##            elevcell[cellnum][comidcell]=riverelev[cellnum][i]
##            totlength[cellnum][comidcell]=reachlength[cellnum][i]
##            weightwidth[cellnum][comidcell]=estwidth[cellnum][i]
##            weightedslope[cellnum][comidcell]=cellslope[cellnum][i]
##            if weightedslope[cellnum][comidcell] < minslope:
##                weightedslope[cellnum][comidcell]=minslope
##    if multsegs>=2:
##        multiple_segs[cellnum]=multsegs
###go through the cells with multiple entries and identify the one with the largest
###conductance - set it as the dominantcomid, used in the printing step to choose bigK or bigKmin
##dominantcomid=dict()
##for cellnum in multiple_segs:
##    maxcond=0.
##    for comidcell in condfactor[cellnum].iterkeys():
##        if condfactor[cellnum][comidcell]>maxcond:
##            maxcond=condfactor[cellnum][comidcell]
##            dominantcomid[cellnum]=comidcell
            
#use hydrosequence numbering to generate lists of segments and reaches in hydrosequence order
#get the hydrosequence and local sequence numbering for SFR2 that was
#generated by RouteRiverCells.py, make dictionaries linking
#cellnumber to both SFR2 sequence number and SFR2 reach number; and also
#make list of cellnumbers with the same SFR2 sequence number
#hydro_orderedcells=[(cells in seq.1),(cells in seq.2)....]
#put them into a default dict and then sort by keys, the sort returns
#a list of key-value tuples -> grab the values which will be the individual
#lists of cellnumbers

ordering=open(ORDER,'r')
ordering.readline()  # skip the header
hydrowork=defaultdict(list)
hydroseq=dict()
uphydroseq=dict()
dnhydroseq=dict()
inv_hydroseq=dict()
hydrocomid=dict()
for line in ordering:
    vals=re.split(",",line)
    cellnum=int(vals[0])
    localcomid=int(vals[1])
    sfrsegment=int(vals[5])
    hydrocomid[sfrsegment]=localcomid
    hydroseq[sfrsegment]=int(vals[2])
    hydrowork[hydroseq[sfrsegment]].append(cellnum)
    inv_hydroseq[int(vals[2])]=sfrsegment
    uphydroseq[sfrsegment]=int(vals[3])
    dnhydroseq[sfrsegment]=int(vals[4])

ordering.close()               

#now go through hydrowork and remove any keys with empty lists of cells
for rawsfrsegment in hydrowork:
    if not hydrowork[rawsfrsegment]:
        print 'raw sfr segment dropped', rawsfrsegment
        del hydrowork[rawsfrsegment]

#now sort the remaining entries of the dictionary in descending (downstream) order and
#then use the map command to convert to a list of list of cells
hydro_ordered=hydrowork.items()
hydro_ordered.sort(key=itemgetter(0),reverse=True)
hydro_orderedcells=map(itemgetter(1),hydro_ordered)  #list of the lists of cells
hydro_orderedindex=map(itemgetter(0),hydro_ordered)  #list of the indexes [hydrosequence] 

print 'done with hydro_ordered step'

#now the raw sfr segments from RouteRiverCells are ordered by hydrosequence
#number and there are dictionaries of hydrosequences (up and downstream) and
#hydrosequence-sfrsegment; sfrsegment-comid;
#need to tie these to the dictionaries of cell-by-cell properties created earlier
#that included the filter for short segments, this could remove cells and
#maybe entire sfr segments.  Allow for renumbering of segments and reaches
#and march through to make final SFR output.

#the shapefile cell_inout has a list of cells and whether each cell
#has one comid (same); a single comid in and a single comid out but
#two or more comids in the cell (onein/oneout); or multiple comids
#either entering or exiting the cell (multiple).  The last case
#could be a convergence or non-intersecting streams.

#loop over hydro_orderedindex and renumber SFR segment and reaches
#if one if missing.  Also tag cells that have been output already and
#do not repeat writing that cell.

#includes a dictionary of the connections for each segment.  The dictionary
#key is the SFR final segment number.  The dictionary item is a list with
#two entries, the first entry is the upstream segment, zero for headwater,
#and the second entry is the downstream segment, zero for ends- outflows to boundary
#or other feature where the SFR segment chain ends.  The lambda in the
#variable definition sets [0,0] as the default for a segment that has not
#yet been assigned a value.

#read in cell_inout
celltype=dict()
cellcomidin=defaultdict(list)
cellcomidout=defaultdict(list)
cellnumpieces=dict()
cellinfo=arcpy.SearchCursor('cell_inout.dbf')
for cell in cellinfo:
    cellnum=int(cell.CELLNUM)
    celltype[cellnum]=cell.TYPE
    cellnumpieces[cellnum]=cell.COUNT
    cellcomidin[cellnum]=cell.COMIDIN
    cellcomidout[cellnum]=cell.COMIDOUT
del cellinfo

#for cells with more than one piece of a COMID, get total length
#and weighted width, comid is a dictionary keyed by cell number
#with all the comids in the cell.  cellcomidin and cellcomidout
#get converted to lists of the comids entering and leaving the cell

weightwidth=dict()
totlength=dict()
weightedslope=dict()
elevcell=dict()
for cell in celltype.iterkeys():
    if not cell in comid:
        print 'mismatch in cell numbering from river elevation and cell_inout'
        print 'check river_w_elevations, that the cellnum entry is not zero'
        print 'for the line with node = %d' % cell
    elif re.match('multiple',celltype[cell]) or re.match('onein/oneout',celltype[cell]):
        if re.search(',',cellcomidin[cell]):
            comidinlist=map(int,re.split(',',cellcomidin[cell]))
        else:
            comidinlist=[int(cellcomidin[cell])]
        if re.search(',',cellcomidout[cell]):
            comidoutlist=map(int,re.split(',',cellcomidout[cell]))
        else:
            comidoutlist=[int(cellcomidout[cell])]
        ttlen=0.
        ww=0.
        ws=0.
        el=0.
        for i in range(0,len(comid[cell])):
            ttlen=ttlen+reachlength[cell][i]
            ww=ww+estwidth[cell][i]*reachlength[cell][i]
            ws=ws+cellslope[cell][i]*reachlength[cell][i]
            el=el+riverelev[cell][i]
        elevcell[cell]=el/len(comid[cell])
        weightwidth[cell]=ww/ttlen
        totlength[cell]=ttlen
        weightedslope[cell]=ws/ttlen
    else:
        totlength[cell]=reachlength[cell][0]
        weightwidth[cell]=estwidth[cell][0]
        weightedslope[cell]=cellslope[cell][0]
        elevcell[cell]=riverelev[cell][0]

SFRfinalsegment=dict()
SFRfinalhydroseq=dict()
SFRfinalreachlist=defaultdict(list)
SFRfinalcomid=dict()
SFRconnect=defaultdict(lambda:[0,0])
ordereduphydro=dict()
ordereddnhydro=dict()
iseg=1
totalreach=0

cellseen=dict()   # keep track of which cells have been encoutered
                  # key is the cell number, value is the segment number
for i in range(0,len(hydro_orderedindex)):
    irch=0
    for j in range(0,len(hydro_orderedcells[i])):
        localcell=hydro_orderedcells[i][j]
        if localcell not in cellseen:
            if (totlength[localcell]<sidelength[localcell]*cutoff):    #skip cells with short lengths
                break
            cellseen[localcell]=iseg
            irch=irch+1
            totalreach=totalreach+1
            SFRfinalsegment[iseg]=hydro_orderedindex[i]
            SFRfinalhydroseq[hydro_orderedindex[i]]=iseg
            SFRfinalcomid[iseg]=hydrocomid[inv_hydroseq[hydro_orderedindex[i]]]
            SFRfinalreachlist[iseg].append(localcell)
            ordereduphydro[iseg]=uphydroseq[inv_hydroseq[hydro_orderedindex[i]]]
            ordereddnhydro[iseg]=dnhydroseq[inv_hydroseq[hydro_orderedindex[i]]]
        else:
            SFRconnect[iseg][1]=cellseen[localcell]
    if irch>0:
        iseg=iseg+1

nss=iseg-1
print "number of segments = %d and number of reaches = %d" % (nss,totalreach)
print 'now going to generate SFR files'

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


printstring=(-totalreach,nss,nsfrpar,nparseg,const,dleak,istcb1,istcb2,isfropt,nstrail,isuzn,nsfrsets)
outfile.write(' '.join(map(str,printstring))+'\n')

for i in range(0,nss):
    iseg=i+1
    localcomid=SFRfinalcomid[iseg]
    alreadyprinted=False
    for j in range(0,len(SFRfinalreachlist[iseg])):
        irch=j+1
        localcell=SFRfinalreachlist[iseg][j]
        printstring=(int(1),  
                        row[localcell],
                        column[localcell],
                        iseg,
                        irch)
        outfile.write(' '.join(map(str,printstring)))
        #assume top of streambed is 1 ft below elevation from NHDPlus
        if localcell in dominantcomid:
            if dominantcomid[localcell]==localcomid:
                floatlist=[totlength[localcell],
                            elevcell[localcell]-1.0,
                            weightedslope[localcell],
                            bedthick,
                            bedK]
                mixedlist=(elevcell[localcell],
                   elevcell[localcell]-1.0,
                   irch,
                   iseg,
                   weightwidth[localcell],
                   totlength[localcell],
                   bedK,
                   bedthick,
                   weightedslope[localcell],
                   roughch)
            else:
                 floatlist=[totlength[localcell],
                            elevcell[localcell]-1.0,
                            weightedslope[localcell],
                            bedthick,
                            bedKmin]
                 mixedlist=(elevcell[localcell],
                   elevcell[localcell]-1.0,
                   irch,
                   iseg,
                   weightwidth[localcell],
                   totlength[localcell],
                   bedKmin,
                   bedthick,
                   weightedslope[localcell],
                   roughch)
        else:
            floatlist=[totlength[localcell],
                            elevcell[localcell]-1.0,
                            weightedslope[localcell],
                            bedthick,
                            bedK]
            mixedlist=(elevcell[localcell],
                   elevcell[localcell]-1.0,
                   irch,
                   iseg,
                   weightwidth[localcell],
                   totlength[localcell],
                   bedK,
                   bedthick,
                   weightedslope[localcell],
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
        query="node=%d"%localcell
        if (i % int(nss/10)==0) and alreadyprinted==False:
            print "%d percent done with shapefile %s"%(int(float(i)/nss*100.),GISSHP)
            alreadyprinted=True
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
#build a dictionary of cells that have a beginning segment
begincell=defaultdict(list)
for i in range(0,nss):
    iseg=i+1
    localbegin=SFRfinalreachlist[iseg][0]
    begincell[localbegin].append(iseg)

maxcells=numrow*numcol

for i in range(0,nss):
    iseg=i+1
    localcomid=SFRfinalcomid[iseg]
    dnhydro=ordereddnhydro[iseg]
    if dnhydro in SFRfinalhydroseq:
        outseg=SFRfinalhydroseq[dnhydro]
    else:
        outseg=int(0)
    if outseg==0:
        #one last check that this is really a terminal segment,
        #check the eight surrounding cells and see if any are
        #SFR cells, then see if it is the start of another segment
        #If there is a match, use the segment for outseg.
        if iseg==9029:
            print 'in eightchecks'
        endcellnumber=SFRfinalreachlist[iseg][-1]
        #skip if row=1, row=numrow; col=1 or col=numcol
        if row[endcellnumber]>1 and row[endcellnumber]<numrow and column[endcellnumber]>1 and column[endcellnumber]<numcol:
            eightchecks=[endcellnumber-1,endcellnumber+1,endcellnumber+numcol,endcellnumber-numcol,
                         endcellnumber+numcol-1,endcellnumber+numcol+1,endcellnumber-numcol-1,endcellnumber-numcol+1]
            #remove cellnumbers that are not feasible
            tcheck=[y for y in eightchecks if y>0]
            tcheck2=[y for y in tcheck if y<maxcells]
            
            eightchecks=tcheck2
            
            for checkcell in eightchecks:
                if checkcell in begincell:
                    for j in range(0,len(begincell[checkcell])):
                        if begincell[checkcell][j]!=iseg:
                            outseg=begincell[checkcell][j]
                            if iseg==9029:
                                print iseg, checkcell, begincell[checkcell], outseg
    #check if a segment should be connected to the downstream end of the outseg -
    #see if the distance between the two beginning cells is > a couple of cells
    #then check if the next outseg (outseg of the outseg) is near the end of the
    #segment and reset the outseg if necessary.
    endcellnumber=SFRfinalreachlist[iseg][-1]
    if outseg !=0:
        nextcellnumber=SFRfinalreachlist[outseg][0]
        diffrow=row[endcellnumber]-row[nextcellnumber]
        diffcol=column[endcellnumber]-column[nextcellnumber]
        if math.fabs(diffrow)>2 or math.fabs(diffcol)>2:
            #see if the other end of the outseg will work
            provcellnumber=SFRfinalreachlist[outseg][-1]
            provrow=row[endcellnumber]-row[provcellnumber]
            provcol=column[endcellnumber]-column[provcellnumber]
            if math.fabs(provrow)<=2 or math.fabs(provcol)<=2:
                nextdwnhydro=ordereddnhydro[outseg]
                if iseg==9029:
                    print outseg, nextdwnhydro
                if nextdwnhydro in SFRfinalhydroseq:
                    outseg=SFRfinalhydroseq[nextdwnhydro]
                    if iseg==9029:
                        print iseg, endcellnumber, provcellnumber, outseg
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
    printstring='{0:.1f}'.format(weightwidth[startcell])  #width1...
    outfile.write(printstring+'\n')
    printstring='{0:.1f}'.format(weightwidth[endcell])  #width2...
    outfile.write(printstring +'\n')
    #write ouput to GWV matrix 2 file
    mlist1=(iseg,icalc,outseg,iupseg,iprior,nstrpts)
    mlist2=(flow,runoff,etsw,pptsw,roughch,roughbk)
    mlist3=(cdepth,fdepth,awdth,bwdth)
    printstring='{0:d},{1:d},{2:d},{3:d},{4:d},{5:d}'.format(*mlist1)
    printstring=printstring+',{0:.2f},{1:.2f},{2:.2f},{3:.2f},{4:.4f},{5:.2f}'.format(*mlist2)
    printstring=printstring+',{0:.2f},{1:.2f},{2:.2f},{3:.2f}'.format(*mlist3)
    mat2out.write(printstring+'\n')

#print a table with reachcode, order, estimated width, Fcode
widthout=open(WIDTH,'w')
widthout.write('cellnum,row,column,comid,stream_order,arbolate_sum,est_width,reach_length,REACHCODE,FCode,Description,Segment,Reach\n')
for i in range(0,nss):
    iseg=i+1
    localcomid=SFRfinalcomid[iseg]
    for j in range(0,len(SFRfinalreachlist[iseg])):
        irch=j+1
        localcell=SFRfinalreachlist[iseg][j]
        if Fcode[localcomid] in Fstring:
            printstring=(localcell,
                     row[localcell],
                     column[localcell],
                     localcomid,
                     streamorder[localcomid],
                     arbolate[localcomid],
                     weightwidth[localcell],
                     totlength[localcell],
                     reachcode[localcomid],
                     Fcode[localcomid],
                     Fstring[Fcode[localcomid]],
                     iseg,
                     irch,)
        else:
            printstring=(cellnum,
                     row[cellnum],
                     column[cellnum],
                     comidcell,
                     streamorder[localcomid],
                     arbolate[localcomid],
                     weightwidth[localcell],
                     totlength[localcell],
                     reachcode[localcomid],
                     Fcode[localcomid],
                     "Unknown",
                     iseg,
                     irch,)
        widthout.write(",".join(map(str,printstring))+'\n')

#close files
widthout.close()    
outfile.close()
mat1out.close()
mat2out.close()

