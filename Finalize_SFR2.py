# FinalizeRiverInput_SFR2.py
# Description: Takes river_w_elevation shapefile after checking RouteStreamNetwork.py
# and assigns final cell elevations.  Also uses routed_cells.txt that is output from
# RouteRiverCells.py
#
# Output file for the SFR2 package
#
# Requirements: os, sys, re, arcpy, defaultdict, itemgetter, math
#
# Author: H.W. Reeves; USGS Michigan Water Science Center
# Date 7/17/12
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
GIS=inputs["GIS"]
MAT1=inputs["MAT1"]
MAT2=inputs["MAT2"]
WIDTH=inputs["WIDTH"]

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
cutoff=0.25

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
    print "cellnum: %s, comid: %s" %(cellnum,comidin)
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
        exit('check elevations, min > max')
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

#make a defaultdict list of meander cells.
#the file ORDER has provisional SFR segment and reaches
#can use this file to make lists of celles within a meander.
    
ordering=open(ORDER,'r')
ordering.readline()  # skip the header
meandercell=defaultdict(list)
segmentcells=defaultdict(list)
cell_segment=defaultdict(list)

for line in ordering:
    vals=re.split(',',line.strip())
    cellnum=int(vals[0])
    provsegment=int(vals[5])
    provreach=int(vals[6])
    segmentcells[provsegment].append(cellnum)
    cell_segment[cellnum].append(provsegment)
ordering.close()

for cellnum in old_fid.iterkeys():
    #sort the list for the cellnum and then see if there are repeats.
    sortedfids=sorted(old_fid[cellnum])
    i=0
    end=len(sortedfids)-1
    for i in range(0,end):
        if (sortedfids[i+1] == sortedfids[i]):
            #check to see if its the dominant COMID for the cell
            tempcomid=old_fid_comid[sortedfids[i]]
            if dominantcomid[cellnum]==tempcomid:
                #assemble a list cells between the two entries
                #that have the same FID in the cell
                tempsegment=cell_segment[cellnum]
                seensegment=dict()
                for i in range(0,len(tempsegment)):
                    if not tempsegment[i] in seensegment:
                        seensegment[tempsegment[i]]=1
                        #find the indices of the segmentcells list that match the cellnumber with a meander
                        j=[indx for indx in range(len(segmentcells[tempsegment[i]])) if segmentcells[tempsegment[i]][indx]==cellnum]
                        #get a slice of the vector from first to last index of the identified list
                        meandercell[cellnum]=segmentcells[tempsegment[i]][j[0]:j[-1]]
            else:
                #meandercell has different dominant comid, skip as a meander
                #print '{1:d} not dominant comid in meander cell {0:d}'.format(cellnum, tempcomid)
                continue

del seensegment, tempsegment
# go through meander cells and identify if any cells within the lists
# are in more than one meander.  If so combine
seencells=defaultdict(list)
for keycellnum, cells in meandercell.iteritems():
    for i in range(len(cells)):
        seencells[cells[i]].append(keycellnum)

alreadymerged=dict()
for cell,meanderkey in seencells.iteritems():
    mergedcells=[]
    if len(meanderkey)>1:
        #cell is in more than one meander list; each one
        #is represented by meander key
        for i in range(0,len(meanderkey)):
            mergedcells.extend(meandercell[meanderkey[i]])
        meandercell[meanderkey[0]]=mergedcells
        

# go through cells and assign elevation and conductance, also
# use information from shapefiles of cells tops and bottoms
# to determine layer of the SFR2 cell
cellseen=dict()
condfactor=dict()
totlength=dict()
weightwidth=dict()
elevcell=dict()
weightedslope=dict()
meandernum=0
for cellnum in meandercell.iterkeys():
    meandernum=meandernum+1
    meanderlist=meandercell[cellnum]
    minelev=9999999
    for cell in meanderlist:
        if cell in cellseen:
            continue
        if min(riverelev[cell]) < minelev:
            minelev=min(riverelev[cell])
        totlength[cell]=0.
        weightwidth[cell]=0.
        sumcondfact=0.
        i=0
        for i in range(0,len(reachlength[cell])):
            sumcondfact=sumcondfact+reachlength[cell][i]*estwidth[cell][i]
            totlength[cell]=totlength[cell]+reachlength[cell][i]
            weightwidth[cell]=weightwidth[cell]+estwidth[cell][i]*reachlength[cell][i]
        reachout=0
        layer=1
        if totlength[cell] > (cutoff*sidelength[cell]):   
            condfactor[cell]=sumcondfact
            elevcell[cell]=minelev
            weightedslope[cell]=minslope   # put a minimum slope in for the meanders
            weightwidth[cell]=weightwidth[cell]/totlength[cell]
            cellseen[cell]=1
#done with meanders
multiple_segs=dict()
cols=[]
for cellnum in row.iterkeys():
    if cellnum in cellseen:
        continue
    elif len(riverelev[cellnum])>1:
        multiple_segs[cellnum]=len(riverelev[cellnum])
    elif reachlength[cellnum][0]>(cutoff*sidelength[cellnum]):
        condfactor[cellnum]=reachlength[cellnum][0]*estwidth[cellnum][0]
        elevcell[cellnum]=riverelev[cellnum][0]
        totlength[cellnum]=reachlength[cellnum][0]
        weightwidth[cellnum]=estwidth[cellnum][0]
        weightedslope[cellnum]=cellslope[cellnum][0]
        if weightedslope[cellnum] < minslope:
            weightedslope[cellnum]=minslope
        cellseen[cellnum]=1

for cellnum in multiple_segs.iterkeys():
    if cellnum in cellseen:
        continue
    sumcondfact=0.
    totlength[cellnum]=0.
    weightwidth[cellnum]=0.
    aveelev=9999999
    sumslope=0.
    for i in range(0,len(riverelev[cellnum])):
        sumcondfact=sumcondfact+reachlength[cellnum][i]*estwidth[cellnum][i]
        totlength[cellnum]=totlength[cellnum]+reachlength[cellnum][i]
        weightwidth[cellnum]=weightwidth[cellnum]+estwidth[cellnum][i]*reachlength[cellnum][i]
        sumslope=sumslope+cellslope[cellnum][i]*reachlength[cellnum][i]
        if riverelev[cellnum][i] < aveelev:
            aveelev=riverelev[cellnum][i]
        
    if totlength[cellnum] > (cutoff*sidelength[cellnum]):
        condfactor[cellnum]=sumcondfact
        elevcell[cellnum]=aveelev
        weightwidth[cellnum]=weightwidth[cellnum]/totlength[cellnum]
        weightedslope[cellnum]=sumslope/totlength[cellnum]
        if weightedslope[cellnum] < minslope:
            weightedslope[cellnum]=minslope
        cellseen[cellnum]=1

#use dominantcomid for each cell and hydrosequence numbering
#to generate lists of segments and reaches in hydrosequence order
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
cellhydros=defaultdict(list)
inv_hydroseq=dict()
for line in ordering:
    vals=re.split(",",line)
    cellnum=int(vals[0])
    localcomid=int(vals[1])
    #if the localcomid is the dominant comid for the cellnumber,
    #then continue, otherwise skip because another comid is used
    #for this cell; also check if the cell was seen above- some
    #cells are dropped because of the length criteria
    if dominantcomid[cellnum]==localcomid and cellnum in cellseen:
        sfrsegment=int(vals[5])
        hydrowork[sfrsegment].append(cellnum)
        hydroseq[sfrsegment]=int(vals[2])
        inv_hydroseq[int(vals[2])]=sfrsegment
        uphydroseq[sfrsegment]=int(vals[3])
        dnhydroseq[sfrsegment]=int(vals[4])
        cellhydros[cellnum].append(hydroseq[sfrsegment])

ordering.close()            

#now go through hydrowork and remove any keys with empty lists of cells
for rawsfrsegment in hydrowork:
    if not hydrowork[rawsfrsegment]:
        print 'raw sfr segment dropped', rawsfrsegment
        del hydrowork[rawsfrsegment]

#now sort the remaining entries of the dictionary in ascending order and
#then use the map command to convert to a list of list of cells
hydro_ordered=hydrowork.items()
hydro_ordered.sort(key=itemgetter(0))
hydro_orderedcells=map(itemgetter(1),hydro_ordered)  #list of the lists of cells
hydro_orderedindex=map(itemgetter(0),hydro_ordered)  #list of the indexes (ordered)

#generate list and dictionary of hydrosequences, in order corresponding
#to hydro_ordered cells
CHK2=open('check2.out','w')
SFRhydrodict=dict()
ordereduphydro=dict()
ordereddnhydro=dict()
iseg=1
for rawindex in hydro_orderedindex:
    SFRhydrodict[hydroseq[rawindex]]=iseg
    ordereduphydro[iseg]=uphydroseq[rawindex]
    ordereddnhydro[iseg]=dnhydroseq[rawindex]
    iseg=iseg+1

CHK2.write('hydroorderedcells\n')
for i in range(0,len(hydro_orderedcells)):
    lastcell=hydro_orderedcells[i][-1]
    tocell=fromCell[lastcell]
    if lastcell in tocell:
        tocell.remove(lastcell)
    match=0
    matchseg=[]
    for j in range(0,len(hydro_orderedindex)):
        for cell in hydro_orderedcells[j]:
            for connectcell in tocell:
                if cell==connectcell:
                    match=1
                    matchseg.append(cell)
    if match==0:
        CHK2.write(",".join(map(str,[i+1,lastcell,tocell]))+", no match\n")
    if match==1:
        CHK2.write(",".join(map(str,[i+1,lastcell,tocell,matchseg]))+"\n")

CHK2.write('up hydrosequence\n')
for k,v in ordereduphydro.items():
    CHK2.write(",".join(map(str,[k,v]))+'\n')
CHK2.write('down hydrosequence\n')
for k,v in ordereddnhydro.items():
    CHK2.write(",".join(map(str,[k,v]))+'\n')

CHK2.write('SFRHYDODICT\n')
for k,v in SFRhydrodict.items():
    CHK2.write(",".join(map(str,[k,v]))+'\n')

CHK2.write('HYDROSEQ\n')
for k,v in hydroseq.items():
    CHK2.write(",".join(map(str,[k,v]))+'\n')

CHK2.write('CELLS-HYDROSEQUENCE\n')
for k,v in cellhydros.iteritems():
    CHK2.write(",".join(map(str,[k,v]))+'\n')
    
CHK2.close()
del hydrowork, hydro_ordered
# have it all now by cellnumber... loop over lists by hydrosequence and generate SFR2 package
try:
    outfile=open(OUT,'w')
except:
    print "could not open output file"
    exit()
outfile.write(r'# SFR2 Package generated by python scrips from grid shapefile and NHDPlus data'+'\n')

#open a file that helps with GIS visualization of SFR 
gisout=open(GIS,'w')
gisout.write('cellnum,row,column,layer,segment,reach\n')

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

#count the segments and reaches- some get dropped because of short lengths...
nstrm=0
cellprint=dict()
iseen=dict()
nss=len(hydro_orderedcells)
for i in range(0,nss):
    for j in range(0,len(hydro_orderedcells[i])):
        if hydro_orderedcells[i][j] not in cellprint:
            nstrm=nstrm+1
            cellprint[hydro_orderedcells[i][j]]=1
          
cellprint=dict()
printstring=(-nstrm,nss,nsfrpar,nparseg,const,dleak,istcb1,istcb2,isfropt,nstrail,isuzn,nsfrsets)
outfile.write(' '.join(map(str,printstring))+'\n')
iseg=0
iseen=dict()
for i in range(0,len(hydro_orderedcells)):
    irch=0
    for j in range(0,len(hydro_orderedcells[i])):
        if hydro_orderedcells[i][j] in cellseen and hydro_orderedcells[i][j] not in cellprint:
            if not i in iseen:
                iseg=iseg+1
                iseen[i]=iseg
            irch=irch+1
            printstring=(int(1),  
                            row[hydro_orderedcells[i][j]],
                            column[hydro_orderedcells[i][j]],
                            iseg,
                            irch)
            outfile.write(' '.join(map(str,printstring)))
            #assume top of streambed is 1 ft below elevation from NHDPlus
            floatlist=[totlength[hydro_orderedcells[i][j]],
                        elevcell[hydro_orderedcells[i][j]]-1.0,
                        weightedslope[hydro_orderedcells[i][j]],
                        bedthick,
                        bedK]
            printstring=' {0:.2f} {1:.2f} {2:.5f} {3:.2f} {4:.2f}'.format(*floatlist)
            outfile.write(printstring + '\n')

            #write ouput to GWV matrix 1 file
            printstring=printstring=(  
                            row[hydro_orderedcells[i][j]],
                            column[hydro_orderedcells[i][j]],
                            int(1),)
            mat1out.write(','.join(map(str,printstring)))
            mixedlist=(elevcell[hydro_orderedcells[i][j]],
                       elevcell[hydro_orderedcells[i][j]]-1.0,
                       irch,
                       iseg,
                       weightwidth[hydro_orderedcells[i][j]],
                       totlength[hydro_orderedcells[i][j]],
                       bedK,
                       bedthick,
                       weightedslope[hydro_orderedcells[i][j]],
                       roughch)
            printstring=',{0:.2f},{1:.2f},{2:d},{3:d},{4:.2f},{5:.2f},{6:.2f},{7:.2f},{8:.5f},{9:.4f}'.format(*mixedlist)
            mat1out.write(printstring+'\n')
            
            gisstring=(hydro_orderedcells[i][j],
                       row[hydro_orderedcells[i][j]],
                       column[hydro_orderedcells[i][j]],
                       int(1),
                       iseg,
                       irch)
            gisout.write(','.join(map(str,gisstring))+'\n')
            cellprint[hydro_orderedcells[i][j]]=1

outfile.write('1 0 0 0\n')   # item 5...

for i in range(0,len(hydro_orderedcells)):
    if i in iseen:
        iseg=iseen[i]
        for j in range(0,len(hydro_orderedcells[i])):
            if hydro_orderedcells[i][j] in cellseen:
                begincell=j
                break
        for j in range(len(hydro_orderedcells[i])-1,-1,-1):
            if hydro_orderedcells[i][j] in cellseen:
                endingcell=j
                break
        #iupseg is for diversions so set it to zero
        iupseg=int(0)    
        if iseg in ordereddnhydro:
            dnhydro=ordereddnhydro[iseg]
            if dnhydro in SFRhydrodict:
                outseg=SFRhydrodict[dnhydro]
            else:
                outseg=int(0)
        else:
            outseg=int(0)
        if outseg==0:
            #one last check that this is really a terminal segment,
            #check the eight surrounding cells and see if any are
            #SFR cells, then see if any of these are values of the
            #dictionary->list  fromCell that gives the cells downstream of
            #the cell that is the key to the dictionary.  If there is
            #a match, then outseg should be the cell listed in fromCell[end of current segment]
            endcellnumber=hydro_orderedcells[i][endingcell]
            eightchecks=[endcellnumber-1,endcellnumber+1,endcellnumber+numcol,endcellnumber-numcol,
                         endcellnumber+numcol-1,endcellnumber+numcol+1,endcellnumber-numcol-1,endcellnumber-numcol+1]
            for checkcell in eightchecks:
                if checkcell in cellhydros:
                    for kk in range(0,len(fromCell[endcellnumber])):
                        if fromCell[endcellnumber][kk]==checkcell:
                            outseg=SFRhydrodict[cellhydros[checkcell][0]]
                        
        printstring=[iseg,icalc,outseg,iupseg]
    ##    if iupseg > 0:
    ##        printstring.append(iprior)
        if icalc==4:
            printstring.append(nstrpts)
        printstring.extend([runoff,etsw,pptsw])
        if icalc==1 or icalc==2:
            printstring.append(roughch)
        outfile.write(' '.join(map(str,printstring))+'\n')
        #for SFR2, K and thick are in previous line (reach-by-reach) only width is needed here...
        #printstring=[Hc1fact,thickm1]
        #printstring.append(elevcell[hydro_orderedcells[i][begincell]])
        printstring='{0:.1f}'.format(weightwidth[hydro_orderedcells[i][begincell]])  #width1...
        outfile.write(printstring+'\n')
        #printstring=[Hc2fact,thickm2]
        #printstring.append(elevcell[hydro_orderedcells[i][endingcell]])
        printstring='{0:.1f}'.format(weightwidth[hydro_orderedcells[i][endingcell]])  #width2...
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
gisout.close()
mat1out.close()
mat2out.close()


