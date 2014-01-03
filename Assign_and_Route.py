# Assign_and_Route.py
# Description: Takes shapefile generated from NHDplus with start (X,Y), end(X,Y),
# and line length added and a cell-by-cell intersection of the NHDplus
#
# assigns elevations to stream segments based on their position along
# the NHDPlus comid
#
# identifies comid-modflow cell mapping
#             Same COMID in and out of a cell
#             One COMID in- One COMID out routed through the cell
#             Multiple COMIDs in a cell that are not all routed
#                 could be a confluence or braid
#                 could be different streams that don't mean in the cell
#
# input shapefiles must have line starts and ends added
#
# This version adds LevelPathID information to the ordered_cells.txt file
# and generalizes the next down function to look at either hydrosequence or
# levelpathID
#
# Requirements: os, sys, arcpy, defaultdict, math
#
# Author: H.W. Reeves; USGS Michigan Water Science Center
# Date 10/22/13
#
import os
import sys
import arcpy
from collections import defaultdict
import math
from operator import itemgetter

# get name of field (useful for case issues and appended field names in joined tables, etc)
def getfield(table,joinname):
    Fields=arcpy.ListFields(table)
    joinname=joinname.lower()
    for field in Fields:
        if joinname in field.name.lower():
            joinname=field.name
            break
    return(joinname)

# find next downstream sequence (either hydrosequence or levelpath)
# that is used in a model from the
# full sequence table - sometimes in developing the network for
# a model some sequences get dropped and this finds the next
# one or assigns a zero value indicating that the current segment is a
# downstream end.  Using a formally recursive function exceeded the
# python recursion depth, so an iterative method using a limit and
# while loop is used.
def nextdown(dwn,dnseq,seqseen):
    # dwn is the current downstream sequence (hydrosequence or levelpathID)
    # dnhydroseq is the full down sequence dictionary, keyed by seq
    # seqseen is a dictionary of sequences used in the model
    # the full dictionary is read from the NHDPlus, PlusFlowlineVAA table
    #
    # this function is usually called if dwn is not in the
    # seqseen dictionary
    limit=1000
    knt=0
    while knt<limit:
        if not dwn in seqseen:
            if not dwn in dnseq:
                nxtdwn=0
                return nxtdwn             #next one not in dnseq, call it the downstream end
            nxtdwn=dnseq[dwn]
            if nxtdwn==0:                 #found that it is downstream end
                return nxtdwn
            elif nxtdwn in seqseen:  #found an existing downstream sequence
                return nxtdwn
            else:
                knt=knt+1
                dwn=nxtdwn
        else:
            return dwn                    #dwn is in seqseen, so just return it
                                          #allows function to be called for any seq

    nxtdwn=0                              #fell through, assign it as a downstream end
    return nxtdwn


# HARD CODE INPUT SHAPEFILES

INTERSECT="river_explode.shp"
OUTTAB="river_elevs.dbf"
outfilename="fix_comids.txt"
# Get input files locations
# Global Input file for SFR utilities (see also for input instructions)
infile="SFR_setup.in"
infile=open(infile,'r').readlines()
inputs=defaultdict()
for line in infile:
    if "#" not in line.split()[0]:
        varname=line.split("=")[0]
        var=line.split("=")[1].split()[0]
        inputs[varname]=var.replace("'","") # strip extra quotes

ELEV=inputs["ELEV"] # 'river_w_elevations.shp'
FLOW=inputs["FLOW"] # 'PlusFlow_merge.dbf'
GRID=inputs["MFgrid"] # 'nacp_mg_cells_chk_node.shp'

#set elevation difference to assign same value to all cells
#so that check does not get messed up by rounding floats
eps=1.0e-02

#set workspace
path=os.getcwd()
arcpy.env.workspace=path
arcpy.env.overwriteOutput=True

#create output DBF table
if arcpy.Exists(OUTTAB):
    arcpy.Delete_management(OUTTAB)
arcpy.CreateTable_management(path,OUTTAB)
arcpy.AddField_management(OUTTAB,"OLDFID","LONG")
arcpy.AddField_management(OUTTAB,"CELLNUM","LONG")
arcpy.AddField_management(OUTTAB,"ELEVMAX","DOUBLE")
arcpy.AddField_management(OUTTAB,"ELEVAVE","DOUBLE")
arcpy.AddField_management(OUTTAB,"ELEVMIN","DOUBLE")
        
#loop over COMIDs explode shapefile and read information into
#hashes
startx=defaultdict(list)
starty=defaultdict(list)
endx=defaultdict(list)
endy=defaultdict(list)
maxsmoothelev=dict()
minsmoothelev=dict()
totallength=dict()
FID=defaultdict(list)
lengthft=defaultdict(defaultdict)
cellnum=defaultdict(defaultdict)
#make a list and some dictionaries to be keyed by cellnumber
cellnumlist=[]
FIDlist=defaultdict(list)
COMIDlist=defaultdict(list)
cellsforcomid=defaultdict(list)

#clean up river explode a bit
keepfields={"X_start":1,
            "Y_start":1,
            "COMID":1,
            "X_end":1,
            "Y_end":1,
            "FID":1,
            "LengthFt":1,
            "LENGTHKM":1,
            "MAXELEVSMO":1,
            "MINELEVSMO":1,
            "node":1,
            "row":1,
            "column":1,
            "delx":1,
            "dely":1,
            "ORIG_FID":1,}
allfields=arcpy.ListFields(INTERSECT)
dropfields=[]

for field in allfields:
    if not field.name in keepfields and not field.required:
        dropfields.append(field.name)
if len(dropfields)>0:
    arcpy.DeleteField_management(INTERSECT,dropfields)
segments=arcpy.SearchCursor(INTERSECT)
for seg in segments:
    comid=int(seg.COMID)
    startx[comid].append(float(seg.X_start))
    starty[comid].append(float(seg.Y_start))
    endx[comid].append(float(seg.X_end))
    endy[comid].append(float(seg.Y_end))
    FID[comid].append(int(seg.FID))
    lengthft[comid][int(seg.FID)]=float(seg.LengthFt)
    cellnum[comid][int(seg.FID)]=int(seg.node)
    lengthkmft=float(seg.LENGTHKM)*1000.0*3.28
    if(lengthkmft<0.01):
        print "check LENGTHKM of comid= ",comid
        exit()
    #the 04 smoothed elevations are in cm
    #convert to ft/  NACP input are in m, convert to ft
    maxsmoothelev[comid]=float(seg.MAXELEVSMO)*3.2808
    minsmoothelev[comid]=float(seg.MINELEVSMO)*3.2808
    cn=int(seg.node)  #cellnumber
    cellnumlist.append(cn)
    FIDlist[cn].append(int(seg.FID))
    COMIDlist[cn].append(comid)
    cellsforcomid[comid].append(cn)
        
del seg
del segments

#read through river_elev shapefile and get list of all COMIDs; then
#read through NHDFlow and pull the from/to information for those COMIDs
comids=dict()
elev=dict()
with arcpy.da.SearchCursor(INTERSECT,("COMID")) as cursor:
    for row in cursor:
        comids[int(row[0])]=1

print 'have comids from %s' % INTERSECT

tempfromcomid=[]
temptocomid=[]
with arcpy.da.SearchCursor(FLOW,("FROMCOMID","TOCOMID")) as cursor:
    for row in cursor:
        tempfromcomid.append(int(row[0]))
        temptocomid.append(int(row[1]))

print 'have the froms and tos from %s'% FLOW
    
fromCOMIDlist=defaultdict(list)
toCOMIDlist=defaultdict(list)
#fromCOMIDlist - key is the 'from', value is the 'to'
#toCOMIDlist - key is the 'to', value is the 'from'

#read in the file with information for clipped segments at the boundary
print "Processing clipped segments along domain boundary..."
CLIP=open("boundaryClipsRouting.txt",'r').readlines()
header = CLIP.pop(0) #header
clipto=dict()
clipfrom=dict()
for line in CLIP:
    vals=line.strip().split(',')
    fromin=int(float(vals[0]))
    toin=int(float(vals[1]))
    if fromin==99999:
        clipto[toin]=1
    if toin==99999:
        clipfrom[fromin]=1

for i in range(0,len(tempfromcomid)):
    fromin=tempfromcomid[i]
    toin=temptocomid[i]
    if not toin in clipto:
        toCOMIDlist[toin].append(fromin)
    if fromin in comids:
        if not fromin in clipfrom:
            if toin in comids:
                fromCOMIDlist[fromin].append(toin)
            else:
                fromCOMIDlist[fromin].append(99999)
        else:
            fromCOMIDlist[fromin].append(99999)
#now have full from->to table as a dictionary of lists
#go through river_elev file COMID by COMID and check
#that next downstream COMID has lower river elevation

            
#............old stuff
fact=0.95
segcounter=0
segelevinfo=defaultdict(defaultdict)
rows=arcpy.InsertCursor(OUTTAB)
noelev=dict()
COMID_orderedFID=dict()
fix_comids_summary = []
fixcomids_flag = False
for comid in FID.iterkeys():
    segcounter=segcounter+1
    #print "COMID = ",comid," segment is ",segcounter
    fidlist=FID[comid]
    start_has_end=dict()
    end_has_start=dict()

    for i in range(0,len(fidlist)):
        haveend=False
        havestart=False
        for j in range(0,len(fidlist)):
            if j==i:
                continue
            diffstartx=startx[comid][i]-endx[comid][j]
            diffstarty=starty[comid][i]-endy[comid][j]
            diffendx=endx[comid][i]-startx[comid][j]
            diffendy=endy[comid][i]-starty[comid][j]
            if(math.fabs(diffstartx)<fact and math.fabs(diffstarty)<fact):
                start_has_end[fidlist[i]]=fidlist[j]
                haveend=True
            if(math.fabs(diffendx)<fact and math.fabs(diffendy)<fact):
                end_has_start[fidlist[i]]=fidlist[j]
                havestart=True
            if haveend and havestart:
                break
        
    #find key in start_has_end that didn't match an end and
    #key in end_has_start that didn't match a start
        
    numstart=0
    numend=0
    startingFID=[]
    endingFID=[]
    for test in fidlist:
        if test not in start_has_end:
            startingFID.append(test)
            numstart=numstart+1
        if test not in end_has_start:
            endingFID.append(test)
            numend=numend+1
    if (numstart!=1 or numend !=1):
        if fixcomids_flag == False:
            outfile = open(outfilename,'w')
            fixcomids_flag = True
        outfile.write("numstart ="+ str(numstart)+" \n")
        outfile.write("numend = "+ str(numend)+" \n")
        outfile.write("starting FIDs: " + ",".join(map(str,startingFID))+"\n")
        outfile.write("ending FIDs: " + ",".join(map(str,endingFID))+"\n")
        outfile.write("manually fix COMID = %d\n" %comid)
        fix_comids_summary.append('%d\n' %comid)
        noelev[comid]=1  #set flag 
        continue

    orderedFID=[]
    orderedFID.append(startingFID[0])
    for i in range(1,len(end_has_start)):
        orderedFID.append(end_has_start[orderedFID[i-1]])
    if orderedFID[-1]!=endingFID[0]:       #don't repeat the last entry FID...
        orderedFID.append(endingFID[0])
    #total length read through lengthkm didn't always match up
    #to the sum of the lengths of the segments (exactly), sumup total length
    totallength=0
    for i in range(0,len(orderedFID)):
        totallength=totallength+lengthft[comid][orderedFID[i]]
    if totallength==0:
        exit('check length ft for FIDs in COMID= %d' % comid)
    slope=(maxsmoothelev[comid]-minsmoothelev[comid])/totallength
    distance=0.
    COMID_orderedFID[comid]=orderedFID
    for i in range(0,len(orderedFID)):
        maxcellrivelev=maxsmoothelev[comid]-slope*distance
        distance=distance+lengthft[comid][orderedFID[i]]
        mincellrivelev=maxsmoothelev[comid]-slope*distance
        avecellrivelev=0.5*(maxcellrivelev+mincellrivelev)
        segelevinfo[comid][cellnum[comid][orderedFID[i]]]=avecellrivelev
        row=rows.newRow()
        row.OLDFID=orderedFID[i]
        row.CELLNUM=cellnum[comid][orderedFID[i]]
        row.ELEVMAX=maxcellrivelev
        row.ELEVAVE=avecellrivelev
        row.ELEVMIN=mincellrivelev
        rows.insertRow(row)
# write out the summary of comids to fix, then close the outfile
if len(fix_comids_summary) > 0:
    print 'Some cells have multiple COMIDs entering and/or leaving.\n See file "fix_comids.txt"'
    outfile.write('#' * 30 + '\nSummary of COMIDS to fix:\n' +
    'Delete these COMIDs from river_explode.shp, \nthen run CleanupRiverCells.py and rerun Assign_and_Route.py\n')
    [outfile.write(line) for line in fix_comids_summary]
    outfile.close()
del row
del rows

#....................new stuff
uniqcells=list(set(cellnumlist))
del cellnumlist

print "shapefile read in"

cellcomidin=defaultdict(list)
cellcomidout=defaultdict(list)
connecttype=dict()

OUT=open('checkrevised.out','w')
for cn in uniqcells:
    #if number of comids in the cell > 1; see if they are connected..
    if len(COMIDlist[cn])>1:
        connected=0
        uniqs=list(set(COMIDlist[cn]))
        yesupstream=[]
        for i in range(0,len(uniqs)):
            comid=uniqs[i]
            tempset=set(uniqs)
            tempset.remove(comid)
            checkset=set(fromCOMIDlist[comid])
            #see if there is a match in tempset and checkset
            matches=list(set(tempset).intersection(checkset))
            if len(matches)>0:
                OUT.write(str(cn)+' There is a downstream comid in the cell, up '+str(comid)+' :down '+",".join(map(str,matches))+'\n')
                for j in range(0,len(matches)):
                    yesupstream.append(matches[j])
            else:
                OUT.write(str(cn)+' no downstream comid in the cell, up '+str(comid)+'\n')
                cellcomidout[cn].append(comid)
        ups=set(yesupstream)
        uniqs=set(COMIDlist[cn])
        noupstream=list(uniqs.difference(ups))
        if len(noupstream)>0:
            OUT.write(str(cn)+' no upstream comid in the cell, '+",".join(map(str,noupstream))+'\n\n')
            cellcomidin[cn].extend(noupstream)
        else:
            OUT.write('\n')
    else:
        OUT.write('only 1 comid in cell: {0:d}\n\n'.format(cn))
        cellcomidin[cn].append(COMIDlist[cn][0])
        cellcomidout[cn].append(COMIDlist[cn][0])
        connecttype[cn]='same'
OUT.close()

cellelevlist=defaultdict(list)

for cn in uniqcells:
    if not cn in connecttype:
        if len(cellcomidin[cn])==1 and len(cellcomidout[cn])==1:
            connecttype[cn]='onein/oneout'
            if cellcomidin[cn][0] in noelev or cellcomidout[cn][0] in noelev:
                cellelevlist[cn]=[999999,999999]
            else:
                cellelevlist[cn]=[segelevinfo[cellcomidin[cn][0]][cn],segelevinfo[cellcomidout[cn][0]][cn]]
        else:
            connecttype[cn]='multiple'
            for i in range(0,len(cellcomidin[cn])):
                if cellcomidin[cn][i] in noelev:
                    cellelevlist[cn].append(99999)
                else:
                    cellelevlist[cn].append(segelevinfo[cellcomidin[cn][i]][cn])
            for i in range(0,len(cellcomidout[cn])):
                if cellcomidout[cn][i] in noelev:
                    cellelevlist[cn].append(99999)
                else:
                    cellelevlist[cn].append(segelevinfo[cellcomidout[cn][i]][cn])
    else:
        comid=cellcomidin[cn][0]
        if not cn in segelevinfo[comid]:
            cellelevlist[cn].append(99999)
        elif comid in noelev:
            cellelevlist[cn].append(99999)
        else:
            cellelevlist[cn].append(segelevinfo[comid][cn])
            

OUTTAB='cell_in_out.dbf'
if arcpy.Exists(OUTTAB):
    arcpy.Delete_management(OUTTAB)
arcpy.CreateTable_management(path,OUTTAB)
arcpy.AddField_management(OUTTAB,"CELLNUM","LONG")
arcpy.AddField_management(OUTTAB,"TYPE","TEXT")
arcpy.AddField_management(OUTTAB,"COUNT","SHORT")
arcpy.AddField_management(OUTTAB,"COMIDIN","TEXT")
arcpy.AddField_management(OUTTAB,"COMIDOUT","TEXT")
arcpy.AddField_management(OUTTAB,"ELEVS","TEXT")
arcpy.AddField_management(OUTTAB,"AVE_ELV","FLOAT")
rows=arcpy.InsertCursor(OUTTAB)
#sort the connecttype dictionary by values.. but get the keys
for (cn,v) in sorted(connecttype.iteritems(), key=itemgetter(1)):   
    row=rows.newRow()
    row.CELLNUM=cn
    row.TYPE=connecttype[cn]
    row.COUNT=len(COMIDlist[cn])
    row.COMIDIN=",".join(map(str,cellcomidin[cn]))
    row.COMIDOUT=",".join(map(str,cellcomidout[cn]))
    row.ELEVS=",".join(map(str,cellelevlist[cn]))
    summ=0
    count=0
    for i in range(0,len(cellelevlist[cn])):
        summ=summ+cellelevlist[cn][i]
        count=count+1
    row.AVE_ELV=summ/count    
    rows.insertRow(row)
del row, rows

#join with grid
if arcpy.Exists("grid_temp"):
    arcpy.Delete_management("grid_temp")
if arcpy.Exists("cell_inout.shp"):
    arcpy.Delete_management("cell_inout.shp")
#use the field names in the join, not table.fieldname
#the longer version gets trucated for the shapefile...
arcpy.env.qualifiedFieldNames=False
#do the join and only keep the cells that get joined, save as shapefile
arcpy.MakeFeatureLayer_management("nacp_mg_cells_chk_node.shp","grid_temp")
arcpy.AddJoin_management("grid_temp","node",OUTTAB,"CELLNUM","KEEP_COMMON")
arcpy.CopyFeatures_management("grid_temp","cell_inout.shp")

#drop fields not needed
keepfields={"CELLNUM":1,
            "TYPE":1,
            "COUNT":1,
            "COMIDIN":1,
            "COMIDOUT":1,
            "ELEVS":1,
            "AVE_ELV":1}
allfields=arcpy.ListFields("cell_inout.shp")
dropfields=[]
for field in allfields:
    if not field.name in keepfields and not field.required:
        dropfields.append(field.name)
arcpy.DeleteField_management("cell_inout.shp",dropfields)

#drop the dbf file
arcpy.Delete_management(OUTTAB)

arcpy.RefreshCatalog(path)

RIVER_ELEVS="river_elevs.dbf"
if arcpy.Exists(ELEV):
    arcpy.Delete_management(ELEV)
if arcpy.Exists("grid_temp"):
    arcpy.Delete_management("grid_temp")
print "\nJoining %s to %s; \nsaving as %s (might take a minute)" %(RIVER_ELEVS,"river_explode.shp",ELEV)
arcpy.MakeFeatureLayer_management("river_explode.shp","grid_temp")
arcpy.AddJoin_management("grid_temp","FID",RIVER_ELEVS,"OLDFID")
arcpy.CopyFeatures_management("grid_temp",ELEV)

#read through river_elev shapefile and get list of all COMIDs; then
#read through NHDFlow and pull the from/to information for those COMIDs
comids=dict()
elev=dict()
with arcpy.da.SearchCursor(ELEV,("COMID","ELEVAVE")) as cursor:
    for row in cursor:
        comids[int(row[0])]=1
        elev[int(row[0])]=float(row[1])

print 'have comids and elevations from %s' % ELEV
print "now checking COMID elevation routing"

seencomid=dict()
flag=0
outfile=open("Flagged_comids.txt",'w')
outfile.write("COMID,DWNCOMID,UPELEV,DWNELEV\n")
for comid in comids.iterkeys():
    if comid in seencomid:
        continue
    if comid not in fromCOMIDlist:
        outfile.write('%d,99999,99999,99999\n' % comid)
    else:
        for dwncomid in fromCOMIDlist[comid]:
            if dwncomid!=99999:
                if elev[dwncomid] > (elev[comid]+0.01):
                    outfile.write("%d,%d,%f,%f\n" %(comid,dwncomid,elev[comid],elev[dwncomid]))
                    flag=1
                seencomid[comid]=1
outfile.close()
if flag==1:
    print "\nCHECK the Flagged_comids.txt file -> at least 1 downstream COMID"
    print "had a greater elevation than its upstream COMID\n"

print "\nNow routing the river cells"

CELLS=inputs["CELLS"]
NHD=inputs["NHD"]  # from intersect.py
VAA=inputs["PlusflowVAA"]

try:
    OUT=open("routed_cells.txt",'w')
except:
    print "could not open output file"
    exit()

OUT.write("fromcell, tocell\n")

#make dictionaris of hydroseq, uphydroseq, dnhydroseq numbers
#hydroseq by COMID, up and dwn hydroseq keyed by hydroseq
hydroseq=dict()
uphydroseq=dict()
dnhydroseq=dict()
levelpathID=dict()
uplevelpath=dict()
dnlevelpath=dict()
with arcpy.da.SearchCursor(VAA,("ComID","Hydroseq","UpHydroseq","DnHydroseq","LevelPathI","UpLevelPat","DnLevelPat")) as cursor:
    for row in cursor:
        comid=int(row[0])
        hydroseq[comid]=int(row[1])
        uphydroseq[int(row[1])]=int(row[2])
        dnhydroseq[int(row[1])]=int(row[3])
        levelpathID[comid]=int(row[4])
        uplevelpathin=int(row[5])
        if uplevelpathin!=levelpathID[comid]:
            uplevelpath[levelpathID[comid]]=uplevelpathin
        dnlevelpathin=int(row[6])
        if dnlevelpathin!=levelpathID[comid]:
            dnlevelpath[levelpathID[comid]]=dnlevelpathin
print 'have dictionary'

#delete any working layers if found
if arcpy.Exists("temp_lyr"):
    arcpy.Delete_management("temp_lyr")
if arcpy.Exists("temp_nhd_lyr"):
    arcpy.Delete_management("temp_nhd_lyr")

#now read in segment level information to order cells for each comid
#using logic from AssignRiverElev.py
count=0
RCH=open('reach_ordering.txt','w')
RCH.write('CELLNUM,COMID, hydroseq, uphydroseq, dnhydroseq, levelpathID, uplevelpath, dnlevelpath, SFRseq, localseq\n')

SFRseq=0
hydroseqused=dict()
levelpathused=dict()
#get a list of the comids/hydroseq in this model
for comid in fromCOMIDlist.iterkeys():
    if comid in hydroseq:
        hydroseqused[comid]=hydroseq[comid]
        levelpathused[comid]=levelpathID[comid]
    else:
        exit('check VAA table for COMID %d' % comid)

del hydroseq    
#sort from high hydroseq to low hydroseq and then
#use the sorted list to loop over comids

sortedtuples=hydroseqused.items()
sortedtuples.sort(key=itemgetter(1), reverse=True)
sortedlist=map(itemgetter(0),sortedtuples)

icount=0
for comid in sortedlist:
    if comid in noelev:
        continue
    orderedFID=COMID_orderedFID[comid]
    for i in range(0,len(orderedFID)-1):
        fromcell=cellnum[comid][orderedFID[i]]
        tocell=cellnum[comid][orderedFID[i+1]]
        if tocell == fromcell:
            continue
        else:
            OUT.write("%d, %d\n" % (fromcell, tocell))
    dwncomid=fromCOMIDlist[comid][0]
    if orderedFID[-1] in cellnum[comid]:
        if orderedFID[0] in cellnum[dwncomid]:
            OUT.write("%d, %d\n" % (fromcell, tocell))

#build a dictionary of the hydroseqences that are used
#to help find next downstream if needed
hydroseqseen=dict()
for comid,seen in hydroseqused.iteritems():
    #if the comid is in noelev, then it gets skipped below, don't count it as seen
    if not comid in noelev:
        hydroseqseen[seen]=1
levelpathseen=dict()
for comid,seen in levelpathused.iteritems():
    #if the comid is in noelev, it gets skipped (same as hydroseq)
    if not comid in noelev:
        levelpathseen[seen]=1

for comid in sortedlist:
    if comid in noelev:
        continue
    icount+=1
    SFRseq=SFRseq+1
    #print SFRseq
    count=count+1
    orderedFID=COMID_orderedFID[comid]
    #now the COMID is ordered 
    #if cellnum[comid][orderedFID[0]]==cellnum[comid][orderedFID[-1]]:
        #comid starts and ends in same cell, doesn't help with cell routing, skip
        #print "starts and ends in same cell %d" % comid
        #continue    
    if(segelevinfo[comid][cellnum[comid][orderedFID[0]]]< segelevinfo[comid][cellnum[comid][orderedFID[-1]]]):
        orderedFID.reverse()

    for i in range(0,len(orderedFID)):
        #sometimes a hydrosequence from the full VAA table gets dropped
        #use a function to find the next downstream hydrosequence number
        #used in the model
        dwn=dnhydroseq[hydroseqused[comid]]
        if not dwn in hydroseqseen:
            dnhydroseq[hydroseqused[comid]]=nextdown(dwn,dnhydroseq,hydroseqseen)
        #same check for the downstream levelpathID
        dwn=dnlevelpath[levelpathID[comid]]
        if not dwn in levelpathseen:
            dnlevelpath[levelpathID[comid]]=nextdown(dwn,dnlevelpath,levelpathseen)
        #now write the information to reach_ordering file
        RCH.write(",".join(map(str,[cellnum[comid][orderedFID[i]],comid,
                                    hydroseqused[comid],uphydroseq[hydroseqused[comid]],dnhydroseq[hydroseqused[comid]],
                                    levelpathID[comid],uplevelpath[levelpathID[comid]],dnlevelpath[levelpathID[comid]],
                                    SFRseq,i+1]))+"\n")
    percentdone=round(100*icount/len(sortedlist),2)
    #print "%s %% done" %(percentdone)
print 'cell routing done, reach_ordering.txt written' 
OUT.close()
RCH.close()
arcpy.RefreshCatalog(path)


    

