# Program to assign layers to SFR cells based on top of streambed - streambed thickness

import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Input files
botsfile='BadRiver_bots_corr.dat' # GWV mat with bottom elevations for all layers
L1top='L1top.DAT' # GWV mat with top elevations for layer 1
SFRmat1='BR_GWVmat1.csv' # SFR matrix 1 from Howard Reeve's scripts
SFRmat2='BR_GWVmat2.csv'

# Settings
thick=3 # Uniform streambed thickness
buff=1 # if STOP - thick - buffer < bottom of cell, move to layer below
Lowerbot=False # if True lowers model bottom where SFR drops below
minimum_slope=1e-6

# Outputfiles
SFRpackage='Columbia.SFR'
newSFRmat1=SFRmat1[:-4]+'_Layers_assigned.csv'
newSFRmat2=SFRmat2[:-4]+'_updated_segnums.csv'
Layerinfo='SFR_layer_assignments.txt'
BotcorPDF='Model_bottom_adjustments.pdf' # PDF showing model bottom before and after adjustments if Lowerbot==True

print "Getting model grid elevations..."
temp=open(L1top).readlines()
ncols=len(temp[0].strip().split())
nrows=len(temp)

topdata=np.fromfile(L1top,sep=' ')
topdata=np.reshape(topdata,(nrows,ncols))
    
bots=np.fromfile(botsfile,sep=' ')
nlayers=int(len(bots)/np.size(topdata))
bots_rs=bots.reshape(nlayers,nrows,ncols) # apparently np convention is l,r,c
SFRinfo=np.genfromtxt(SFRmat1,delimiter=',',names=True,dtype=None)

print "Assigning layers to SFR cells..."
below_bottom=open('below_bot.csv','w')
below_bottom.write('SFRbot,ModelBot,Land_surf,cellnum,segment\n')
below_botdepths=defaultdict() # list row,column locations where SFR goes below model bottom
nbelow=0
New_Layers=[]
for i in range(len(SFRinfo)):
    r=int(SFRinfo['row'][i])
    c=int(SFRinfo['column'][i])
    l=int(SFRinfo['layer'][i])
    STOP=float(SFRinfo['top_streambed'][i])
    cellbottoms=list(bots_rs[:,r-1,c-1])
    for b in range(nlayers):
        SFRbot=STOP-thick-buff
        if (SFRbot)<cellbottoms[b]:
            if b+1 < nlayers:
                continue
            else:
                print 'Streambottom elevation=%s, Model bottom=%s at row %s, column %s, cellnum %s' %(SFRbot,cellbottoms[-1],r,c,(r-1)*ncols+c)
                print 'Land surface is %s' %(topdata[r-1,c-1])
                below_bottom.write('%s,%s,%s,%s,%s\n' %(SFRbot,cellbottoms[-1],topdata[r-1,c-1],(r-1)*ncols+c,SFRinfo['segment'][i]))
                below_botdepths[(r-1,c-1)]=cellbottoms[-1]-SFRbot # difference between SFR bottom and model bottom
                nbelow+=1
                New_Layers.append(b+1)
        else:
            New_Layers.append(b+1)
            break
below_bottom.close()
    
New_Layers=np.array(New_Layers)

botsnew=np.ndarray(shape=np.shape(topdata))
if Lowerbot:
    print "\n\nAdjusting model bottom to accomdate SFR cells that were below..."
    print "see %s\n" %(BotcorPDF)
    for r in range(nrows):
        for c in range(ncols):
            if (r,c) in below_botdepths.keys():
                print bots_rs[-1,r,c]
                botsnew[r,c]=bots_rs[-1,r,c]-below_botdepths[(r,c)]
                print botsnew[r,c]
            else:
                botsnew[r,c]=bots_rs[-1,r,c]

    outarray=botsfile[:-4]+'_SFRcorr.dat'
    bots_rs=np.append(bots_rs[:-1,:,:],botsnew)
    bots_rs=np.reshape(bots_rs,(nlayers,nrows,ncols))
    
    with file(outarray, 'w') as outfile:
        for slice_2d in bots_rs:
            np.savetxt(outfile,slice_2d,fmt='%.2f')
    outfile.close()
    
    # show bottom before and after corrections
    outpdf=PdfPages(BotcorPDF)
    for mat in [bots_rs[-1,:,:],botsnew]:
        plt.figure()
        plt.imshow(mat)
        plt.colorbar()
        outpdf.savefig()
    outpdf.close()

# histogram of layer assignments
freq=np.histogram(list(New_Layers),range(nlayers+2)[1:])

# Before writing output, enforce segment number continuity
# (get rid of any numbering gaps that may have resulted from deleting segments)
n=0
current_segnum=1
newsegment_numbers=defaultdict(list)
for segnum in SFRinfo['segment']:
    if segnum==current_segnum:
        newsegment_numbers[segnum].append(current_segnum-n)
        continue
    else:
        if segnum==current_segnum+1:
            newsegment_numbers[segnum]=[segnum-n]
            current_segnum=segnum
        else:
            newsegment_numbers[segnum]=[current_segnum+1-n]
            current_segnum=segnum
            n+=1

# flatten the dictionary of new segment numbers into single list; replace in SFRinfo
newsegnums=[item for sublist in sorted(newsegment_numbers.values()) for item in sublist]
SFRinfo['segment']=np.array(newsegnums)

# write new SFRmat1 file
print "\nWriting output to %s..." %(newSFRmat1)
ofp=open(newSFRmat1,'w')
ofp.write(','.join(SFRinfo.dtype.names)+'\n')

for i in range(len(SFRinfo)):
    line=list(SFRinfo[i])[0:2]+[New_Layers[i]]+list(SFRinfo[i])[3:]
    line=','.join(map(str,line))
    ofp.write(line+'\n')
ofp.close()


# update SFRmat2 (with new segment numbers)
newsegment_numbers[0]=[0] # add zero in (for terminal segments)
SFRinfo2=np.genfromtxt(SFRmat2,names=True,delimiter=',',dtype=None)
newinsegs=[]
newoutsegs=[]
for line in SFRinfo2:
    oldinseg=line[0]
    oldoutseg=line[2]
    newinsegs.append(newsegment_numbers[oldinseg][0])
    newoutsegs.append(newsegment_numbers[oldoutseg][0])
SFRinfo2['segment']=np.array(newinsegs)
SFRinfo2['outseg']=np.array(newoutsegs)
# write new SFRmat2 file
ofp=open(newSFRmat2,'w')
ofp.write(','.join(SFRinfo2.dtype.names)+'\n')
for line in SFRinfo2:
    ofp.write(','.join(map(str,line))+'\n')
ofp.close()


# writeout results to SFR package file
icalc=1
nreaches=len(SFRinfo)
nseg=np.max(SFRinfo['segment'])
ofp=open(SFRpackage,'w')
ofp.write("#SFRpackage file generated by Assign_Layers.py\n")
ofp.write('%s %s 0 0 128390.4 0.0001 50 53 1 0 0 0\n' %(-1*nreaches,nseg,))
for i in range(len(SFRinfo)):
    slope=SFRinfo['bed_slope'][i]
    if slope<=minimum_slope: # one last check for negative or zero slopes
        slope=minimum_slope
    ofp.write('%s %s %s %s %s %s %s %s %s %s\n' %(New_Layers[i],SFRinfo['row'][i],SFRinfo['column'][i],SFRinfo['segment'][i],SFRinfo['reach'][i],SFRinfo['length_in_cell'][i],SFRinfo['top_streambed'][i],slope,thick,SFRinfo['bed_K'][i]))
ofp.write('%s 0 0 0\n' %(nseg))
for i in range(len(SFRinfo2)):
    segment=SFRinfo2['segment'][i]
    seg_Mat1_inds=np.where(SFRinfo['segment']==segment)
    seginfo=SFRinfo[seg_Mat1_inds]
    ofp.write('%s %s %s 0 %s 0.0 0.0 0.0 %s\n' %(segment,icalc,SFRinfo2['outseg'][i],SFRinfo2['flow'][i],SFRinfo2['roughch'][i]))
    ofp.write('%s\n' %(seginfo['width_in_cell'][0]))
    ofp.write('%s\n' %(seginfo['width_in_cell'][-1]))
ofp.close()
    
# writeout info on Layer assignments
ofp=open(Layerinfo,'w')
ofp.write('Layer\t\tNumber of assigned reaches\n')
print '\nLayer assignments:'
for i in range(nlayers):
    ofp.write('%s\t\t%s\n' %(freq[1][i],freq[0][i]))
    print '%s\t\t%s\n' %(freq[1][i],freq[0][i])
ofp.close()

if not Lowerbot:
    if nbelow>0:
        print "Warning %s SFR streambed bottoms were below the model bottom. See below_bots.csv" %(nbelow)
print "Done!"
