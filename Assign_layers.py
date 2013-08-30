# Program to assign layers to SFR cells based on top of streambed - streambed thickness

import numpy as np
from collections import defaultdict

# Input files
botsfile='Columbia_bots_corr2.DAT' # GWV mat with bottom elevations for all layers
L1top='L1TOP.DAT' # GWV mat with top elevations for layer 1
SFRmat1='SFR_GWVmat1.txt' # SFR matrix 1 from Howard Reeve's scripts

# Settings
thick=3 # Uniform streambed thickness
buff=0 # if STOP - thick - buffer < bottom of cell, move to layer below

# Outputfiles
SFRout=SFRmat1[:-4]+'_Layers_assigned.DAT'
Layerinfo='SFR_layer_assignments.txt'

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
New_Layers=[]
for i in range(len(SFRinfo)):
    r=int(SFRinfo['row'][i])
    c=int(SFRinfo['column'][i])
    l=int(SFRinfo['layer'][i])
    STOP=float(SFRinfo['top_streambed'][i])
    cellbottoms=list(bots_rs[:,r-1,c-1])
    for b in range(nlayers):
        if (STOP-thick-buff)<cellbottoms[b]:
            if b+1 < nlayers:
                continue
            else:
                print 'Streambottom elevation=%s, Model bottom=%s at row %s, column %s, cellnum %s' %(STOP-thick-buff,cellbottoms[-1],r,c,(r-1)*ncols+c)
                print 'Land surface is %s' %(topdata[r-1,c-1])
                below_bottom.write('%s,%s,%s,%s,%s\n' %(STOP-thick-buff,cellbottoms[-1],topdata[r-1,c-1],(r-1)*ncols+c,SFRinfo['segment'][i]))
        else:
            New_Layers.append(b+1)
            break
below_bottom.close()
    
New_Layers=np.array(New_Layers)

# histogram of layer assignments
freq=np.histogram(list(New_Layers),range(nlayers+2)[1:])

print "Writing output to %s..." %(SFRout)
ofp=open(SFRout,'w')
ofp.write(','.join(SFRinfo.dtype.names)+'\n')

for i in range(len(SFRinfo)):
    line=list(SFRinfo[i])[0:2]+[New_Layers[i]]+list(SFRinfo[i])[3:]
    line=' '.join(map(str,line))
    ofp.write(line+'\n')
ofp.close()

ofp=open(Layerinfo,'w')
ofp.write('Layer\t\tNumber of assigned reaches\n')
print '\nLayer assignments:'
for i in range(nlayers):
    ofp.write('%s\t\t%s\n' %(freq[1][i],freq[0][i]))
    print '%s\t\t%s\n' %(freq[1][i],freq[0][i])
ofp.close()
print "\nsee also %s"
print "Done!"