# Program to assign layers to SFR cells based on top of streambed - streambed thickness

import numpy as np
from collections import defaultdict

# hard coded model dimensions
rows=800
columns=800 
# layers will be based on length of input GWV matrix and number of cells

# Input files
botsfile='BR_L1L5bot.DAT' # GWV mat with bottom elevations for all layers
l1topfile='L1top.DAT' # GWV mat with top elevations for layer 1
SFRmat1='BR_GWVmat1.DAT' # SFR matrix 1 from Howard Reeve's scripts

# Settings
thick=3 # Uniform streambed thickness
buff=0

# Outputfiles
SFRout=SFRmat1[:-4]+'_corr.DAT'
Layerinfo='SFR_layer_assignments.txt'

# load in GWV matricies
bots=np.fromfile(botsfile,sep=" ")
layers=len(bots)/(rows*columns)
bots_rs=bots.reshape(layers,rows,columns) # apparently np convention is l,r,c
SFRinfo=np.genfromtxt(SFRmat1,skiprows=1,dtype=None)


New_Layers=[]
for i in range(len(SFRinfo)):
    
    r=SFRinfo[i][0]
    c=SFRinfo[i][1]
    l=SFRinfo[i][2]
    STOP=SFRinfo[i][4]
    cellbottoms=list(bots_rs[:,r-1,c-1])
    
    for b in range(layers):
        if (STOP-thick-buff)<cellbottoms[b]:
            continue
        else:
            New_Layers.append(b+1)
            break
        
New_Layers=np.array(New_Layers)

# histogram of layer assignments
freq=np.histogram(list(New_Layers),range(layers+2)[1:])


ofp=open(SFRout,'w')
ofp.write('row,column,layer,stage,top_streambed,reach,segment,width_in_cell,length_in_cell,bed_K,bed_thickness,bed_slope,bed_roughness\n')

for lines in range(len(SFRinfo)):
    line=list(SFRinfo[lines])[0:2]+[New_Layers[lines]]+list(SFRinfo[lines])[3:]
    line=' '.join(map(str,line))
    ofp.write(line+'\n')
ofp.close()

ofp=open(Layerinfo,'w')
ofp.write('Layer\t\tNumber of assigned reaches\n')
for i in range(layers):
    ofp.write('%s\t\t%s\n' %(freq[1][i],freq[0][i]))
ofp.close()