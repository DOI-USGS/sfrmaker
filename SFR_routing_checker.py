import numpy as np
import matplotlib.pyplot as plt
import sys
import xml.etree.ElementTree as ET


        

try:
    import shapefile
    shapefiles_imported = True
except: 
    shapefiles_imported = False
    print ('\n'
       'WARNING: Unable to import the "shapefile" module. \n\n'
       'Shapefile generation has been disabled for this run, \n'
       'regardless of whether it was requested. \n'
       'The ascii output and PDF images will still be generated. \n\n')

    
class model_rc_conversion:
    def __init__(self,infile):
        self.spcfile = infile
        self.activeriv = 1.0
        self.upstream_routed = 2.0
        
    def rc2world_coords(self,row,col):
        xoff = self.xoffset + np.sum(self.deltax[0:col])+self.deltax[0]/2.0
        yoff = self.yoffset - np.sum(self.deltay[0:row])-self.deltay[0]/2.0
        x_low_left = self.xoffset
        y_low_left = self.yoffset - np.sum(self.deltay)
        xoff_to_rot = xoff - x_low_left
        yoff_to_rot = yoff - y_low_left
        cloc = np.array([[xoff_to_rot],[yoff_to_rot]])
        cloc_adj = np.dot(cloc.T,self.rot_matrix)[0]
        return cloc_adj[0] + x_low_left,cloc_adj[1]+y_low_left
        
    def spc_reader(self):
        # reads a PEST-compatible model.spc file
        indat = open(self.spcfile,'r').readlines()
        header = indat.pop(0).strip().split()
        self.nrows = int(header[0])
        self.ncols = int(header[1])
        coords = indat.pop(0).strip().split()
        self.xoffset = float(coords[0])
        self.yoffset = float(coords[1])
        self.mod_rot_deg = float(coords[2])
        self.mod_rot_rad = np.pi*self.mod_rot_deg/180.0
        th = self.mod_rot_rad # quick shorthand to make the rotation matrix
        self.rot_matrix = np.array([[np.cos(th), np.sin(th)],
                                    [-np.sin(th), np.cos(th)]])
        spacing = []
        for line in indat:
            if '*' not in line:
                spacing.extend(line.strip().split())
            else:
                tmp = line.strip().split()
                tmp2 = line.strip().split('*')
                print tmp
                print tmp2
                starlox = [i for i,x in enumerate(line.strip()) if x=='*']
                for i in np.arange(len(starlox)):
                    # check for the line balls on the left
                    if i==0:
                        mults = int(tmp2[i])
                    else:
                        mults = int(tmp2[i].split()[1])
                    # check for lienballs on the right
                    if i==len(starlox):
                        cval = float(tmp2[-1])
                    else:
                        cval = float(tmp2[i+1].split()[0])
                    for cn in np.arange(mults):
                        spacing.extend([cval])
        spacing = np.array(spacing,dtype=float)
        self.deltax = spacing[0:self.ncols]
        self.deltay = spacing[self.ncols:]
        
def tf2flag(intxt):
    # converts text written in XML file to True or Fale flag
    if intxt.lower()=='true':
        return True
    else:
        return False
# ######
#  MAIN
# ######

parfilename = sys.argv[1]
inpardat = ET.parse(parfilename)
inpars = inpardat.getroot()    
infile = inpars.findall('.//sfr_file')[0].text
modspecfile = inpars.findall('.//spc_file')[0].text
make_shapefiles = inpars.findall('.//make_shapefiles')[0].text
make_PDFs = inpars.findall('.//make_pdffiles')[0].text
make_all_layers = inpars.findall('.//all_layers')[0].text
make_shapefiles = tf2flag(make_shapefiles)
make_PDFs = tf2flag(make_PDFs)
make_all_layers = tf2flag(make_all_layers)
origseg = int(inpars.findall('.//orig_seg')[0].text)
origreach = int(inpars.findall('.//orig_reach')[0].text)


print 'Calculating routing from seg: {0:d} reach: {1:d}'.format(origseg, origreach)
model_spc_data= model_rc_conversion(modspecfile)

# read in the model.spc file
model_spc_data.spc_reader()

x,y = model_spc_data.rc2world_coords(450,500)

indat = open(infile,'r').readlines()

# first get rid of the first lines with comments
indat2 = []
for line in indat:
    if not '#' in line:
        indat2.append(line.strip().split())
del indat

numreaches = np.abs(int(indat2[0][0]))
numsegs = np.abs(int(indat2[0][1]))
reachdata = np.array(indat2[1:numreaches+1])[:,0:5].astype(int)
np.savetxt('reachdata.dat',reachdata,fmt='%12d')
# kludgey reading of the bottom part of the SFR file for segments
tmp = indat2[numreaches+2:]
del indat2
# only grab first stress period
segdata = np.array(tmp[::3][0:numsegs])[:,0:3].astype(int)
del tmp
np.savetxt('segdata.dat',segdata,fmt='%12d')

nrows = np.max(reachdata[:,1])
ncols = np.max(reachdata[:,2])
nlays = np.max(reachdata[:,0])


# now we have the reach data section and the segment data section
# reach data columns are:
# [0] layer  [1] row  [2] column  [3] segment  [4] reach route (1 is farthest upstream)
# segment data columns are:
# [0] segment  [1] icalc  [2] outseg (if outseg==0, routed out of model)



startind = np.where((reachdata[:,3]==origseg) &
                    (reachdata[:,4]==origreach))[0]
starting_reach = startind[0]
starting_seg = reachdata[starting_reach,3]

# find all the upstream parents of our starting segment
parents = np.zeros(0)
livesegs  = [starting_seg]
segcount = 0
while 1 and segcount < numsegs:
    segcount += 1
    tmp = np.zeros(0)
    for cseg in livesegs:
        upinds = np.where(segdata[:,2]==cseg)[0]
        circular = np.intersect1d(parents,segdata[upinds,0])
        if len(circular) > 0:
            print 'WARNING! -- possible circular routing!'
        parents = np.append(parents,segdata[upinds,0])
        tmp = np.append(tmp,segdata[upinds,0])
    livesegs = tmp
    if len(tmp) < 1:
        break
    del tmp
if segcount == numsegs:
    print 'WARNING! -- possible circular routing!'
np.savetxt('parents.dat',parents,fmt='%12d')

plotting = np.zeros((nlays,nrows,ncols))
allsegs = np.zeros_like(plotting)
allreaches = np.zeros_like(plotting)


# first burn in all the streams
for i in np.arange(nlays):
    inds = np.where(reachdata[:,0]==i+1)
    plotting[i,reachdata[inds,1]-1,reachdata[inds,2]-1] = model_spc_data.activeriv
    for cind in inds:
        allsegs[i,reachdata[cind,1]-1,reachdata[cind,2]-1] = reachdata[cind,3]
        allreaches[i,reachdata[cind,1]-1,reachdata[cind,2]-1] = reachdata[cind,4]
    # plot the upstream reaches within the current segment
    
    inds = np.where((reachdata[:,0]==i+1) & 
                    (reachdata[:,3]==origseg) & 
                    (reachdata[:,4]<origreach))[0]
    plotting[i,reachdata[inds,1]-1,reachdata[inds,2]-1] = model_spc_data.upstream_routed
    for cind in inds:
        allsegs[i,reachdata[cind,1]-1,reachdata[cind,2]-1] = reachdata[cind,3]
        allreaches[i,reachdata[cind,1]-1,reachdata[cind,2]-1] = reachdata[cind,4]

# plot up all full segments upstream from the evaluated segment
for cseg in parents:
    inds = np.where(reachdata[:,3]==cseg)
    tmp_reaches = np.atleast_2d(np.squeeze(reachdata[inds,:]))
    for i in np.arange(nlays):
        inds2 = np.where(tmp_reaches[:,0]==i+1)[0]
        if len(inds2) > 0:
            plotting[i,tmp_reaches[inds2,1]-1,tmp_reaches[inds2,2]-1] = model_spc_data.upstream_routed
            for cind in inds2:
                allsegs[i,reachdata[cind,1]-1,reachdata[cind,2]-1] = reachdata[cind,3]
                allreaches[i,reachdata[cind,1]-1,reachdata[cind,2]-1] = reachdata[cind,4]
            
if make_PDFs:    
    for i in np.arange(nlays):
        plt.figure()
        plt.hold=True
        plt.imshow(plotting[i,:,:],interpolation='nearest',cmap="jet")
        plt.xlabel('Columns')
        plt.ylabel('Rows')
        plt.title('Streams in Layer %d' %(i+1))
        plt.colorbar()
        plt.savefig('layer%d.pdf' %(i+1))

if make_shapefiles:
    if shapefiles_imported == False:
        print 'requested shapefile generation but shapefile package could not be imported'
    else:
        pshape_all=shapefile.Writer(shapefile.POINT)
        pshape_all.field('SFR_status')
        pshape_all.field('segment')
        pshape_all.field('reach')
        for clay in np.arange(nlays):
            print 'making shapefile for layer %d' %(clay+1)
            if make_all_layers:
                pshape = shapefile.Writer(shapefile.POINT) 
                pshape.field('SFR_status')
                pshape.field('segment')
                pshape.field('reach')
            
            inds = np.where(plotting[clay,:,:]>0)
            for cind,r in enumerate(inds[0]):
                c=inds[1][cind]
                plotx,ploty = model_spc_data.rc2world_coords(r,c)
#                plotx *= .3048
#                ploty *= .3048

                if make_all_layers:                
                    pshape.point(plotx,ploty)
                pshape_all.point(plotx,ploty)
                if plotting[clay,r,c]==2:
                    if make_all_layers: 
                        pshape.record('upstream',allsegs[clay,r,c],allreaches[clay,r,c])
                    pshape_all.record('upstream',allsegs[clay,r,c],allreaches[clay,r,c])
                    
                elif plotting[clay,r,c]==1:
                    if make_all_layers: 
                        pshape.record('active',allsegs[clay,r,c],allreaches[clay,r,c])
                    pshape_all.record('active',allsegs[clay,r,c],allreaches[clay,r,c])
                else:
                    print 'point should be empty!'
            if make_all_layers: 
                pshape.save('Layer%d' %(clay+1))
        pshape_all.save('All_Layers')