# brings in Layer1 top elevations and computes differences with streambed top

import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def stopcomp(L1top,MAT1,outfile):
    # L1top= ascii matrix export of MODFLOW top elevations (n columns x n rows; no wrapping!)
    # Mat1= output from HWR's SFR_utilities script with SFR reach information
    
    print "Getting model grid information..."
    temp=open(L1top).readlines()
    ncols=len(temp[0].strip().split())
    nrows=len(temp)
    
    topdata=np.fromfile(L1top,sep=' ')
    topdata=np.reshape(topdata,(nrows,ncols))
    
    L1tops=dict()
    for r in range(nrows):
        for c in range(ncols):
            cellnum=r*ncols+c+1
            L1tops[cellnum]=topdata[r,c]
    
    
    MAT1data=np.genfromtxt(MAT1, delimiter=',',names=True,dtype=None)
    infile_data=defaultdict(list)
    for i in range(len(MAT1data)):
        row,column=MAT1data['row'][i],MAT1data['column'][i]
        cellnum=(row-1)*ncols+column
        infile_data['MSEG'].append(MAT1data['segment'][i])
        infile_data['STOP'].append(MAT1data['top_streambed'][i])
        infile_data['TOPNEW'].append(L1tops[cellnum])
        infile_data['DIF'].append(L1tops[cellnum]-MAT1data['top_streambed'][i])
    diffs=infile_data['DIF']
    
    # write output file comparing L1top and SFR elevations
    ofp=open(outfile,'w')
    ofp.write('MSEG,STOP,TOPNEW,DIF\n')
    for i in range(len(infile_data['MSEG'])):
        ofp.write("%s,%s,%s,%s\n" %(infile_data['MSEG'][i],infile_data['STOP'][i],infile_data['TOPNEW'][i],infile_data['DIF'][i]))
    ofp.close()    
    
    return(diffs)
    
    
    # comparative summary statistics
def plot_diffstats(diffs,title):
    outfile=PdfPages(title+'.pdf')
    std,mx,mn=np.std(diffs),np.max(diffs),np.min(diffs)
    floating=[d for d in diffs if d<0]
    incised=[d for d in diffs if d>0]
    nfloating, nincised=len(floating),len(incised)
    meanfloating,meanincised=np.mean(floating),np.mean(incised)
    
    # plot histogram
    hist,bins=np.histogram(diffs,bins=100)
    cuttoff=round(len(diffs)*0.001) # only plot bins with 0.1% or more of segments
    inds=np.where(hist>cuttoff)
    widths=1*(bins[1]-bins[0])
    plt.figure()
    plt.bar(bins[inds],hist[inds]/float(len(diffs)),width=widths)
    plt.title(title)
    plt.xlabel('Stream bottom elevation, feet below land surface (negative=floating)')
    plt.ylabel('Portion of stream reaches')
    plt.title(title)
    outfile.savefig()
    
    # CDF
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.hist(diffs,bins=1000,normed=True,cumulative=True,histtype='step')
    xlim=[bins[inds][0],bins[inds][-1]]
    ax.set_xlim(xlim)
    ax.set_title(title)
    ax.set_xlabel('Stream bottom elevation, feet below land surface (negative=floating)')
    ax.set_ylabel('Portion of stream reaches')
    plt.grid(True)
    summary='SummaryStatistics:\nMax floating: %s\nMax incised: %s\nStdev: %s\nNum floating: %s\nNum incised: %s\nMean floating: %s\nMean incised: %s' %(mx,mn,std,nfloating,nincised,meanfloating,meanincised)
    ax.text(0.95,.05,summary,verticalalignment='bottom',horizontalalignment='right',transform=ax.transAxes,bbox={'facecolor':'white','pad':10})    
    outfile.savefig()
    outfile.close()



'''
L1tops=dict()
for r in range(nrows):
    for c in range(ncols):
        cellnum=r*ncols+c+1
        L1tops[cellnum]=topdata[r,c]
    
    DEM=open(L1top,'r')
    gwv=dict()
    for i in range(0,nrow):
        line=DEM.readline()
        vals=re.split('\s+',line.strip())
        for j in range(0,ncol):
            cellnum=(i*ncol)+(j+1)
            gwv[cellnum]=float(vals[j])
    
    GWVmat1_data=np.genfromtxt(MAT1, delimiter=',',names=True,dtype=None)
    
    infile_data=defaultdict(list)
    for i in range(len(GWVmat1_data)):
        row,column=GWVmat1_data['row'][i],GWVmat1_data['column'][i]
        cellnum=(row-1)*ncol+column
        infile_data['MSEG'].append(GWVmat1_data['segment'][i])
        infile_data['STOP'].append(GWVmat1_data['top_streambed'][i])
        infile_data['TOPNEW'].append(gwv[cellnum])
        infile_data['DIF'].append(gwv[cellnum]-GWVmat1_data['top_streambed'][i])
    ofp=open('STOP_comparison.csv','w')
    ofp.write('MSEG,STOP,TOPNEW,DIF\n')
    for i in range(len(infile_data['MSEG'])):
        ofp.write("%s,%s,%s,%s\n" %(infile_data['MSEG'][i],infile_data['STOP'][i],infile_data['TOPNEW'][i],infile_data['DIF'][i]))
    ofp.close()
'''

    