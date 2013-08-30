import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def get_sinkers(infile,segment_col):
    indata=np.genfromtxt(infile,delimiter=',',dtype=None,names=True)
    segments=np.unique(indata[segment_col])
    return segments  
    
    
def plot_profiles(segs2plot,seg_distdict,TOPNEWdict,STOP2dict,Bottomsdict,pdffile):
#segs2plot=[2,1081,324,1188,1060,1201,1076,341,1113,1131,949,523,408,469,1064,200,300,400,500,600,700,800,900,1000,1100,1200]
#segs2plot=segments[::int(np.floor(len(segments)/num_segs2plot))]
    segs2plot=get_sinkers('below_bot.csv','segment')
    segs2plot=sorted(segs2plot)
    
    pdf=PdfPages(pdffile)
    
    print "saving plots of selected segments to " + pdffile
    for seg in segs2plot:
        
        fig=plt.figure()
        ax1=plt.subplot(111)
        ax1.grid(True)
        p1=ax1.plot(seg_distdict[seg],TOPNEWdict[seg],'b',label='land surface')
        #p2=ax1.plot(seg_distdict[seg],STOP1dict[seg],'r',label='sfr_utilities')
        if len(seg_distdict[seg])==len(STOP2dict[seg]):
            p3=ax1.plot(seg_distdict[seg],STOP2dict[seg],'g',label='fix_w_DEM')
        p4=ax1.plot(seg_distdict[seg],Bottomsdict[seg],'k',label='model bottom')
        handles, labels = ax1.get_legend_handles_labels()
        ax1.legend(handles,labels)
        ax1.set_title('segment ' + str(seg))
        plt.xlabel('distance along segment (ft.)')
        ax1.set_ylabel('Elevation (ft)')
    
        pdf.savefig(fig)
    pdf.close()
    plt.close('all')
