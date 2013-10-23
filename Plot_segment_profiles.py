import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pdb

def get_sinkers(infile,segment_col):
    indata=np.genfromtxt(infile,delimiter=',',dtype=None,names=True)
    segments=np.unique(indata[segment_col])
    return segments  
    
    
def plot_profiles(segs2plot,seg_distdict,profiles,profilenames,pdffile,**kwargs):
    
    # segs2plot= list of segments to plot
    # seg_distdict= list of distances along segments
    # profiles= list of dictionaries containing profiles for each segment
    # profilenames= list of names, one for each type of profile e.g., model top, STOP post-fix_w_DEM, etc.
    # pdffile= name for output pdf
    
    try:
        Bottomsdict=kwargs['Bottoms']
        Bottoms=True
    except KeyError:
        Bottoms=False
    try:
        plot_slopes=kwargs['plot_slopes']
        Slopesdict=kwargs['slopes']
        Reach_lengthsdict=kwargs['reach_lengths']
    except KeyError:
        plot_slopes=False
    
    # function to reshape distances and elevations to plot actual cell elevations
    def reshape_seglist(seg_dict,distance):
        if distance:
            seg_list=[0]
        else:
            seg_list=[]
        for i in range(len(seg_dict)):
            seg_list.append(seg_dict[i])
            seg_list.append(seg_dict[i])
        if distance:
            seg_list=seg_list[:-1] # trim last, since last seg distances denotes end of last reach
        return seg_list
        
    pdf=PdfPages(pdffile)
    print "saving plots of selected segments to " + pdffile
    for seg in segs2plot:
        print seg
        # reshape distances and elevations to plot actual cell elevations
        seg_distances=reshape_seglist(seg_distdict[seg],True)
        profiles2plot=[]
        for i in range(len(profiles)):
            profile=reshape_seglist(profiles[i][seg],False)
            profiles2plot.append(profile)
            
        if Bottoms:
            seg_Bots=reshape_seglist(Bottomsdict[seg],False)
            profiles2plot.append(seg_Bots)
            profilenames.append('model bottom')
        
        if plot_slopes:
            slopes=reshape_seglist(Slopesdict[seg],False)
            reachlengths=reshape_seglist(Reach_lengthsdict[seg],False)


        fig=plt.figure()
        if plot_slopes:
            ((ax1, ax2)) = fig.add_subplot(2,1,sharex=True,sharey=False)
        else:
            ax1=fig.add_subplot(1,1,1)
        ax1.grid(True)
        colors=['b','g','r','k']

        for i in range(len(profiles2plot)):
            ax1.plot(seg_distances,profiles2plot[i],color=colors[i],label=profilenames[i])

        handles, labels = ax1.get_legend_handles_labels()
        ax1.legend(handles,labels,loc='best')
        ax1.set_title('segment ' + str(seg))
        plt.xlabel('distance along segment (ft.)')
        ax1.set_ylabel('Elevation (ft)')
        
        # adjust limits to make all profiles visible
        ymax,ymin=np.max(profiles2plot),np.min(profiles2plot)
        ax1.set_ylim([ymin-10,ymax+10])
        
        # plot segment slopes if desired
        if plot_slopes:
            ax2.grid(True)
            ax2.plot(seg_distances,slopes,color='0.75',label='streambed slopes')
            ax2.set_ylabel('Streambed slope')
            ax3=ax2.twinx()
            ax3.plot(seg_distdict[seg],reachlengths,'b',label='reach length')
            ax3.set_ylabel('reach length (ft)')
            handles, labels = ax2.get_legend_handles_labels()
            ax2.legend(handles,labels)
            ax3.legend(loc=0) 
        pdf.savefig(fig)
    pdf.close()
    plt.close('all')
