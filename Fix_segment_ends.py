import numpy as np
from collections import defaultdict
import pdb

def get_upSEGs(GWVmat2):
    upSEGs=defaultdict(list)
    for i in range(len(GWVmat2)):
        segment=GWVmat2['outseg'][i]
        upSEGs[segment].append(GWVmat2['segment'][i])
    return upSEGs

def get_dnSTOP(segment,GWVmat1,dnSEG,searchdist):
    # NOTE: according to SFR manual, p2, trib inflows are only allowed in first reach
    # so this function is unecessary
    
    #look in 8 nearest cells for dnSEG
    
    seg_indsGWVmat1=list(np.where(GWVmat1['segment']==segment)[0])
    
    # start at row column of last (min. el) reach in segment
    row,col=GWVmat1['row'][seg_indsGWVmat1[-1]],GWVmat1['column'][seg_indsGWVmat1[-1]]
    rows,columns=[row-searchdist,row,row+searchdist],[col-searchdist,col,col+searchdist]
    
    # find indicies of cells in GWVmat1 (with some ugly reformatting)
    rowinds=tuple([a[0] for a in [np.where(GWVmat1['row']==i) for i in rows]])
    rowinds=np.concatenate(rowinds,axis=0)
    colinds=tuple([a[0] for a in [np.where(GWVmat1['column']==i) for i in columns]])
    colinds=np.concatenate(colinds,axis=0)
    intersect_rc=np.intersect1d(rowinds,colinds)
    GWVmat1_ind=np.intersect1d(np.where(GWVmat1['segment']==dnSEG)[0],intersect_rc)
    
    # get STOP of downstream reach (if multiple, use minimum)

    dnSTOP=np.min(GWVmat1[GWVmat1_ind]['top_streambed'])
    return dnSTOP

def check4backwards_routing(Seg_ends,GWVmat2):
    nbackwards=0
    rec_fname='backwards_routed_segments.txt'
    ofp=open(rec_fname,'w')
    ofp.write('segment,upSEGs,downSEG,reach_elevations\n')
    for segment in Seg_ends.keys():
        STOPs=Seg_ends[segment]
        if STOPs[-1]>STOPs[0]:
            nbackwards+=1
            upSEGs=get_upSEGs(GWVmat2)
            upSEGz=upSEGs[segment]
            dnSEG=GWVmat2['outseg'][segment-1]
            ofp.write('%s,%s,%s,%s\n' %(segment,' '.join(map(str,upSEGz)),dnSEG,' '.join(map(str,STOPs))))
    ofp.close()
    return nbackwards,rec_fname
        
        
    
def fixends(segment,STOP_elevs,TOP_elevs,upSEGs,Seg_ends,GWVmat1,GWVmat2,headwater_incise,outfile_object,nbackwards):
    headwater=False
    if segment not in GWVmat2['outseg']:
        headwater=True
               
    def fixend(segment,TOP_elev,STOP_elev,Seg_ends,upSTOP,dnSTOP,incise,firstlast,headwater,terminal,nbackwards):
        backwards=False
        
        # if original STOP is between land surface and an incised depth parameter, leave it alone
        # if incise is zero, this option will probably be passed
        if TOP_elev>=STOP_elev>=(TOP_elev-incise):
            newelev=STOP_elev
            
        # otherwise, if an ideal elevation for SFR STOP is between min/max elevs of up/down segements, set it
        elif upSTOP>=TOP_elev-incise>=dnSTOP:
            newelev=TOP_elev-incise
            
        # else, this means that the properly route elevation is not ideal (not within the interval between ls and incise depth)
        # if there is no gradient, set STOP= upSTOP
        # otherwise, figure out if previous STOP (from SFR_utilities) is above or below, and then start at either end of the interval # and adjust until new STOP is between upSTOP and dnSTOP
        
        # if terminal segment and above were false, just set to 1 ft below upSTOP
        elif terminal:
            newelev=upSTOP-1
            
        elif upSTOP==dnSTOP:
            newelev=upSTOP

        else:
            print "fiddling with %s reach..." %(firstlast)
            print "upSTOP is %s, downSTOP is %s" %(upSTOP,dnSTOP)
            newelev=STOP_elev
            
            Above=False
            Below=False
            
            if newelev>upSTOP:
                Above=True
                newelev=TOP_elev-incise
            elif newelev<dnSTOP:
                Below=True
                newelev=TOP_elev
                
            n=0.1
            overshoot=False
            while Above or Below:
                if n>100: # if upSTOP<dnSTOP, might be able to add an adjustment if headwater
                    backwards=True
                    break
                increment=1.0/n
                print newelev
                # this happens if the new elev overshoots (see below)
                if overshoot:
                    n*=10
                    overshoot=False
                    continue
                if Above:
                    newelev=newelev-increment
                    if newelev<upSTOP:
                        Above=False
                    if newelev<dnSTOP:
                        Above=False
                        Below=True
                        overshoot=True
                        continue

                elif Below:
                    newelev=newelev+increment
                    if newelev>dnSTOP:
                        Below=False
                    if newelev>upSTOP:
                        Below=False
                        Above=True
                        overshoot=True

        if firstlast=='last' and upSTOP<dnSTOP:
            backwards=True
            nbackwards+=1
        if backwards:
            if headwater: # Let it alone; probably floating- start can't be adjusted downward because next segment is too high
                newelev=STOP_elev
            if dnSEG not in Seg_ends:
                newelev=STOP_elev
            # if still backwards, could try searching for othersegments within 2 or more cells
            else:
                pdb.set_trace()
                newelev=9999
        return newelev,nbackwards
  
    terminal=False
    # Get max elevation of downstream segment
    dnSEG=GWVmat2['outseg'][segment-1]
    if dnSEG==0: # this is a terminal segment
        dnSTOP=TOP_elevs[-1]
        terminal=True
    else:
        #dnSTOP=get_dnSTOP(segment,GWVmat1,dnSEG) # this Fn would apply if routing could go to any reach
        seg_indsGWVmat=list(np.where(GWVmat1['segment']==dnSEG)[0])
        dnSTOP=GWVmat1['top_streambed'][seg_indsGWVmat][0]
    
    # Get min elevation of upstream segment (land surface if headwater)
    if headwater:
        first,nbackwards=fixend(segment,TOP_elevs[0],STOP_elevs[0],Seg_ends,TOP_elevs[0],dnSTOP,headwater_incise,'first',headwater,terminal,nbackwards)        
        upSEGz='headwater'
        upSTOP=TOP_elevs[0]
    else:
        upSEGz=upSEGs[segment]
        try:
            upSTOPs=[]
            for segs in upSEGz:
                upSTOPs.append(Seg_ends[segs][-1])
            upSTOP=np.min(upSTOPs)
        except:
            print "Fixed 2nd order segment before 1st order"
            upSTOPs=[]
            for segs in upSEGz:
                seg_indsGWVmat=list(np.where(GWVmat1['segment']==segs)[0])
                upSTOPs.append(GWVmat1['top_streambed'][seg_indsGWVmat][-1])
            upSTOP=np.min(upSTOPs)
        first,nbackwards=fixend(segment,TOP_elevs[0],STOP_elevs[0],Seg_ends,upSTOP,dnSTOP,0,'first',headwater,terminal,nbackwards)

    if len(STOP_elevs)>1:
        last,nbackwards=fixend(segment,TOP_elevs[-1],STOP_elevs[-1],Seg_ends,upSTOP,dnSTOP,0,'last',headwater,terminal,nbackwards)
    else:
        last=first
        
    # in case the first and last reaches are between the adjacent STOPs, but backwards
    if last>first:
        backwards=True
        n=10.0
        
        overshoot=False
        while backwards:
            increment=1.0/n
            
            if overshoot:
                first-=increment
                if first<upSTOP:
                    overshoot=False
                    
            if first<last:
                first+=increment
                if first>upSTOP:
                    overshoot=True
                    n*=10
                    continue
            else:
                backwards=False
   
    outfile_object.write('%s,%s,%s,%s,%s,%s,%s,%s,%s\n' %(segment,STOP_elevs[0],STOP_elevs[-1],'_'.join(map(str,upSEGz)),dnSEG,upSTOP,dnSTOP,first,first))
    
    print 'Reach 1 adjusted elev is %s, Last reach adjusted elev is %s' %(first,last)
    return first,last,nbackwards

