import numpy as np
from collections import defaultdict
import pdb

def get_upSEGs(GWVmat2):
    upSEGs=defaultdict(list)
    for i in range(len(GWVmat2)):
        segment=GWVmat2['outseg'][i]
        upSEGs[segment].append(GWVmat2['segment'][i])
    return upSEGs
    
def get_upSTOP(segment,upSEGs,Seg_ends,GWVmat1):
    upSEGz=upSEGs[segment]
    try: # if upstream segment has already been updated, get segment ends
        upSTOPs=[]
        for segs in upSEGz:
            upSTOPs.append(Seg_ends[segs][-1])
        upSTOP=np.min(upSTOPs)
    except: # otherwise get them from the original Mat1 file
        print "Fixed 2nd order segment before 1st order"
        upSTOPs=[]
        for segs in upSEGz:
            seg_indsGWVmat=list(np.where(GWVmat1['segment']==segs)[0])
            upSTOPs.append(GWVmat1['top_streambed'][seg_indsGWVmat][-1])
        upSTOP=np.min(upSTOPs)
    return upSTOP,upSEGz

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

def check4backwards_routing(Seg_ends,GWVmat2,tolerance):
    GWVmat1=[]
    nbackwards=0
    rec_fname='backwards_routed_segments.txt'
    ofp=open(rec_fname,'w')
    ofp.write('Segments left with upstream segments that are higher than downstream segments:')
    ofp.write('segment,upSEGs,downSEG,upSTOP,reach_elevations,dnSTOP\n')
    upSEGs=get_upSEGs(GWVmat2)
    for segment in Seg_ends.keys():
        STOPs=Seg_ends[segment]
        dnSEG_ind=np.where(GWVmat2['segment']==segment)[0]
        dnSEG=GWVmat2['outseg'][dnSEG_ind][0]
        
        # if segment is terminal, continue
        if dnSEG==0:
            continue
        dnSTOP=Seg_ends[dnSEG][-1]
        
        # if segment isn't a headwater, compare upstream/downstream segment elevations
        if segment in upSEGs:        
            upSTOP,upSEGz=get_upSTOP(segment,upSEGs,Seg_ends,GWVmat1)
        
            if dnSTOP>upSTOP+tolerance:
                nbackwards+=1
                ofp.write('%s,{%s},%s,%s,{%s},%s\n' %(segment,','.join(map(str,upSEGz)),dnSEG,upSTOP,','.join(map(str,STOPs)),dnSTOP))
    ofp.close()
    return nbackwards,rec_fname
    
def fix_last_reaches(segment,Seg_ends,GWVmat2):
    
    dnSEG_ind=np.where(GWVmat2['segment']==segment)[0]
    dnSEG=GWVmat2['outseg'][dnSEG_ind][0]
    lastSTOP=Seg_ends[segment][-1]

    if dnSEG==0:
        return lastSTOP
    dnSEG_reach1STOP=Seg_ends[dnSEG][0]
    if (lastSTOP-dnSEG_reach1STOP)>1:
        lastSTOP=dnSEG_reach1STOP+1
    return lastSTOP
        
def find_nearby_SFR(GWVmat1_data,segnum,firstlast,maxsearch):
    
    # determine cellnum based on GWVmat1, segment and first/last    
    GWVmat1_inds=np.where(GWVmat1_data['segment']==segnum)[0]
    
    if firstlast=='first':
        reach_inds=np.where(GWVmat1_data[GWVmat1_inds]['reach']==1)[0]
    elif firstlast=='last':
        lastreach=np.max(GWVmat1_data['reach'][GWVmat1_inds])
        reach_inds=np.where(GWVmat1_data['reach'][GWVmat1_inds]==lastreach)
    
    row=GWVmat1_data[GWVmat1_inds][reach_inds]['row']
    column=GWVmat1_data[GWVmat1_inds][reach_inds]['column']
    
    # start in upper left corner of box defined by maxsearch
    ul_row=row-maxsearch
    ul_col=column-maxsearch
    neighbors=[] # list of row,col tuples for surrounding maxsearch cells
    for c in range(maxsearch*2+1):
        for r in range(maxsearch*2+1):
            crow=ul_row+r
            ccol=ul_col+c
            if crow==row and ccol==column:
                continue # exclude current cell
            else:
                neighbors.append((crow,ccol))

    # search through neighbors and find nearby SFR cells
    nearbySFR=[]
    for n in neighbors:
        try:
            SFR_index=np.argwhere((GWVmat1_data['row']==n[0]) & (GWVmat1_data['column']==n[1]))[0][0]
            tempnum=GWVmat1_data['segment'][SFR_index]
            nearbySFR.append(tempnum)
        except IndexError:
            continue
    
    # toss neighbors that are in the same segment        
    nearbySFR=np.unique([s for s in nearbySFR if s<>segnum])
    #print 'nearby SFR cells: %s' %(nearbySFR)
    return nearbySFR


def find_outSEG(segnum,GWVmat1,nearbySFR,outseg_old):
    tolerance=1.0 # new outSEG must be this much lower than old for the routing to be switched
    
    # get streambed top elev in current segment
    inds=np.where(GWVmat1['segment']==segnum)
    current_elev=np.min(GWVmat1['top_streambed'][inds])
    
    # get end elevations of nearby segments (could instead use reaches)
    dnElev=9999
    for segnum in nearbySFR:
        inds=np.where(GWVmat1['segment']==segnum)
        elev=np.max(GWVmat1['top_streambed'][inds])
        #print 'elev= %s' %(elev)
        # make sure that dnElev is < current elev
        if elev<dnElev and elev<=current_elev:
            dnElev=elev
            outSEG=segnum

    #if dnElev==9999:
        #outSEG=
    
    # only pick new outSEG if appreciably lower than existing
    inds=np.where(GWVmat1['segment']==outseg_old)
    old_dnELEV=np.min(GWVmat1['top_streambed'][inds])
    if dnElev>old_dnELEV-tolerance:
        outSEG=outseg_old
    return outSEG
            
    
def fixends(segment,STOP_elevs,TOP_elevs,upSEGs,Seg_ends,GWVmat1,GWVmat2,headwater_incise,outfile_object,warningsfile_object,nbackwards,tolerance):
    headwater=False
    if segment not in GWVmat2['outseg']:
        headwater=True
        
    def fixend(segment,TOP_elev,STOP_elev,Seg_ends,upSTOP,dnSTOP,incise,firstlast,headwater,terminal,warningsfile_object,nbackwards):
        
        # if original STOP is between land surface and an incised depth parameter, 
        # and if it is between up/down segment elevations, leave it alone
        # if incise is zero, this option will probably be passed
        if TOP_elev>=STOP_elev>=(TOP_elev-incise) and upSTOP>STOP_elev>dnSTOP:
            newelev=STOP_elev
            
        # otherwise, if an ideal elevation for SFR STOP is between min/max elevs of up/down segements, set it
        # when readjusting; make sure that
        elif upSTOP>=TOP_elev-incise>=dnSTOP:
            newelev=TOP_elev-incise
            
        # else, this means that the properly route elevation is not ideal (not within the interval between ls and incise depth)
        # if the segment is terminal, set new elevation to a specified increment below minimum elev in UP segment
        # if there is no gradient, set STOP= upSTOP
        # otherwise, figure out if previous STOP (from SFR_utilities) is above or below, and then start at either end of the interval # and adjust until new STOP is between upSTOP and dnSTOP
        
        # if terminal segment and above were false, just set to 1 ft below upSTOP
        elif terminal:
            newelev=upSTOP-1
            
        elif upSTOP==dnSTOP:
            newelev=upSTOP

        elif headwater and upSTOP<dnSTOP:
            newelev=upSTOP # upSTOP was already set to model TOP elevation - incise parameter
        
        # if adjusting the last reach and segment below is higher than segment above:
        elif firstlast=='last' and upSTOP<dnSTOP:
            if headwater:
                newelev=upSTOP # the upstream STOP should already have been appropriately adjusted to TOP - incise.
                if dnSEG in Seg_ends: # this means that dnSEG will not be subsequently adjusted; add a warning (backwards segment)
                    warningsfile_object.write('%s,%s,%s,headwater; reach 1 set to Land surface - %s ft.; downstream segment is higher\n' %(segment,upSTOP,dnSTOP,incise))
                    nbackwards+=1
            # if dnSEG not yet been adjusted, leave elevation alone
            if dnSEG not in Seg_ends: 
                newelev=STOP_elev
                nbackwards+=1
            # otherwise, if downstream segment has been adjusted, set elevation equal reach 1 and record
            else:
                newelev=upSTOP
                warningsfile_object.write('%s,%s,%s,not a headwater; downstream segment has already been adjusted and is higher than upstream segment\n' %(segment,upSTOP,dnSTOP))
                nbackwards+=1
                
        else: # upSTOP>dnSTOP, meaning there is room to adjust the segment end elevations
            print "fiddling with %s reach..." %(firstlast)
            print "upSTOP is %s, downSTOP is %s" %(upSTOP,dnSTOP)
            newelev=STOP_elev # original streambed top elevation from SFR_utilities
            
            # if elevation is floating above land surface:
            # since it failed the tests above, it 'has' to float; however try to minimiz floating (within constraints of upSTOP/dnSTOP)
            # if the adjustment places the elevation above the upSTOP, last algorithm should correct
            if newelev>TOP_elev:
                print "land surface is %s; adjusting %s reach elevation to minimize floating..." %(TOP_elev,firstlast) 
                if dnSTOP>TOP_elev:                
                    newelev=dnSTOP+0.1
                elif dnSTOP>TOP_elev-incise: 
                    newelev=TOP_elev-0.1
                else:
                    newelev=TOP_elev-incise
                    
            Above=False
            Below=False

            # if elevation is above minimum in UP segment
            if newelev>upSTOP:
                Above=True
                newelev=TOP_elev-incise
            # or if elevation is below maximum in DOWN segment
            if newelev<dnSTOP:
                Above=False
                Below=True
                newelev=TOP_elev
            
            # algorithm to maintain downstream routing with adjustments of some higher order segments prior to lower order
            # (such adjustments may disrupt the downstream monotonicity from before)
            n=0.1
            overshoot=False
            while Above or Below:
                if n>100: 
                    #backwards=True
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


        return newelev,nbackwards
  
    terminal=False
    
    # Get max elevation of downstream segment
    seg_indsGWVmat2=np.where(GWVmat2['segment']==segment)[0][0]
    dnSEG=GWVmat2['outseg'][seg_indsGWVmat2]
    if dnSEG==0: # this is a terminal segment
        dnSTOP=TOP_elevs[-1]
        terminal=True
    if dnSEG>0: # if downstream segment has already been updated, get ends for downstream segment
        try:
            dnSTOP=Seg_ends[dnSEG][0]
        except: # otherwise get them from the original Mat1 file
            seg_indsGWVmat=list(np.where(GWVmat1['segment']==dnSEG)[0])
            dnSTOP=GWVmat1['top_streambed'][seg_indsGWVmat][0]
        
    # Get min elevation of upstream segment (land surface if headwater)
    if headwater:
        upSEGz='headwater'
        upSTOP=TOP_elevs[0]-headwater_incise
        if TOP_elevs[0]<dnSTOP and dnSEG in Seg_ends.keys():
            warningsfile_object.write('%s,%s,%s,L1top elevation is below downstream segment, which was already adjusted\n' %(segment,upSTOP,dnSTOP))
            first=TOP_elevs[0]
        else:
            #first,nbackwards=fixend(segment,TOP_elevs[0],STOP_elevs[0],Seg_ends,upSTOP,dnSTOP,headwater_incise,'first',headwater,terminal,warningsfile_object,nbackwards)        
            first=upSTOP
            
    if not headwater:
        upSTOP,upSEGz=get_upSTOP(segment,upSEGs,Seg_ends,GWVmat1)
        
        # adjust reach 1 elevation based on upstream/downstream segment elevations, land surface, and incise parameter
        first,nbackwards=fixend(segment,TOP_elevs[0],STOP_elevs[0],Seg_ends,upSTOP,dnSTOP,0,'first',headwater,terminal,warningsfile_object,nbackwards)

    # if downstream segment is > upstream segment and downstream was already updated, quit and write to error file    
    if dnSTOP>upSTOP+tolerance and dnSEG in Seg_ends.keys():
        warningsfile_object.write('%s,%s,%s,downstream segment is higher than upstream segment; downstream segment already adjusted\n' %(segment,upSTOP,dnSTOP))
        print "Warning, backwards-routed segment! see Fix_ends warning file for details"
        nbackwards+=1
        return dnSTOP,dnSTOP,nbackwards
    
    
    # if there is more than one reach in the segment, adjust the last reach elevation based on same criteria
    if len(STOP_elevs)>1:
        if dnSTOP>first: # this shouldn't happen
            first=upSTOP # set first reach equal to min elev in upSEG, should be above downstop (otherwise function would have ended)
        last,nbackwards=fixend(segment,TOP_elevs[-1],STOP_elevs[-1],Seg_ends,first,dnSTOP,0,'last',headwater,terminal,warningsfile_object,nbackwards)
    else:
        last=first
        
    # in case the first and last reaches are between the adjacent STOPs, but backwards
    if last>first:
        print 'last reach higher than first! last=%s, first=%s' %(last,first)
        print 'upSTOP=%s,dnSTOP=%s' %(upSTOP,dnSTOP)
        print 'upSEGs=%s' %(upSEGz)
        print 'Top elevations: %s' %(TOP_elevs)
        print 'terminal=%s' %(terminal)
        backwards=True
        n=1.0
        
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
   
    outfile_object.write('%s,%s,%s,%s,%s,%s,%s,%s,%s\n' %(segment,STOP_elevs[0],STOP_elevs[-1],'_'.join(map(str,upSEGz)),dnSEG,upSTOP,dnSTOP,first,last))
    
    print 'Reach 1 adjusted elev is %s, Last reach adjusted elev is %s' %(first,last)
    return first,last,nbackwards

