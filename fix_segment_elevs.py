# Fix segment elevs via linear interpolation
# So that elevs are monotonically downhill from ELEVMAX of first reach to ELEVMIN of last reach
#
# Inputs:
# feature (feature layer containing stream segments)
# COMID
# first_cellnum (cellnum containing reach 1)
# last_cellnum (cellnum containing last reach)
# max_elev (max elevation of first reach)
# min_elev (min elevation of last reach)

import arcpy



def fixsegelevs(feature,COMID,first_cellnum,last_cellnum,max_elev,min_elev):
    if arcpy.Exists("segment"):
        arcpy.Delete_management("segment")
    if arcpy.Exists("first_reach"):
        arcpy.Delete_management("first_reach")
    arcpy.SelectLayerByAttribute_management(feature,"CLEAR_SELECTION")
    # create temp layer of reaches in segment
    theQuery="COMID= "+ str(COMID)
    arcpy.SelectLayerByAttribute_management(feature,"NEW_SELECTION",theQuery)
    arcpy.MakeFeatureLayer_management(feature,"segment")
    num_reaches=int(arcpy.GetCount_management("segment").getOutput(0))
    selected=arcpy.SearchCursor("segment")

    dh=min_elev-max_elev
    # get total length of segment
    dist=0
    for rows in selected:
        clength=rows.getValue("LengthFt")
        dist+=clength
    slope=dh/dist
    print "dh=%s, segment length=%s, slope=%s" %(dh,dist,slope)
    del selected
    
    print "assigning elevations for %s reaches" %(num_reaches)
    # create temp layer containing first reach; set max elevation in that reach to max_elev
    theQuery="node= "+ str(first_cellnum)
    arcpy.SelectLayerByAttribute_management("segment","NEW_SELECTION",theQuery)
    arcpy.MakeFeatureLayer_management("segment","first_reach")
    selected=arcpy.UpdateCursor("first_reach")
    for rows in selected:
        print "reach 1",
        length=rows.LengthFt
        rows.ELEVMAX=max_elev
        rows.ELEVAVE=max_elev+slope*0.5*length
        rows.ELEVMIN=max_elev+slope*length
        prev_min=max_elev+slope*length
        print "%s,%s,%s" %(rows.ELEVMAX,rows.ELEVAVE,rows.ELEVMIN)
        selected.updateRow(rows)
        x_end,y_end=rows.X_end,rows.Y_end
    del selected

    # then for each reach in segment, find next that has start coords equal to previous end coords
    for i in range(num_reaches-1):
        arcpy.SelectLayerByAttribute_management("segment","CLEAR_SELECTION")
        selected=arcpy.UpdateCursor("segment")
        for rows in selected:
            x_start,y_start=rows.X_start,rows.Y_start    
            if x_start==x_end and y_start==y_end:
                print "reach %s" %(i+2),
                length=rows.LengthFt
                rows.ELEVMAX=prev_min
                rows.ELEVAVE=prev_min+slope*0.5*length
                rows.ELEVMIN=prev_min+slope*length
                prev_min=prev_min+slope*length
                print "%s,%s,%s" %(rows.ELEVMAX,rows.ELEVAVE,rows.ELEVMIN)
                selected.updateRow(rows)
                x_end,y_end=rows.X_end,rows.Y_end
                del selected
                break
            else:
                continue
    