__author__ = 'aleaf'

import sys
#sys.path.append('../GIS_utils')
sys.path.append('D:/ATLData/Documents/GitHub/GIS_utils')
import GISio

import pandas as pd
import numpy as np
import os
import fiona
from shapely.geometry import Point, LineString, shape, mapping
import re
import discomb_utilities as disutil
import shutil
#import smooth_streambed as sm
from collections import defaultdict

# input files
path = 'D:/ATLData' # path to prepend input files below (for switching between PC an Mac)
MFgrid = path + '/LittlePlover/grid/LPRgrid_atl.shp'
MFdomain = path + '/LittlePlover/input/LPR_model_nearfield_revised_line.shp' # must be a single part line
MFgrid_node_attribute = 'node' # name of attribute field containing node numbers
DEM = path + '/LittlePlover/input/DEM10mwtm_ft'

# if specifying multiple features to merge, they must not overlap!
stream_linework = [path + '/LittlePlover/input/LPR_ditches.shp',
                   path + '/LittlePlover/input/canal_flowlines.shp',
                   path + '/LittlePlover/input/flowlines_clipped_mk.shp',
                   path + '/LittlePlover/input/hay_meadow_crk_etc.shp']

# settings
preproc = True # whether or not to run preprocessing routine
from_scratch = True # True/False: whether a new SFR dataset is being created (instead of adding to existing data)
DEM_z_conversion_factor = 1
reach_cutoff = 1.0 # model units
max_distance_from_SFR = 2 * np.sqrt(2 * 100 ** 2) # max distance that new linework can be away from existing SFR cells
max_distance_to_model_boundary = 650 # end points within this distance of model domain edge will be considered outlets
routing_tol = 100 # search distance when looking for adjacent upstream / downstream segment ends / starts for routing segments
width_in_cell = 5 # constant width for all new SFR reaches
bed_slope = 1e-4

# pull generic settings from Mat1
bed_K = 1
bed_thickness = 3
bed_roughness = 0.037


existing_SFR_linework = None # (must be broken up by cell and contain cellnum)
existing_sfr_shp = None
existing_sfr_node_attribute = None # name of attribute field containing node numbers (also applies to existing linework)
existing_sfr_elevation_attribute = None
Mat1file = None
Mat2file = None
DISfile = path + '/LittlePlover/input/LPR_Mod_refined30m.dis'
nreaches = 0

if not from_scratch:
    # read in SFR tables
    Mat1 = pd.read_csv(Mat1file)
    Mat2 = pd.read_csv(Mat2file)
    nsegments = int(np.max(Mat1['segment']))
    nreaches = len(Mat1['segment'])



# intermediate files
stream_cells = path + '/LittlePlover/intermediate/new_stream_cells.shp'
stream_cells_dissolve = path + '/LittlePlover/intermediate/new_stream_cells_dissolve.shp'
stream_fragments = path + '/LittlePlover/intermediate/new_streams_fragments.shp'
stream_fragments_points = path + '/LittlePlover/intermediate/new_streams_fragments_points.shp'

# output files
Mat1_updated = path + '/LittlePlover/SFRoutput/Mat1.csv'
Mat2_updated = path + '/LittlePlover/SFRoutput/Mat2.csv'
stream_cells_updated = path + '/LittlePlover/SFRoutput/SFR_cells.shp'
SFR_linework_updated = path + '/LittlePlover/SFRoutput/SFR_lines.shp'


def shp2df(shp):
    '''
    Read shapefile into Pandas dataframe
    shp = shapefile name
    IDfield =
    '''
    print "loading attributes from {} into pandas dataframe...".format(shp)
    shp_obj = fiona.open(shp, 'r')

    attributes_dict = {}
    knt = 0
    length = len(shp_obj)
    for line in shp_obj:

        props = line['properties']
        attributes_dict[line['id']] = props
        knt += 1
        print '\r{:d}%'.format(100*knt/length),
    print '\n'
    # convert to pandas dataframe, join in centroids, sort by FID
    shp_df = pd.DataFrame.from_dict(attributes_dict, orient='index')

    return shp_df

# if adding to existing SFR dataset
if not from_scratch:

    print "creating pandas DataFrame of existing SFR cell information..."
    # use dictionaries for speed, convert later
    existing_SFR = fiona.open(existing_sfr_shp, 'r')
    existing_sfr_cells = {}
    centroids = {}
    knt = 0
    ncells = len(existing_SFR)
    for polygon in existing_SFR:
        # get the geometry, calculate the centroid and scrape out the coordinates with an re.
        # ugly but looks like the way to do it
        geometry = shape(polygon['geometry'])
        regexp = '(?<=\().+?(?=\))' # return everything between the ()
        centroid = map(float, re.findall(regexp, geometry.centroid.wkt)[0].split())
        FID = polygon['id']
        attributes = polygon['properties']
        existing_sfr_cells[FID] = attributes
        centroids[FID] = centroid
        knt +=1
        print '\r{:d}%'.format(100*knt/ncells),

    # convert dictionaries to pandas dataframe, join in centroids, sort by FID
    existing_SFR_df = pd.DataFrame.from_dict(existing_sfr_cells, orient='index')
    df = pd.DataFrame.from_dict(centroids, orient='index')
    df.columns = ['X', 'Y']
    existing_SFR_df = existing_SFR_df.astype(int)
    existing_SFR_df = existing_SFR_df.join(df).sort()
    existing_SFR_df.index = existing_SFR_df.index.astype(int)
else:
    existing_SFR_df = None

ofp = open('addSFR.log', 'w')
ofp.write('Log file for adding streams to SFR network\n\n')

def preprocessing(stream_linework):

    import arcpy

    # initialize the arcpy environment
    arcpy.env.workspace = os.getcwd()
    arcpy.env.overwriteOutput = True
    arcpy.env.qualifiedFieldNames = False
    arcpy.CheckOutExtension("spatial") # Check spatial analyst license




    print "\nperforming spatial join of linework to grid... "

    # if a list of features is provided, merge them together; put in same place as first input file
    if not isinstance(stream_linework, basestring):
        print "merging: "
        for lw in stream_linework:
            print lw
        merged_linework = os.path.join(os.path.split(stream_linework[0])[0], 'input_linework.shp')
        arcpy.Merge_management(stream_linework, merged_linework)
        stream_linework = merged_linework

    arcpy.SpatialJoin_analysis(MFgrid, stream_linework,
                               stream_cells,
                               "JOIN_ONE_TO_MANY",
                               "KEEP_COMMON")

    print "\nDissolving river cells on cell number to isolate unique cells...\n"
    arcpy.Dissolve_management(stream_cells, stream_cells_dissolve, MFgrid_node_attribute)

    print "Exploding new stream linework to grid cells using Intersect and Multipart to Singlepart..."
    arcpy.Intersect_analysis([stream_cells_dissolve, stream_linework], "tmp_intersect.shp")
    arcpy.MultipartToSinglepart_management("tmp_intersect.shp", stream_fragments)

    # make a new feature layer from exploded shape file
    arcpy.MakeFeatureLayer_management(stream_fragments, 'stream_fragments')

    if not from_scratch:
        print "Removing linework that overlaps with existing SFR cells..."
        arcpy.MakeFeatureLayer_management(existing_sfr_shp, 'existing_sfr_cells')
        # use "WITHIN"; "INTERSECT" was also deleting fragments that touched existing cells
        arcpy.SelectLayerByLocation_management('stream_fragments', "INTERSECT", 'existing_sfr_cells', "", "NEW_SELECTION")
        arcpy.DeleteFeatures_management('stream_fragments')
        #arcpy.CopyFeatures_management('stream_fragments', stream_fragments) # save layer to shapefile


    print "Adding in stream geometry..."
    #set up list and dictionary for fields, types, and associated commands
    fields = ('X_start', 'Y_start', 'X_end', 'Y_end', 'LengthFt')
    types = {'X_start': 'DOUBLE',
             'Y_start': 'DOUBLE',
             'X_end': 'DOUBLE',
             'Y_end': 'DOUBLE',
             'LengthFt': 'DOUBLE'}
    commands = {'X_start': "float(!SHAPE.firstpoint!.split()[0])",
                'Y_start': "float(!SHAPE.firstpoint!.split()[1])",
                'X_end': "float(!SHAPE.lastpoint!.split()[0])",
                'Y_end': "float(!SHAPE.lastpoint!.split()[1])",
                'LengthFt': "float(!SHAPE.length!)"}

    #add fields for start, end, and length
    for fld in fields:
        arcpy.AddField_management('stream_fragments', fld, types[fld])


    #calculate the fields
    for fld in fields:
        print "\tcalculating %s(s)..." % (fld)
        arcpy.CalculateField_management('stream_fragments', fld, commands[fld], "PYTHON")

    ofp.write('\n' + 25*'#' + '\nRemoving reaches with lengths less than or equal to %s...\n' % reach_cutoff)
    print "\nRemoving reaches with lengths less than or equal to %s..." % reach_cutoff
    table = arcpy.UpdateCursor('stream_fragments')
    count = 0
    for reaches in table:
        if reaches.getValue('LengthFt') <= reach_cutoff:
            print "cellnum: %d" % (reaches.getValue(MFgrid_node_attribute)),
            ofp.write("cellnum: %d" % (reaches.getValue(MFgrid_node_attribute)))
            table.deleteRow(reaches)
            count += 1
    print "\nremoved %s reaches with lengths <= %s\n" % (count, reach_cutoff)
    ofp.write("removed %s reaches with lengths <= %s\n" % (count, reach_cutoff))


    # create a unique ID for each new stream fragment using FID and number of existing reaches
    arcpy.AddField_management('stream_fragments', "FragID", "LONG")
    arcpy.CalculateField_management('stream_fragments', "FragID", "!FID! + {0:d}".format(nreaches), "PYTHON")
    #arcpy.CopyFeatures_management('stream_fragments', stream_fragments) # save layer to shapefile

if preproc:
    preprocessing(stream_linework)

print "creating dataframe of new SFR cell information..."

stream_fragments_shp = fiona.open(stream_fragments, 'r')
new_streamcells_df = False
new_streamcells = {}
knt = 0
nfrags = len(stream_fragments_shp)
for line in stream_fragments_shp:

    attributes = line['properties']
    new_streamcells[attributes['FragID']] = attributes
    knt += 1
    print '\r{:d}%'.format(100*knt/nfrags),

# convert to pandas dataframe, join in centroids, sort by FID
new_streamcells_df = pd.DataFrame.from_dict(new_streamcells, orient='index')
new_streamcells_df.index = new_streamcells_df['FragID']


print "\nestablishing segments"
tol = reach_cutoff/np.sqrt(2) # tol below considers x and y distances individually (not the diagonal)
new_streamcells_df['Segment'] = np.zeros((len(new_streamcells_df)))
new_streamcells_df['Reach'] = np.zeros((len(new_streamcells_df)))

# if adding to existing SFR dataset, start segment number after last existing seg
if not from_scratch:
    seg = np.max(existing_SFR_df['segment'])
    last_existing_seg = seg
else:
    seg = 0
    last_existing_seg = 0

n_newcells = len(new_streamcells_df)

# make vectors of start and end locations for segments
# not sure if this is faster than simply indexing/slicing in pandas or not
xstarts, ystarts = new_streamcells_df['X_start'].values, new_streamcells_df['Y_start'].values
xends, yends = new_streamcells_df['X_end'].values, new_streamcells_df['Y_end'].values

assignedFragIDs = [] # keep track of frag IDs that have already been considered to minimize tests with np.where

for FragID in new_streamcells_df.index:
    '''
    This loop is pretty slow!, probably not ideal for large stream networks (>100 segments)
    '''

    # ignore fragments that have already been assigned to segments
    if FragID in assignedFragIDs:
        continue

    # check if there is anything upstream
    xstart, ystart = new_streamcells_df.ix[FragID, ['X_start', 'Y_start']].values

    up = np.where((np.abs(xends - xstart) < tol ) & (np.abs(yends - ystart) < tol))[0]

    # test for upstream divergence (if one is found, make a new segment)
    up_div = np.where((np.abs(xstarts - xstart) < tol ) & (np.abs(ystarts - ystart) < tol))[0]

    # if there is nothing upstream, or if the current fragment is at a confluence, record segment and then look downstream
    at_end = False
    reach = 1
    if len(up) == 0 or len(up) > 1 or len(up_div) > 1:
        seg += 1
        print "\r{:d}".format(seg),

        if seg == 49:
            j=2

        if seg > 20 and 20 not in new_streamcells_df.Segment.values:
            j=2

        #print "\nsegment: {}, reaches:".format(seg),
        new_streamcells_df.ix[FragID, 'Segment'] = seg
        new_streamcells_df.ix[FragID, 'Reach'] = reach
        assignedFragIDs.append(FragID)

        while not at_end:
            reach += 1
            # find downstream reach, if none, stop
            xend, yend = new_streamcells_df.ix[FragID, ['X_end', 'Y_end']].values

            # get indices of all reaches that have starts within tol of end
            down = np.where((np.abs(xstarts - xend) < tol) &
                            (np.abs(ystarts - yend) < tol))[0]

            # if there are no reaches (outlet) or multiple (divergence), at segment end
            if len(down) != 1:
                at_end = True
                break

            # get index of downstream reach in streamcells dataframe
            down = new_streamcells_df.index[down][0]

            # first check to see if the next reach down is below a confluence, stop if it is
            xstart_d, ystart_d = new_streamcells_df.ix[down, ['X_start', 'Y_start']].values
            up_of_down = np.where((np.abs(xends - xstart_d) < tol) &
                                  (np.abs(yends - ystart_d) < tol))[0]
            if len(up_of_down) > 1:
                at_end = True
                break

            # then check to see if the next reach down is an outlet, record segment, then stop if it is
            # this is needed because if the current reach is an outlet, it doesn't get assigned to a segment
            xend_d, yend_d = new_streamcells_df.ix[down, ['X_end', 'Y_end']].values
            down_of_down = np.where((np.abs(xstarts - xend_d) < tol) &
                                    (np.abs(ystarts - yend_d) < tol))[0]

            if len(down_of_down) == 0:
                new_streamcells_df.ix[down, 'Segment'] = seg
                new_streamcells_df.ix[down, 'Reach'] = reach
                assignedFragIDs.append(down) # record fragments that were already assigned to segments
                #print "{} ".format(reach),
                at_end = True
                break

            elif down in assignedFragIDs:
                break
            # otherwise record segment number and keep going
            else:
                new_streamcells_df.ix[down, 'Segment'] = seg
                new_streamcells_df.ix[down, 'Reach'] = reach
                assignedFragIDs.append(down)
                FragID = down
                #print "{} ".format(reach),
    else:
        continue

# segments were established in above loop based on confluences (where multiple line fragments shared same end or start)
# and based on line fragment ends that were not touching other line fragments
# however, with linework from multiple data sources, segment ends might not touch
# for example, if a segment from one data source makes a "T" with a segment from another datasource, at a distance
# greater than the tight tolerance used in establishing segment above, a separate segment will be made for the trib,
# but the main stem will not be subdivided.
# Go through and look for connecting segments with a looser "routing" tolerance
print 'establishing more segments by subdividing at confluences involving interior reaches...'
new_streamcells_df['Outseg'] = np.zeros((len(new_streamcells_df)))

segments = np.unique(new_streamcells_df.Segment).astype(int)
nsegs = np.max(segments)
for seg in segments:

    # get info for just current segment, sort on reach, and then record end coordinates
    df = new_streamcells_df[new_streamcells_df['Segment'] == seg].sort('Reach')
    xend, yend = df.iloc[-1]['X_end'], df.iloc[-1]['Y_end']

    # find closest segment, excluding current
    othersegs = new_streamcells_df.ix[new_streamcells_df.Segment != seg]

    # apparently when working with a dataframe that is a view on another dataframe (and not a copy)
    # np.argmin will return the index number of the larger dataframe (not the position within the view)
    out_idx = np.argmin(np.sqrt((othersegs['X_start'] - xend)**2 +
                            (othersegs['Y_start'] - yend)**2))

    outreach = new_streamcells_df.ix[out_idx, 'Reach']

    dist = np.sqrt((new_streamcells_df.ix[out_idx, 'X_start'] - xend)**2 +
                   (new_streamcells_df.ix[out_idx, 'Y_start'] - yend)**2)

    # check if the closest segment is within the routing tolerance
    if dist < routing_tol:

        # if it is, set it as the outseg of the current segment
        outsegnum = new_streamcells_df.ix[out_idx, "Segment"]
        new_streamcells_df.Outseg[new_streamcells_df.Segment == seg] = new_streamcells_df.ix[out_idx, "Segment"]

        # now check if the closest reach is reach 1; if it isn't we need to break up that segment
        if outreach != 1:

            # assign new segment to reaches in the outseg upstream of the outreach
            new_streamcells_df.Segment[(new_streamcells_df['Segment'] == outsegnum) &
                                       (new_streamcells_df['Reach'] < outreach)] = nsegs + 1
            nsegs += 1
            print nsegs

            # renumber subsequent reaches in outseg by subtracting (outreach - 1)
            new_streamcells_df.Reach[(new_streamcells_df['Segment'] == outsegnum) &
                                     (new_streamcells_df['Reach'] >= outreach)] -= (outreach - 1)


new_streamcells_df.to_csv('out_segments.csv')

new_streamcells_df = pd.read_csv('out_segments.csv')

print "\nrouting new segments..."
tol = 1.0
new_segs = map(int, np.unique(new_streamcells_df[new_streamcells_df['Segment'] > 0].Segment))
new_streamcells_df['Elevmin'] = np.zeros((len(new_streamcells_df)))
new_streamcells_df['Elevmax'] = np.zeros((len(new_streamcells_df)))

# load active area boundary and make into shapely LineString object
b = fiona.open(MFdomain)
MFbound = LineString(b.next()['geometry']['coordinates'])

for seg in new_segs:

    if seg == 20:
        j=2
    # get info for just current segment, sort on reach, and then record end coordinates
    data = new_streamcells_df[new_streamcells_df['Segment'] == seg].sort('Reach')

    # find outseg with matching start coordinates and assign to outseg column
    out = np.where((np.abs(new_streamcells_df['X_start'] - data.iloc[-1]['X_end']) < tol) &
                    (np.abs(new_streamcells_df['Y_start'] - data.iloc[-1]['Y_end']) < tol))[0]

    if len(out) > 0:

        out_idx = new_streamcells_df.index[out][0]
        new_streamcells_df.Outseg[new_streamcells_df['Segment'] == seg] = new_streamcells_df.ix[out_idx, "Segment"]

    # if no outseg within new stream segments, search for nearby existing SFR cells
    elif len(out) == 0 and not from_scratch:
        closest = np.argmin(np.sqrt((existing_SFR_df['X'] - data.iloc[-1]['X_end'])**2 + (existing_SFR_df['Y'] - data.iloc[-1]['Y_end'])**2))
        dist = np.sqrt((existing_SFR_df.iloc[closest]['X'] - data.iloc[-1]['X_end'])**2 + (existing_SFR_df.iloc[closest]['Y'] - data.iloc[-1]['Y_end'])**2)

        if dist <= max_distance_from_SFR:
            # if closest existing SFR segment is first-order (no upsegs), route to next downstream segment
            # this helps avoid problems of artificially routing water to first order streams that are dry
            if existing_SFR_df.iloc[closest]['upseg'] == 0:
                outseg = existing_SFR_df.iloc[closest]['outseg']
            else:
                outseg = existing_SFR_df.iloc[closest]['segment']
            new_streamcells_df.Outseg[new_streamcells_df['Segment'] == seg] = outseg
            # set minimum elevation for segment from existing SFR
            out_elevation = existing_SFR_df.iloc[closest][existing_sfr_elevation_attribute]
            new_streamcells_df.Elevmin[new_streamcells_df['Segment'] == seg] = out_elevation

    # or, test for proximity to grid boundary (outlet condition)
    else:
        # make a buffer of max_distance_to_model_boundary around segment endpoint
        endpoint = Point(data.iloc[-1]['X_end'], data.iloc[-1]['Y_end'])
        endpoint_buff = endpoint.buffer(max_distance_to_model_boundary)

        # test for intersection with active domain bound
        if endpoint_buff.intersects(MFbound):
            new_streamcells_df.Outseg[new_streamcells_df['Segment'] == seg] = 0
            new_streamcells_df.Elevmin[new_streamcells_df['Segment'] == seg] = 0
        else:
            continue

# for now, leave in unrouted segments
'''
# add FragID field identifying unique river_explode segments if one doesn't exist
arcpy.MakeFeatureLayer_management(stream_fragments, "river_explode")
Fields = [f.name.lower() for f in arcpy.ListFields("river_explode")]
if "fragid" not in Fields:
    arcpy.AddField_management(stream_fragments, "FragID", "LONG")
    arcpy.CalculateField_management(stream_fragments, "FragID", "!FID!", "PYTHON")

This elevations section no longer needed, since Elevations.py can assign elevations using Mat1 and 2

print "\nassigning elevations from DEM..."

print "Intersecting \n{} \nfragments with \n{}...".format(stream_fragments, DEM)
# convert end vertices of river_explode Fragments to points
arcpy.FeatureVerticesToPoints_management(stream_fragments, stream_fragments_points)



# extract DEM values at locations of points
arcpy.sa.ExtractMultiValuesToPoints(stream_fragments_points, [[DEM, "DEM_elev"]])

# read points back into pandas dataframe
new_streams_points_df = shp2df(stream_fragments_points)

# assign DEM elevations to each FragID following intersect_DEM preprocessing
new_streamcells_df['DEM_elev'] = np.zeros(len(new_streamcells_df))
for FragID in new_streamcells_df.FragID:
    # assign minimum elevation from all vertices to fragment
    ind = new_streamcells_df[new_streamcells_df.FragID == FragID].index[0]
    new_streamcells_df.ix[ind, "DEM_elev"] = \
        np.min(new_streams_points_df[new_streams_points_df['FragID'] == FragID].DEM_elev) * DEM_z_conversion_factor

# now determine min/max elevs for each segment
for seg in new_segs:

    outseg = int(new_streamcells_df.Outseg[seg])
    upsegs = np.unique(new_streamcells_df.Segment[new_streamcells_df.Outseg == seg])

    elevations = list(new_streamcells_df[new_streamcells_df['Segment'] == seg].sort('Reach').DEM_elev)
    outseg_elevations = new_streamcells_df[new_streamcells_df['Segment'] == outseg].sort('Reach').DEM_elev

    # set elevmax/elevmin for headwater segments
    if len(upsegs) == 0:
        # if outseg is existing SFR, and existing Elevmin <= starting elevation, set Elevmax to starting elevation
        if outseg <= last_existing_seg:
            if new_streamcells_df[new_streamcells_df['Segment'] == seg].Elevmin.iloc[0] <= elevations[0]:
                new_streamcells_df.Elevmax[new_streamcells_df['Segment'] == seg] = elevations[0]
                continue
        # otherwise set Elevmin from minimum; set elevmax equal to Elevmin
            else:
                new_streamcells_df.Elevmin[new_streamcells_df['Segment'] == seg] = np.min(elevations)
                new_streamcells_df.Elevmax[new_streamcells_df['Segment'] == seg] = np.min(elevations)

    # for non-headwater streams
    else:
        # if outseg is existing SFR, and existing Elevmin <= minimum upseg elevation, set Elevmax to minimum upseg elevation
        upseg_mins = []
        for upseg in upsegs:
            upseg_mins.append(np.min(new_streamcells_df[new_streamcells_df['Segment'] == upseg].sort('Reach').DEM_elev))
        upseg_min = np.min(upseg_mins)
        new_streamcells_df.Elevmax[new_streamcells_df['Segment'] == seg] = upseg_min
        if outseg <= last_existing_seg:
            if new_streamcells_df[new_streamcells_df['Segment'] == seg].Elevmin.iloc[0] <= upseg_min:
                continue
            # unless upseg minimum is below outseg, then set equal to outseg (flat segment)
            else:
                new_streamcells_df.Elevmax[new_streamcells_df['Segment'] == seg] = \
                    new_streamcells_df[new_streamcells_df['Segment'] == seg].Elevmin
        # otherwise segment is an interior new segment; set elevmax from upseg minimum, and elevmin from outseg max
        else:
            outseg_max = np.max(outseg_elevations)
            new_streamcells_df.Elevmax[new_streamcells_df['Segment'] == seg] = upseg_min
            new_streamcells_df.Elevmin[new_streamcells_df['Segment'] == seg] = outseg_max
            if upseg_min < outseg_max:
                print "Warning, segment {} has backwards elevations".format(seg)
                ofp.write("Warning, segment {} has backwards elevations".format(seg))
                ofp.write("outseg: {}, maxelev: {}; upsegs: {}, minelev: {}\n".format(outseg, outseg_max, list(upsegs), upseg_min))

new_streamcells_df['DEM_elev_smoothed'] = np.zeros(len(new_streamcells_df))
seg_data = new_streamcells_df[new_streamcells_df['Segment'] == seg].sort('Reach')

# now smooth interior elevations
# first build dictionaries of raw elevations, reach lengths, and segment min/max elevations
lengths = defaultdict(list)
raw_elevations = defaultdict(list)
min_max_elevations = defaultdict(list)
for seg in new_segs:

    lengths[seg] = list(new_streamcells_df[new_streamcells_df['Segment'] == seg].sort('Reach').LengthFt)
    raw_elevations[seg] = list(new_streamcells_df[new_streamcells_df['Segment'] == seg].sort('Reach').DEM_elev)
    min_max_elevations[seg] = [np.min(new_streamcells_df[new_streamcells_df['Segment'] == seg].sort('Reach').Elevmin),
                               np.max(new_streamcells_df[new_streamcells_df['Segment'] == seg].sort('Reach').Elevmax)]

# smooth interior elevations
elevations_sm, slopes = sm.connect_downhill(new_segs, lengths, raw_elevations, min_max_elevations, ofp)
elevations_sm = raw_elevations
'''
# read in elevations from dis file
DX, DY, NLAY, NROW, NCOL, i = disutil.read_meta_data(DISfile)
topdata, i = disutil.read_nrow_ncol_vals(DISfile, NROW, NCOL, np.float, i)

print "writing {}...".format(Mat1_updated)
# lookup row column from node number
X, Y = np.meshgrid(np.arange(1, NCOL+1), np.arange(1, NROW+1))
cr = np.vstack([X.ravel(), Y.ravel()])

Mat1_new_dict = {}
# build dict of new_segs (parallel structure to exisiting SFR dataframe)
for cell in new_streamcells_df.index:
    if new_streamcells_df.ix[cell, 'Segment'] > 0:
        node = new_streamcells_df.ix[cell, MFgrid_node_attribute]
        c, r = cr[:, node-1] # headache! -1 needed to translate from base=1 to base=0

        #sbtop = topdata[r, c] # get streambed top from model top
        seg, reach = int(new_streamcells_df.ix[cell, 'Segment']), int(new_streamcells_df.ix[cell, 'Reach'])
        #sbtop = elevations_sm[seg][reach-1]
        #slope = slopes[seg][reach-1]
        sbtop = 0
        slope = 0

        reach_props = {'row': r,
                       'column': c,
                       'layer': 1,
                       'stage': sbtop + 1,
                       'top_streambed': sbtop,
                       'reach': reach,
                       'segment': seg,
                       'width_in_cell': width_in_cell,
                       'length_in_cell': new_streamcells_df.ix[cell, 'LengthFt'],
                       'bed_K': bed_K,
                       'bed_thickness': bed_thickness,
                       'bed_slope': slope,
                       'bed_roughness': bed_roughness}
        Mat1_new_dict[new_streamcells_df.ix[cell, 'FragID']] = reach_props

# append to Mat1 and write out a new one
Mat1_new = pd.DataFrame.from_dict(Mat1_new_dict, orient='index')
if not from_scratch:
    Mat1 = Mat1.append(Mat1_new)
    Mat1 = Mat1.sort(['segment', 'reach'])
    Mat1.to_csv(Mat1_updated, index=False)
else:
    Mat1_new.to_csv(Mat1_updated, index=False)

print "writing {}...".format(Mat2_updated)
new_streamcells_df.to_csv('out_segments2.csv')

Mat2_new_dict = {}
for seg in new_segs:
    if new_streamcells_df.Segment[new_streamcells_df['Segment'] == seg].iloc[0] > 0:
        #segment	icalc	outseg	iupseg	iprior	nstrpts	flow	runoff	etsw	pptsw	roughch	roughbk	cdepth	fdepth	awdth	bwdth

        seg_props = {'segment': int(seg),
                     'icalc': 1,
                     'outseg': int(new_streamcells_df.Outseg[new_streamcells_df['Segment'] == seg].iloc[0]),
                     'iupseg': 0,
                     'iprior': 0,
                     'nstrpts': 0,
                     'flow': 0,
                     'runoff': 0,
                     'etsw': 0,
                     'pptsw': 0,
                     'roughch': bed_roughness,
                     'roughbk': 0,
                     'cdepth': 0,
                     'fdepth': 0,
                     'awdth': 0,
                     'bwdth': 0}

        Mat2_new_dict[seg] = seg_props

Mat2_new = pd.DataFrame.from_dict(Mat2_new_dict, orient='index')
if not from_scratch:
    Mat2.index = Mat2.index + 1 # get rid of zero-indexing
    Mat2 = Mat2.append(Mat2_new)
    Mat2.to_csv(Mat2_updated, index=False)
else:
    Mat2_new.to_csv(Mat2_updated, index=False)

# if adding to existing SFR dataset, combine new info with existing
# otherwise it is easier just to run buildSFRshapefile2 with new Mat 1 and Mat2

if not from_scratch:
    print "making new shapefiles {} containing existing and new streamcells and linework...".format(stream_cells_updated)
    ofp.write("making new shapefile {} containing existing and new streamcells and linework...\n".format(stream_cells_updated))
    temp_shp = os.path.join(os.path.split(stream_cells_updated)[0], 'streamcells_tmp.shp')
    arcpy.Union_analysis([existing_sfr_shp, stream_cells_dissolve], temp_shp, "ALL")

    # open shapefile back up and consolidate node numbers
    knt = 0
    nSFRcells = len(Mat1_new)

    def remap_node_number(node_atrb1, node_atrb2, inshp, outshp, type, prj, error_logfile):
        with fiona.collection(inshp, "r") as input:

            schema = {'geometry': type,
                      'properties': {'node': 'int'}}

            with fiona.collection(outshp, "w", "ESRI Shapefile", schema) as output:
                for node in input:
                    node_num1 = node['properties'][node_atrb1]
                    node_num2 = node['properties'][node_atrb2]

                    # pick a node number
                    if node_num1 == 0 or node_num1 == node_num2:
                        node_num = node_num2
                    elif node_num2 == 0:
                        node_num = node_num1
                    else:
                        error_logfile.write("Warning! node number conflict. MFgrid node number: {}, "
                                  "Existing SFR node number: {}\n".format(node_num1, node_num2))

                    print "\rnode {:d}".format(node_num),

                    output.write({'properties': {'node': node_num},
                                  'geometry': mapping(shape(node['geometry']))})
        # copy over prj file
        shutil.copyfile(prj, "{}.prj".format(outshp[:-4]))

    # consolidate node numbers for cells
    remap_node_number(MFgrid_node_attribute, existing_sfr_node_attribute, temp_shp, stream_cells_updated, 'Polygon',
                      "{}.prj".format(existing_sfr_shp[:-4]), ofp)

    # do same for linework
    arcpy.Merge_management([existing_SFR_linework, stream_fragments], temp_shp)
    remap_node_number(MFgrid_node_attribute, existing_sfr_node_attribute, temp_shp, SFR_linework_updated, 'LineString',
                      "{}.prj".format(existing_sfr_shp[:-4]), ofp)



    ofp.close()

    print "Done"


