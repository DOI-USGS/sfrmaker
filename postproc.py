__author__ = 'aleaf'
'''

'''
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from shapely.geometry import Polygon, LineString
import flopy
from sfr_plots import SFRshapefile
import GISio
import GISops


# Functions
def header(infile):
    ofp = open(infile)
    knt=0
    while True:
        try:
            line = ofp.readline().strip().split()
            int(line[0])
            break
        except:
            knt +=1
            continue
    return knt


class SFRdata(object):

    def __init__(self, sfrobject=None, Mat1=None, Mat2=None, sfr=None, node_column=False,
                 mfpath='', mfnam=None, mfdis=None, landsurfacefile=None, GIS_mult=1, to_meters_mult=0.3048,
                 Mat2_out=None, xll=0.0, yll=0.0, prj=None, proj4=None, epsg=None,
                 minimum_slope=1e-4, streamflow_file=None):
        """
        base object class for SFR information in the SFRmaker postproc module.

        Parameters
        ----------
        sfrobject: SFRdata instance
            Instantiates SFRdata with attributes from another SFRdata instance

        Mat1: dataframe or filepath


        GIS_mult: float
            Multiplier to go from model units to GIS units

        """
        # if an sfr object is supplied, copy all of its attributes over
        # need to clean up the init method!!
        self.outpath = ''
        self.elevs = None

        if sfrobject is not None:
            self.__dict__ = sfrobject.__dict__.copy()

        else:
            if Mat1 is not None:
                # read Mats 1 and 2
                # allow the dataframes to be passed directly instead of reading them in from csvs
                if isinstance(Mat1, pd.DataFrame):
                    self.m1 = Mat1
                    self.m2 = Mat2
                else:
                    self.Mat1 = Mat1
                    self.Mat2 = Mat2
                    self.m1 = pd.read_csv(Mat1)
                    self.m2 = pd.read_csv(Mat2)
                    self.outpath = os.path.split(Mat1)[0]
                self.m1.sort(['segment', 'reach'], inplace=True)
                self.m2.sort('segment', inplace=True)
                self.m2.index = map(int, self.m2.segment.values)
            elif sfr is not None:
                self.read_sfr_package()
            else:
                # is this the right kind of error?
                raise AssertionError("Please specify either Mat1 and Mat2 files or an SFR package file.")

            # Discretization
            self.elevs_by_cellnum = {} # dictionary to store the model top elevation by cellnumber
            self.cell_geometries = {}

            if mfdis is not None:
                self.mfpath = mfpath
                self.basename = mfdis[:-4]
                if mfnam is None:
                    mfnam = self.basename + '.nam'
                self.mfdis = os.path.join(self.mfpath, mfdis)
                self.mfnam = os.path.join(self.mfpath, mfnam)

                # node numbers for cells with SFR
                if not node_column:
                    self.node_column = 'node'
                    self.gridtype = 'structured'
                    self.read_dis2()
                    self.m1['node'] = (self.dis.ncol * (self.m1['row'] - 1) + self.m1['column']).astype('int')
                else:
                    self.m1.rename(columns={node_column, 'node'}, inplace=True)
                    self.node_column = 'node' # so that some code doesn't break for the time being
            else:
                self.mfpath = ''
                self.mfnam = None
                self.mfdis = None
                self.basename = 'MF'

            self.node_attribute = 'node' # compatibility with sfr_plots and sfr_classes

            self.landsurfacefile = landsurfacefile
            self.landsurface = None
            self.xll = xll
            self.yll = yll
            self.to_m = to_meters_mult # multiplier to convert model units to meters (used for width estimation)
            self.to_km = self.to_m / 1000.0
            self.GIS_mult = GIS_mult
            self.minimum_slope = minimum_slope

            # coordinate projection of model grid
            self.prj = prj # projection file
            self.proj4 = proj4
            self.epsg = epsg
            if proj4 is None and prj is not None:
                self.proj4 = GISio.get_proj4(prj)

            # streamflow file for plotting SFR results
            if mfnam is not None and streamflow_file is None:
                streamflow_file = mfnam[:-4] + '_streamflow.dat'
            self.streamflow_file = streamflow_file

            self.segments = sorted(np.unique(self.m1.segment))
            self.Mat1_out = Mat1
            self.Mat2_out = Mat2_out

            # assign upstream segments to Mat2
            self.m2['upsegs'] = [self.m2.segment[self.m2.outseg == s].tolist() for s in self.segments]

            # check for circular routing
            c = [s for s in self.m2.segment if s == self.m2.outseg[s]]
            if len(c) > 0:
                raise ValueError('Warning! Circular routing in segments {}.\n'
                                 'Fix manually in Mat2 before continuing'.format(', '.join(map(str, c))))
    @ property
    def shared_cells(self):
        return np.unique(self.m1.ix[self.m1.node.duplicated(), 'node'])

    def read_sfr_package(self):
        """method to read in SFR file
        """
        pass

    def read_dis2(self, mfdis=None, mfnam=None):
        """read in model grid information using flopy
        """
        if mfdis is not None:
            self.mfdis = mfdis
            self.mfnam = mfnam
            self.mfpath = os.path.split(self.mfdis)[0]

        print self.mfdis
        self.m = flopy.modflow.Modflow(model_ws=self.mfpath)
        self.nf = flopy.utils.mfreadnam.parsenamefile(self.mfnam, {})
        self.dis = flopy.modflow.ModflowDis.load(self.mfdis, self.m, self.nf)
        self.elevs = np.zeros((self.dis.nlay + 1, self.dis.nrow, self.dis.ncol))
        self.elevs[0, :, :] = self.dis.top.array
        self.elevs[1:, :, :] = self.dis.botm.array

        # check if there Quasi-3D confining beds
        if np.sum(self.dis.laycbd.array) > 0:
            print 'Quasi-3D layering found, skipping confining beds...'
            layer_ind = 0
            for l, laycbd in enumerate(self.dis.laycbd.array):
                self.elevs[l + 1, :, :] = self.dis.botm.array[layer_ind, :, :]
                if laycbd == 1:
                    print '\tbetween layers {} and {}'.format(l+1, l+2)
                    layer_ind += 1
                layer_ind += 1
        else:
            self.elevs[1:, :, :] = self.dis.botm.array

        # make dictionary of model top elevations by cellnum
        for c in range(self.dis.ncol):
            for r in range(self.dis.nrow):
                cellnum = r * self.dis.ncol + c + 1
                self.elevs_by_cellnum[cellnum] = self.elevs[0, r, c]

    def read_DIS(self, MFdis):
        """read discretisation file using disutil
        """
        self.DX, self.DY, self.NLAY, self.NROW, self.NCOL, i = disutil.read_meta_data(MFdis)

        # get layer tops/bottoms
        self.layer_elevs = np.zeros((self.NLAY+1, self.NROW, self.NCOL))
        for c in range(self.NLAY + 1):
            tmp, i = disutil.read_nrow_ncol_vals(MFdis, self.NROW, self.NCOL, 'float', i)
            self.layer_elevs[c, :, :] = tmp

        # make dictionary of model top elevations by cellnum
        self.elevs_by_cellnum = {}
        for c in range(self.NCOL):
            for r in range(self.NROW):
                cellnum = r*self.NCOL + c + 1
                self.elevs_by_cellnum[cellnum] = self.layer_elevs[0, r, c]

    def get_cell_geometries(self):

        print 'computing cell geometries...'
        self.centroids = self.dis.get_node_coordinates()
        self.get_cell_centroids()

        for i, dfr in self.m1.iterrows():
            r, c = dfr.row - 1, dfr.column - 1
            cn = r * self.dis.ncol + c + 1
            cx, cy = dfr.centroids
            dx = self.dis.delr[c] * self.GIS_mult
            dy = self.dis.delc[r] * self.GIS_mult  # this is confusing, may be a bug in flopy

            # calculate vertices for parent cell
            x0 = cx - 0.5 * dx
            x1 = cx + 0.5 * dx
            y0 = cy - 0.5 * dy
            y1 = cy + 0.5 * dy

            self.cell_geometries[cn] = Polygon([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)])
        self.m1['geometry'] = [self.cell_geometries[cn] for cn in self.m1.node]

    def get_cell_centroids(self):
        """Adds column of cell centroids coordinates (tuples) to Mat1,
        from flopy DIS object.
        """

        self.centroids = self.dis.get_node_coordinates()

        centroids = []
        for i, r in self.m1.iterrows():
            r, c = r.row - 1, r.column - 1
            cx = self.centroids[1][c] * self.GIS_mult + self.xll
            cy = self.centroids[0][r] * self.GIS_mult + self.yll
            centroids.append((cx, cy))

        self.m1['centroids'] = centroids

    def map_outsegs(self):
        '''
        from Mat2, returns dataframe of all downstream segments (will not work with circular routing!)
        '''
        outsegsmap = pd.DataFrame(self.m2.outseg)
        outsegs = self.m2.outseg
        max_outseg = np.max(outsegsmap[outsegsmap.columns[-1]])
        knt = 2
        nsegs = len(self.m2)
        while max_outseg > 0:
            outsegsmap['outseg{}'.format(knt)] = [outsegs[s] if s > 0 and s < 999999 else 0
                                                    for s in outsegsmap[outsegsmap.columns[-1]]]
            max_outseg = np.max(outsegsmap[outsegsmap.columns[-1]].values)
            if max_outseg == 0:
                break
            knt +=1
            if knt > 200:
                # subset outsegs map to all rows outseg number > 0 at iteration 1000
                circular_segs = outsegsmap[outsegsmap[outsegsmap.columns[-1]] > 0]

                # only retain one instance of each outseg number at iteration 1000
                circular_segs.drop_duplicates(subset=outsegsmap.columns[-1], inplace=True)

                # cull the dataframe again to remove duplicate instances of routing circles
                circles = []
                duplicates = []
                for i in circular_segs.index:
                    repeat_start_ind = np.where(circular_segs.ix[i, :] ==
                                                circular_segs.ix[i, circular_segs.columns[-1:]]
                                                .values[0])[0][-2:][0]
                    circular_seq = circular_segs.ix[i, circular_segs.columns[repeat_start_ind:]].tolist()

                    if set(circular_seq) not in circles:
                        circles.append(set(circular_seq))
                    else:
                        duplicates.append(i)
                circular_segs.drop(duplicates, axis=0, inplace=True)

                rf = 'Circular_routing_outsegs_table.csv'
                circular_segs.to_csv(rf)
                return '{} instances of circular routing encountered! See {} for details.'\
                        .format(len(circular_segs), rf)

        self.outsegs = outsegsmap

        # create new column in Mat2 listing outlets associated with each segment
        self.m2['Outlet'] = [r[r != 0][-1] if len(r[r != 0]) > 0
                             else i for i, r in self.outsegs.iterrows()]

        # assign the outlets to each reach listed in Mat1
        self.m1['Outlet'] = [self.m2.ix[seg, 'Outlet'] for seg in self.m1.segment]

    def map_confluences(self, landsurfacefile=None, landsurface_column=None):

        if landsurfacefile is None and self.landsurfacefile is not None:
            landsurfacefile = self.landsurfacefile

        if hasattr(self, 'Elevations'):
            self.Elevations.map_confluences()
        else:
            self.Elevations = Elevations(sfrobject=self, landsurfacefile=landsurfacefile,
                                         landsurface_column=landsurface_column)
            self.Elevations.map_confluences()

        self.m1 = self.Elevations.m1
        self.m2 = self.Elevations.m2
        self.confluences = self.Elevations.confluences

    def consolidate_conductance(self, bedKmin=1e-8):
        """For model cells with multiple SFR reaches, shift all conductance to widest reach,
        by adjusting the length, and setting the lengths in all smaller collocated reaches to 1,
        and the K-values in these to bedKmin

        Parameters
        ----------

        bedKmin : float
            Hydraulic conductivity value to use for collocated SFR reaches in a model cell that are not the dominant
            (widest) reach. This is used to effectively set conductance in these reaches to 0, to avoid circulation
            of water between the collocated reaches.

        Returns
        -------
            Modifies the SFR Mat1 table dataframe attribute of the SFRdata object. A new SFR package file can be
            written from this table.

        Notes
        -----
            See the ConsolidateConductance notebook in the Notebooks folder.
        """
        # Calculate SFR conductance for each reach
        def cond(X):
            c = X['bed_K'] * X['width_in_cell'] * X['length_in_cell'] / X['bed_thickness']
            return c

        self.m1['Cond'] = self.m1.apply(cond, axis=1)

        #shared_cells = np.unique(self.m1.ix[self.m1.node.duplicated(), 'node'])

        # make a new column that designates whether a reach is dominant in each cell
        # dominant reaches include those not collocated with other reaches, and the longest collocated reach
        self.m1['Dominant'] = [True] * len(self.m1)

        for c in self.shared_cells:

            # select the collocated reaches for this cell
            df = self.m1[self.m1.node == c].sort('width_in_cell', ascending=False)

            # set all of these reaches except the largest to not Dominant
            self.m1.loc[df.index[1:], 'Dominant'] = False

        # Sum up the conductances for all of the collocated reaches
        # returns a series of conductance sums by model cell, put these into a new column in Mat1
        Cond_sums = self.m1[['node', 'Cond']].groupby('node').agg('sum').Cond
        self.m1['Cond_sum'] = [Cond_sums[c] for c in self.m1.node]

        # Calculate a new length for widest reaches, set length in secondary collocated reaches to 1
        # also set the K values in the secondary cells to bedKmin
        self.m1['old_length_in_cell'] = self.m1.length_in_cell

        def consolidate_lengths(X):
            if X['Dominant']:
                lnew = X['Cond_sum'] * X['bed_thickness'] / (X['bed_K'] * X['width_in_cell'])
            else:
                lnew = 1.0
            return lnew

        self.m1['length_in_cell'] = self.m1.apply(consolidate_lengths, axis=1)
        self.m1['bed_K'] = [r['bed_K'] if r['Dominant'] else bedKmin for i, r in self.m1.iterrows()]

    def smooth_segment_ends(self, landsurfacefile=None, landsurface_column=None,
                            report_file='smooth_segment_ends.txt'):
        if hasattr(self, 'Elevations'):
            self.Elevations.smooth_segment_ends(report_file=report_file)
        else:
            self.Elevations = Elevations(sfrobject=self, landsurfacefile=landsurfacefile,
                                         landsurface_column=landsurface_column)
            self.Elevations.smooth_segment_ends(report_file=report_file)
        self.m2 = self.Elevations.m2

    def smooth_interior_elevations(self, landsurfacefile=None, landsurface_column=None,
                                   report_file='smooth_segment_interiors.txt'):
        """Allow elevations smoothing to be run on SFRdata instance (Composition)
        """
        if hasattr(self, 'Elevations'):
            self.Elevations.smooth_segment_interiors(report_file=report_file)
        else:
            self.Elevations = Elevations(sfrobject=self, landsurfacefile=landsurfacefile,
                                         landsurface_column=landsurface_column)
            self.Elevations.smooth_segment_interiors(report_file=report_file)
        self.m1 = self.Elevations.m1

    def reset_model_top_2streambed(self, minimum_thickness=1, outdisfile=None, outsummary=None):
        if hasattr(self, 'Elevations'):
            self.Elevations.reset_model_top_2streambed(minimum_thickness=minimum_thickness,
                                                       outdisfile=outdisfile, outsummary=outsummary)
        else:
            self.Elevations = Elevations(sfrobject=self)
            self.Elevations.reset_model_top_2streambed(minimum_thickness=minimum_thickness,
                                                       outdisfile=outdisfile, outsummary=outsummary)
        self.__dict__ = self.Elevations.__dict__.copy()

    def incorporate_field_elevations(self, shpfile, elevs_field, distance_tol):
        if hasattr(self, 'Elevations'):
            self.Elevations.incorporate_field_elevations(shpfile=shpfile, elevs_field=elevs_field,
                                                         distance_tol=distance_tol)
        else:
            self.Elevations = Elevations(sfrobject=self)
            self.Elevations.incorporate_field_elevations(shpfile=shpfile, elevs_field=elevs_field,
                                                         distance_tol=distance_tol)
        self.__dict__ = self.Elevations.__dict__.copy()

    def plot_stream_profiles(self,  outpdf='complete_profiles.pdf', units='ft', add_profiles={}):
        """Runs the Outsegs.plot_outsegs() method
        """
        if hasattr(self, 'Outsegs'):
            self.Outsegs.plot_outsegs(outpdf=outpdf, units=units, add_profiles=add_profiles)
        else:
            self.Outsegs = Outsegs(sfrobject=self)
            self.Outsegs.plot_outsegs(outpdf=outpdf, units=units, add_profiles=add_profiles)
        self.__dict__ = self.Outsegs.__dict__.copy()

    def plot_routing(self, outpdf='routing.pdf'):
        """Runs the Outsegs.plot_routing() method
        """
        if hasattr(self, 'Outsegs'):
            self.Outsegs.plot_routing(outpdf=outpdf)
        else:
            self.Outsegs = Outsegs(sfrobject=self)
            self.Outsegs.plot_routing(outpdf=outpdf)
        self.__dict__ = self.Outsegs.__dict__.copy()

    def estimate_stream_widths(self, Mat2_out=None):
        self.Widths = Widths(sfrobject=self)
        self.Widths.estimate_from_arbolate()
        self.__dict__ = self.Widths.__dict__.copy()

    def renumber_sfr_cells_from_polygons(self, intersect_df=None, intersect_shapefile=None, intersect_prj=None,
                                         sfr_shapefile=None,
                                         node_attribute=None, GIS_mult=None):
        """
        Subdivides SFR segments where they intersect polygon features supplied in a dataframe or a shapefile.
        Contiguous sequences of SFR reaches (within a segment)
        Parameters
        ----------
        intersect_df : DataFrame
            Contains a 'geometry' column of shapely polygon objects that intersect SFR cells. Contiguous sequences
            of SFR reaches intersecting a polygon will be made into new segments.
        :param intersect_shapefile:
        :param sfr_shapefile:
        :param node_attribute:
        :param GIS_mult:
        :return:
        """

        self.Spatial = Spatial(sfrobject=self)
        dfi = self.Spatial.intersect_with_SFR_cells(intersect_df=intersect_df, intersect_shapefile=intersect_shapefile,
                                                    intersect_prj=intersect_prj,
                                                    sfr_shapefile=sfr_shapefile,
                                                    node_attribute=node_attribute, GIS_mult=GIS_mult)

        '''
        sfr_cells = pd.read_csv('waterbodies_intersected.csv').sfr_cells.tolist()
        sfrc = []
        for c in sfr_cells:
            try:
                sfrc.append(map(int, c.replace('[','').replace(']','').split(',')))
            except:
                sfrc.append([])
        '''
        self.Segments = Segments(sfrobject=self)
        self.Segments.renumber_SFR_cells(dfi['sfr_cells'].tolist())
        #self.Segments.renumber_SFR_cells(sfrc)
        self.__dict__ = self.Segments.__dict__.copy()

    def renumber_sfr_cells_from_nodes(self, nodes_list):
        self.Segments = Segments(sfrobject=self)
        self.Segments.renumber_SFR_cells(nodes_list)
        self.__dict__ = self.Segments.__dict__.copy()

    def run_diagnostics(self, model_domain=None):
        """Run diagnostic suite on Mat1 and Mat2, including:
        * segment numbering
        * circular routing
        * collocated SFR reaches with non-zero conductances
        * elevations that increase in downstream direction,
          or that are inconsistent with layer tops/bottoms
        * check for outlets to stream network that

        Parameters:


        """
        from diagnostics import diagnostics
        self.diagnostics = diagnostics(sfrobject=self)
        self.diagnostics.check_numbering()
        self.diagnostics.check_routing()
        self.diagnostics.check_overlapping()
        self.diagnostics.check_elevations()
        self.diagnostics.check_outlets(model_domain=model_domain)

    def segment_reach2linework_shapefile(self, lines_shapefile, new_lines_shapefile=None,
                                         iterations=2):
        """Assigns segment and reach information to a shapefile of SFR reach linework,
        using intersections (coincident starts/ends) with neighboring SFR reaches in the case of
        model cells that have multiple co-located SFR reaches. This approach is more difficult
        in instances of multiple adjacent model cells, each with collocated SFR reaches, because
        the segment and reach information for geometries in the neighboring cells is also unknown.
        The method can overcome this by iterating. Right now the number of iterations is set to 2,
        which appeared to work fine for a small to medium -sized SFR package with few instances of this
        problem. But the method may not work, or more iterations may be required, for SFR packages in large
        models with large cell sizes (and therefore likely an increased number of cells with multiple SFR reaches).

        The best approach to this problem is to simply write a linework shapefile with segment and reach information
        when Mat1 is written for the first time. This will be implemented soon.

        Parameters
        ----------
        lines_shapefile: string
            Path to shapefile of SFR linework exploded by cell (e.g. one linework geometry for each SFR reach).

        new_lines_shapefile: string
            Path to new linework shapefile that will be written

        iterations: int
            With each iteration, the method will loop through the model cells containing multiple SFR reaches,
            and attempt to assign geometries using intersections with neighboring cells where the geometries
            are associated with a segment and reach. On this first iteration, the only geometries that are known
            are those associated with reaches that are not collocated. On subsequent iterations, reach geometries
            in collocated cells are added.

        Returns
        -------
        dataframe containing node number, segment, reach, and geometry for each SFR reach.
        """

        print "Adding segment and reach information to linework shapefile..."
        print "(may take quite a while for large models...)"
        df = self.m1.copy()
        df_lines = GISio.shp2df(lines_shapefile)
        prj = lines_shapefile[:-4] + '.prj'
        if new_lines_shapefile is None:
            new_lines_shapefile = lines_shapefile[:-4] + '_seg_reach.shp'

        df_lines_geoms = df_lines.geometry.values
        geom_idx = [df_lines.index[df_lines.node == n] for n in df.node.values]
        df['geometry'] = [df_lines_geoms[g[0]] if len(g) == 1 else LineString() for g in geom_idx]

        node_col = 'node'
        # then assign geometries for model cells with multiple SFR reaches
        # use length to figure out which geometry goes with which reach
        shared_cells = np.unique(df.ix[df.node.duplicated(), self.node_attribute])
        non_collocated_geometries = df.geometry.tolist()
        ni = iterations * len(shared_cells)
        for i in range(iterations):
            for n in shared_cells:
                print '\r{:.0f}%%'.format(100* (i*n + n/ni)),
                # make dataframe of all reaches within the model cell
                reaches = pd.DataFrame(df_lines.ix[df_lines[node_col] == n, 'geometry'].copy())
                reaches['length'] = [g.length for g in reaches.geometry]
                #reaches.sort('length', inplace=True)

                # streamflow results for that model cell
                dfs = df.ix[df[self.node_attribute] == n]
                #dfs.sort('length', inplace=True)

                # for each collocated reach in reaches dataframe,
                # this inner loop may be somewhat inefficient,
                # but the number of collocated reaches is small enough for it not to matter
                assigned = []
                geoms_entirely_within_node = []
                for rg in reaches.geometry:

                    # get list of neighboring reach geometries that touch the current collocated geometry
                    intersects_ind = df[np.array([rg.intersects(g) for g in df.geometry]) &
                                        np.array([nd != n for nd in df.node])].index

                    # handle collocated geometries that do not touch other segments
                    #if len(intersects_ind) == 0:
                    #    geoms_entirely_within_node.append(rg)
                    #    continue

                    intersects_df = df.ix[intersects_ind]
                    intersects_segments = intersects_df.segment.tolist()
                    dfs_segment = dfs.ix[dfs.segment.isin(intersects_segments)]

                    # handle collocated geometries that do not touch other segments
                    if len(intersects_ind) == 0 or len(dfs_segment) == 0:
                        geoms_entirely_within_node.append(rg)
                        continue

                    # get the index of the collocated reach that is closest to the intersecting reach number(s)
                    # (handles cases where a segment meanders out of and then back into a cell)
                    # only consider reaches that have the same segment as those intersected
                    closest_reach_in_dfs_index = dfs_segment.index[np.argmin([np.sum(np.abs(r - intersects_df.reach))
                                                                   for r in dfs_segment.reach])]

                    # assign the the current collocated geometry to the collocated reach
                    df.loc[closest_reach_in_dfs_index, 'geometry'] = rg
                    assigned.append(closest_reach_in_dfs_index)

                # in case there was a geometry entirely contained within the cell
                # (hopefully there is only one)
                if len(geoms_entirely_within_node) == 1:
                    ind = dfs.index[~dfs.index.isin(assigned)][0]
                    df.loc[ind, 'geometry'] = geoms_entirely_within_node[0]

        self.linework_geoms = df[['segment', 'reach', 'node', 'geometry']].sort(['segment', 'reach'])
        GISio.df2shp(self.linework_geoms, new_lines_shapefile, prj=prj)

    def write_shapefile(self, outshp='SFR.shp', xll=None, yll=None, epsg=None, proj4=None, prj=None):

        for kwd in [xll, yll, epsg, proj4, prj]:
            if kwd is not None:
                self.__dict__[kwd] = kwd

        if 'geometry' not in self.m1.columns:
            try:
                self.get_cell_geometries()
            except:
                print 'No cell geometries found. Please add discretization information using read_dis2(), ' \
                'and then compute cell geometries by running get_cell_geometries().'
        GISio.df2shp(self.m1, shpname=outshp, epsg=self.epsg, proj4=self.proj4, prj=self.prj)

    def write_streamflow_shapefile(self, streamflow_file=None, lines_shapefile=None, node_col=None):
        """Write out SFR Package results from a MODFLOW run.

        Parameters
        ----------
        streamflow_file: string
            Text file (often <model name>_streamflow.dat) containing table of SFR results from a MODFLOW run.

        lines_shapefile: string
            Path to shapefile containing linework geometries for each SFR reach

        node_col: string
            Name of attribute field in lines_shapefile containing the model cell numbers for each SFR reach.

        """

        if streamflow_file is not None:
            self.streamflow_file = streamflow_file

        self.Streamflow = Streamflow(sfrobject=self)
        self.Streamflow.read_streamflow_file(streamflow_file=streamflow_file)
        self.Streamflow.write_streamflow_shp(lines_shapefile=lines_shapefile, node_col=node_col)

    def write_tables(self, basename='SFR'):
        m1, m2 = self.m1.copy(), self.m2.copy()
        for col in ['geometry', 'centroids']:
            if col in m1.columns:
                m1.drop(col, axis=1, inplace=True)

        for col in ['upsegs']:
            if col in m2.columns:
                m2.drop(col, axis=1, inplace=True)

        m1.to_csv('{}mat1.csv'.format(basename), index=False)
        m2.to_csv('{}mat2.csv'.format(basename), index=False)
        print 'Mat1 and 2 saved to {fname}mat1.csv and {fname}mat2.csv'.format(fname=basename)

    def write_sfr_package(self, basename='SFR', tpl=False,
                          minimum_slope=1e-4,
                          bedKmin=1e-6,
                          nsfrpar=0,
                          nparseg=0,
                          const=128390.4,
                          dleak=0.001,
                          istcb1=50,
                          istcb2=66,
                          isfropt=1,
                          nstrail=10,
                          isuzn=1,
                          nsfrsets=30,
                          global_stream_depth=1,
                          iface=None):
        """
        Method to write an SFR package file from the Mat 1 and 2 (m1 and m2) attributes

        Parameters
        ----------

        iface: (int)
            An optional keyword that indicates that an IFACE value will be read for each reach and
            written to the budget file so that MODPATH can track particles to an exit face (or from an entry
            face) of a stream cell.  If "IFACE" is found at the end of record 1c, then IFACE values are read
            from the end of each item 2 record.
        """

        m1, m2 = self.m1.copy(), self.m2.copy()

        if tpl:
            ext = '.tpl'
        else:
            ext = '.sfr'
        outfile = basename + ext

        print 'writing {}'.format(outfile)
        ofp = open(outfile, 'w')

        nreaches = len(m1)
        nseg = len(m2)

        if tpl:
            ofp.write("ptf ~\n")
        ofp.write("#SFRpackage file generated by SFRmaker\n")
        record_1c = '{0:d} {1:d} {2:d} {3:d} {4:e} {5:e} {6:d} {7:d} {8:d} {9:d} {10:d} {11:d}'.format(
            -1*nreaches,
            nseg,
            nsfrpar,
            nparseg,
            const,
            dleak,
            istcb1,
            istcb2,
            isfropt,
            nstrail,
            isuzn,
            nsfrsets
        )
        if iface is not None:
            record_1c += ' IFACE'


        ofp.write(record_1c + '\n')

        for i, r in m1.iterrows():

            slope = r.bed_slope if r.bed_slope > minimum_slope else minimum_slope
            bedK = '~SFRc~' if tpl and r.bed_K > bedKmin else '{:e}'.format(r.bed_K)

            item2_record = '{0:.0f} {1:.0f} {2:.0f} {3:.0f} {4:.0f} {5:e} {6:e} {7:e} {8:e} {9:s}\n'.format(
                r.layer,
                r.row,
                r.column,
                r.segment,
                r.reach,
                r.length_in_cell,
                r.top_streambed,
                slope,
                r.bed_thickness,
                bedK)

            if iface is not None:
                item2_record += ' {:.0f}'.format(iface)

            ofp.write(item2_record + '\n')

        ofp.write('{0:.0f} 0 0 0\n'.format(nseg))
        for i, r in m2.iterrows():
            ofp.write('{0:.0f} {1:.0f} {2:.0f} 0 {3:e} 0.0 0.0 0.0 {4:e}\n'.format(
                r.segment,
                r.icalc,
                r.outseg,
                r.flow,
                r.roughch
                ))

            width = m1.ix[m1.segment == r.segment, 'width_in_cell'].values
            if width[0] > 0 and r.icalc == 1:
                ofp.write('{0:e}\n'.format(width[0]))
                ofp.write('{0:e}\n'.format(width[-1]))
            elif r.icalc == 0:
                ofp.write('{0:e} {1:e}\n'.format(width[0], global_stream_depth))
                ofp.write('{0:e} {1:e}\n'.format(width[-1], global_stream_depth))
            else:
                print 'icalc values >1 not supported.'
                return

        ofp.close()
        print 'Done'


class Elevations(SFRdata):

    def __init__(self, sfrobject=None, Mat1=None, Mat2=None, sfr=None, node_column=False,
                 mfpath=None, mfnam=None, mfdis=None, to_meters_mult=0.3048,
                 minimum_slope=1e-4, landsurfacefile=None, landsurface_column=None, smoothing_iterations=0):
        """
        Smooth streambed elevations outside of the context of the objects in SFR classes
        (works off of information in Mat1 and Mat2; generates updated versions of these files
        """
        SFRdata.__init__(self, sfrobject=sfrobject, Mat1=Mat1, Mat2=Mat2, sfr=sfr, node_column=node_column,
                 mfpath=mfpath, mfnam=mfnam, mfdis=mfdis, landsurfacefile=landsurfacefile, to_meters_mult=to_meters_mult,
                 minimum_slope=minimum_slope)

        if self.landsurfacefile and self.landsurface is None:
            self.landsurface = np.fromfile(self.landsurfacefile, sep=' ') # array of elevations to use (sorted by cellnumber)

            print 'assigning elevations in {} to landsurface column in Mat1...'.format(landsurfacefile)
            # assign land surface elevations based on node number
            self.m1['landsurface'] = [self.landsurface[n-1] for n in self.m1[self.node_column]]
            self.m1['top_streambed'] = self.m1['landsurface'] # might consider getting rid of landsurface column and only working with top_streambed

        elif not landsurfacefile and landsurface_column is not None:
            print 'assigning elevations in Mat1 {} column to landsurface column...'.format(landsurface_column)
            self.m1['landsurface'] = self.m1[landsurface_column]

        self.smoothing_iterations = smoothing_iterations

    def smooth_segment_ends(self, report_file='smooth_segment_ends.txt'):
        '''
        smooth segment end elevations so that they decrease monotonically down the stream network
        '''
        '''
        self.m2 = self.m2.replace(999999, 0)
        print '\nSmoothing segment ends...\n'
        # open a file to report a summary of the elevation adjustments
        self.ofp = open(report_file, 'w')
        self.ofp.write('Segment end smoothing report:\n'
                  'segment,max_elev,min_elev,downstream_min_elev\n')

        # set initial max / min elevations for segments based on max and min land surface elevations along each segment
        self.seg_maxmin = np.array([(np.max(self.m1.landsurface[self.m1.segment == s]),
                      np.min(self.m1.landsurface[self.m1.segment == s])) for s in self.segments])

        # determine constraining upstream min elevations for each segment in Mat2
        upstreamMin = np.array([np.min([self.seg_maxmin[useg-1][1] for useg in self.m2.ix[s, 'upsegs']])
                                  if len(self.m2.ix[s, 'upsegs']) > 0
                                  else self.seg_maxmin[s-1][0] for s in self.m2.index])

        # if the upstream minimum elevation is lower than the max elevation in the segment,
        # reset the max to the upstream minimum
        self.seg_maxmin = np.array([(upstreamMin[s-1], self.seg_maxmin[s-1][1])
                                    if upstreamMin[s-1] < self.seg_maxmin[s-1][0]
                                    else self.seg_maxmin[s-1] for s in self.segments])

        # check for higher minimum elevations in each segment
        bw = dict([(i, maxmin) for i, maxmin in enumerate(self.seg_maxmin) if maxmin[0] < maxmin[1]])

        # smooth any backwards elevations
        if len(bw) > 0:

            # creates dataframe of all outsegs for each segment
            self.map_outsegs()

            # iterate through the mapped outsegs, replacing minimum elevations until there are no violations
            #self.fix_backwards_ends(bw)

        print 'segment ends smoothing finished in {} iterations.\n' \
              'segment ends saved to {}\n' \
              'See {} for report.'\
            .format(self.smoothing_iterations, self.Mat2[:-4] + '_elevs.csv', report_file)

        # populate Mat2 dataframe with bounding elevations for upstream and downstream segments, so they can be checked
        self.m2['upstreamMin'] = [np.min([self.seg_maxmin[useg-1][1] for useg in self.m2.ix[s, 'upsegs']])
                                  if len(self.m2.ix[s, 'upsegs']) > 0
                                  else self.seg_maxmin[s-1][0] for s in self.m2.index]
        self.m2['upstreamMax'] = [np.min([self.seg_maxmin[useg-1][0] for useg in self.m2.ix[s, 'upsegs']])
                                  if len(self.m2.ix[s, 'upsegs']) > 0
                                  else self.seg_maxmin[s-1][0] for s in self.m2.index]
        self.m2['Max'] = self.seg_maxmin[:, 0]
        self.m2['Min'] = self.seg_maxmin[:, 1]
        self.m2['downstreamMax'] = [self.seg_maxmin[s-1][0] for s in self.m2.outseg]
        self.m2['downstreamMin'] = [self.seg_maxmin[s-1][1] for s in self.m2.outseg]

        elevations_summary = self.m2[['segment','upstreamMax','upstreamMin','Max','Min','downstreamMax','downstreamMin']]
        elevations_summary.to_csv(self.ofp)
        self.ofp.close()

        # save new copy of Mat2 with segment end elevations
        self.m2 = self.m2.drop(['upstreamMax', 'upstreamMin', 'downstreamMax', 'downstreamMin', 'upsegs'], axis=1)
        self.m2.to_csv(self.Mat2[:-4] + '_elevs.csv', index=False)
        '''

    def replace_downstream(self, bw, level, ind):
        '''
        replace minimum elevations in segments with either the max or min elevation in downstream segment
        bw = dict with keys that are indices of segements to modify
        level = column of outsegs to reference in outsegs table
        ind = 0 (replace with downstream max elev) or 1 (min elev)
        get downstream max elevations for level from outsegs at that level (for segments with backwards elevs)
        '''
        # make list of downstream elevations (max or min, depending on ind)
        # if segment is an outlet, use minimum elevation of segment
        downstream_elevs = [self.seg_maxmin[self.outsegs.ix[s, level] - 1][ind]
                            if self.outsegs.ix[s, level] > 0
                            else np.min(self.seg_maxmin[s-1]) for s in self.segments]

        # assign any violating minimum elevations to downstream elevs from above
        bw_inds = bw.keys()
        #print 'segment {}: {}, downseg: {}'.format(bw_inds[0]+1, self.seg_maxmin[bw_inds[0]], downstream_elevs[bw_inds[0]])
        self.seg_maxmin[bw_inds] = np.array([(self.seg_maxmin[i][0], downstream_elevs[i]) for i in bw_inds])
        for i in bw_inds:
            self.ofp.write('{},{},{}\n'.format(i+1, self.seg_maxmin[i], downstream_elevs[i]))


    def replace_upstream(self):
        '''
        replace maximum elevations in segments with minimum elevations in upstream segments,
        in segments where the max is above the upstream minimum
        '''
        # make list of min upstream elevations
        # if segment is a headwater, use maximum elevation of segment
        upstreamMin = np.array([np.min([self.seg_maxmin[useg-1] for useg in self.m2.ix[s, 'upsegs']])
                       if len(self.m2.ix[s, 'upsegs']) > 0
                       else self.seg_maxmin[s-1][0] for s in self.segments])

        # if the upstream minimum elevation is lower than the max elevation in the segment,
        # reset the max to the upstream minimum
        self.seg_maxmin = np.array([(upstreamMin[s-1], self.seg_maxmin[s-1][1])
                                    if upstreamMin[s-1] < self.seg_maxmin[s-1][0]
                                    else self.seg_maxmin[s-1] for s in self.segments])


    def fix_backwards_ends(self, bw):

        knt = 1
        for level in self.outsegs.columns:

            # replace minimum elevations in backwards segments with downstream maximums
            self.replace_downstream(bw, level, 0)

            # check for higher minimum elevations in each segment
            bw = dict([(i, maxmin) for i, maxmin in enumerate(self.seg_maxmin) if maxmin[0] < maxmin[1]])

            # if still backwards elevations, try using downstream min elevations
            if len(bw) > 0:

                # replace minimum elevations in backwards segments with downstream minimums
                self.replace_downstream(bw, level, 1)

            # in segments where the max elevation is higher than the upstream minimum,
            # replace max elev with upstream minimum
            self.replace_upstream()

            # check again for higher minimum elevations in each segment
            bw = dict([(i, maxmin) for i, maxmin in enumerate(self.seg_maxmin) if maxmin[0] < maxmin[1]])

            # stop if there are no longer any backwards segments, otherwise continue
            if len(bw) == 0:
                break

            knt += 1

        self.smoothing_iterations = knt


    def interpolate(self, istart, istop):
        '''
        istart = index of previous minimum elevation
        istop = index of current minimum elevation
        dx = interpolation distance
        dS = change in streambed top over interpolation distance
        '''
        dx = self.cdist[istart] - self.cdist[istop]
        dS = self.minelev - self.segelevs[istop]

        dist = 0

        if dS == 0:
            slope = 0
        else:
            slope = dS/dx

        for i in np.arange(istart+1, istop+1):

            dist += self.cdist[i-1] - self.cdist[i]
            self.sm.append(self.minelev - dist * slope)

            self.ofp.write('{:.0f},{:.0f},{:.2f},{:.2f},{:.2f},{:.2e},{:.2f}\n'
                           .format(self.seg, i, self.segelevs[i-1], self.minelev, dist, slope, self.sm[-1]))

        # reset the minimum elevation to current
        self.minelev = self.sm[-1]


    def smooth_segment_interiors(self, report_file='smooth_segment_interiors.txt'):

        print '\nSmoothing segment interiors...\n'

        # open a file to report on interpolations
        self.ofp = open(report_file, 'w')
        self.ofp.write('segment, reach, land_surface, minelev, dist, slope, sb_elev\n')

        try:
            self.m2['Max'], self.m2['Min']
            m2minmax = True
        except:
            m2minmax = False
            print "Max, Min elevation columns not found in Mat2" \
                  "Run smooth_segment_ends() first."
            return

        for seg in self.segments:
            self.seg = seg

            # make a view of the Mat 1 dataframe that only includes the segment
            df = self.m1[self.m1.segment == seg].sort('reach')

            # start with land surface elevations along segment
            self.segelevs = df.landsurface.values

            start, end = self.m2.ix[seg, ['Max', 'Min']]

            # get start and end elevations from Mat2; if they are equal, continue
            if start == end:
                self.sm = [start] * len(self.segelevs)
                self.m1.loc[self.m1.segment == seg, 'top_streambed'] = self.sm
                continue
            else:
                self.sm = [start]
                self.segelevs[-1] = end

            # calculate cumulative distances at cell centers
            lengths = df.length_in_cell.values
            self.cdist = np.cumsum(lengths) - 0.5 * lengths

            self.minelev = start
            minloc = 0
            nreaches = len(self.segelevs)

            for i in range(nreaches)[1:]:

                # if the current elevation is equal to or below the minimum
                if self.segelevs[i] <= self.minelev:

                    # if the current elevation is above the end, interpolate from previous minimum to current
                    # only interpolate to current reach if it is below the previous minimum
                    # (that is, do not allow flat stretches unless the segment minimum elevation is reached)
                    if self.minelev > self.segelevs[i] >= end:
                        self.interpolate(minloc, i)
                        minloc = i

                    # otherwise if it is below the end, interpolate from previous minimum to end
                    elif self.segelevs[i] <= end:
                        self.interpolate(minloc, nreaches-1)
                        break
                else:
                    continue

            # update Mat1 with smoothed streambed tops
            self.m1.loc[self.m1.segment == seg, 'top_streambed'] = self.sm

        # assign slopes to Mat 1 based on the smoothed elevations
        self.calculate_slopes()

        # save updated Mat1
        self.m1.to_csv(self.Mat1[:-4] + '_elevs.csv', index=False)
        self.ofp.close()
        print 'Done, updated Mat1 saved to {}.\nsee {} for report.'.format(self.Mat1[:-4] + '_elevs.csv', report_file)


    def calculate_slopes(self):
        '''
        assign a slope value for each stream cell based on streambed elevations
        this method is run at the end of smooth_segment_interiors, after the interior elevations have been assigned
        '''
        print 'calculating slopes...'
        for s in self.segments:

            # calculate the right-hand elevation differences (not perfect, but easy)
            diffs = self.m1[self.m1.segment == s].top_streambed.diff()[1:].tolist()

            if len(diffs) > 0:
                diffs.append(diffs[-1]) # use the left-hand difference for the last reach

            # edge case where segment only has 1 reach
            else:
                diffs = self.m2.Min[s] - self.m2.Max[s]

            # divide by length in cell; reverse sign so downstream = positive (as in SFR package)
            slopes = diffs / self.m1[self.m1.segment == s].length_in_cell * -1

            # assign to Mat1
            self.m1.loc[self.m1.segment == s, 'bed_slope'] = slopes
        
        # enforce minimum slope
        ns = [s if s > self.minimum_slope else self.minimum_slope for s in self.m1.bed_slope.tolist()]
        self.m1['bed_slope'] = ns


    def reset_model_top_2streambed(self, minimum_thickness=1, outdisfile=None, outsummary=None):
        """Make the model top elevation consistent with the SFR streambed elevations;
        Adjust other layers downward (this puts all SFR cells in layer1)

        Parameters
        ----------
        minimum_thickness : float
            Minimum layer thickness to enforce when adjusting underlying layers to accommodate
            changes to the model top

        outdisfile : string
            Name for new MODFLOW discretization file

        outsummary : string (optional)
            If specified, will save a summary (in SFR Mat1 style format)
            of adjustments made to the model top

        """
        if outdisfile is None:
            outdisfile = self.mfdis[:-4] + '_adjusted_to_streambed.dis'

        # output file summarizing adjustments made to model top
        if outsummary is None:
            outsummary = self.mfdis[:-4] + '_adjustments_to_model_top.csv'

        # make sure that streambed thickness is less than minimum_thickness,
        # otherwise there is potential for MODFLOW altitude errors
        # (with streambed top == model top, the streambed bottom will be below the layer 1 bottom)
        if np.min(self.m1.bed_thickness >= minimum_thickness):
            self.m1.bed_thickness = 0.9 * minimum_thickness

        # make a vector of highest streambed values for each cell containing SFR (for collocated SFR cells)
        self.m1['lowest_top'] = self.m1.top_streambed

        #shared_cells = np.unique(self.m1.ix[self.m1.node.duplicated(), 'node'])
        for c in self.shared_cells:

            # select the collocated reaches for this cell
            df = self.m1[self.m1.node == c].sort('top_streambed', ascending=False)

            # make column of lowest streambed elevation in each cell with SFR
            self.m1.loc[df.index, 'lowest_top'] = np.min(df.top_streambed)

        # make a new model top array; assign highest streambed tops to it
        newtop = self.dis.top.array.copy()
        newtop[self.m1.row.values-1, self.m1.column.values-1] = self.m1.lowest_top.values

        # Now straighten out the other layers, removing any negative thicknesses
        # do layer 1 first
        newbots = self.dis.botm.array.copy()
        conflicts = newbots[0, :, :] > newtop - minimum_thickness
        newbots[0, conflicts] = newtop[conflicts] - minimum_thickness

        for i in range(self.dis.nlay - 1):
            conflicts = newbots[i+1, :, :] > (newbots[i, :, :] - 1)
            newbots[i+1, conflicts] = newbots[i, conflicts] - minimum_thickness

        # make a dataframe that shows the largest adjustments made to model top
        self.m1['top_height'] = self.m1.model_top - self.m1.top_streambed
        adjustments = self.m1[self.m1.top_height > 0].sort(columns='top_height', ascending=False)
        adjustments.to_csv(outsummary)

        # update the model top in Mat1
        self.m1['model_top'] = self.m1.lowest_top

        # update the elevs array
        self.elevs[0, :, :] = newtop
        self.elevs[1:, :, :] = newbots

        # update the layer in Mat1 to 1 for all SFR cells
        self.m1['layer'] = 1

        # write a new DIS file
        new_m = flopy.modflow.mf.Modflow(model_ws=os.path.split(outdisfile)[0],
                                         modelname=os.path.split(outdisfile)[1][:-4])
        newdis = flopy.modflow.ModflowDis(new_m, nlay=self.dis.nlay, nrow=self.dis.nrow, ncol=self.dis.ncol,
                                          delr=self.dis.delr, delc=self.dis.delc, top=newtop, botm=newbots)

        self.mfdis = outdisfile
        newdis.write_file()


    def incorporate_field_elevations(self, shpfile, elevs_field, distance_tol):
        """Update landsurface elevations for SFR cells from nearby field measurements
        """
        # read field measurements into a dataframe
        df = GISio.shp2df(shpfile)

        self.get_cell_centroids()

        # make vectors of all x and y of SFR cell centroids
        X, Y = np.array([c[0] for c in self.m1.centroids]), np.array([c[1] for c in self.m1.centroids])

        # determine the closest SFR reach to each field measurement, within distance_tol
        # if there are no SFR measurements within distance_tol,
        closest_SFRmat1_inds = [np.argmin(np.sqrt((g.x - X)**2 + (g.y - Y)**2)) for g in df.geometry]
        closest_SFR_distances = [np.min(np.sqrt((g.x - X)**2 + (g.y - Y)**2)) for g in df.geometry]

        # update the landsurface column in mat1 with the field measurements within distance_tol
        # the landsurface column is used in interpolating streambed elevations
        for i, mat1_ind in enumerate(closest_SFRmat1_inds):
            print "closest Mat1 ind: {}".format(mat1_ind)
            print "distance: {}".format(closest_SFR_distances[i])
            if closest_SFR_distances[i] <= distance_tol:
                self.m1.loc[mat1_ind, 'landsurface'] = df.iloc[i][elevs_field]

                # reset the minimum elevation for the segment if the field elevation is lower
                segment = self.m1.loc[mat1_ind, 'segment']
                if 'Min' in self.m2.columns and df.iloc[i][elevs_field] < self.m2.ix[segment, 'Min']:
                    self.m2.loc[segment, 'Min'] = df.iloc[i][elevs_field]

    def map_confluences(self):

        # if min/max elevations not in Mat2, compute from Mat1
        if 'Max' not in self.m2.columns:
            print 'Setting initial segment min/max elevations from Mat1...'
            top_streambed = self.m1.top_streambed.values
            m1segments = self.m1.segment.values
            m2segments = self.m2.segment.values
            self.m2['Max'] = [np.max(top_streambed[m1segments == s]) for s in m2segments]
            self.m2['Min'] = [np.min(top_streambed[m1segments == s]) for s in m2segments]
        m2upsegs = self.m2.upsegs.tolist()

        # setup dataframe of confluences
        # confluences are where segments have upsegs (no upsegs means the reach 1 is a headwater)
        confluences = self.m2.ix[np.array([len(u) for u in m2upsegs]) > 0, ['segment', 'upsegs']].copy()
        confluences['node'] = [0] * len(confluences)
        confluences['elev'] = [0] * len(confluences)
        nconfluences = len(confluences)
        print 'Mapping {} confluences and updating segment min/max elevations in Mat2...'.format(nconfluences)
        for i, r in confluences.iterrows():
            # get node/cellnum for downstream segment reach 1 (this is where the confluence is located)
            node = self.m1.ix[(self.m1.segment == i) & (self.m1.reach == 1), 'node']
            confluences.loc[i, 'node'] = node.values[0]

            # include land surface and top of streambed columns to ensure that minimums are found
            # (e.g. if land surface was updated and top_streambed wasn't)
            elevation_columns = ['top_streambed', 'landsurface']
            for c in elevation_columns:
                if c not in self.m1.columns:
                    elevation_columns.remove(c)
            if len(elevation_columns) == 0:
                raise IndexError("top_streambed and landsurface columns not found in Mat1!")

            # confluence elevation is the minimum of the ending segments minimums, starting segments maximums
            endsmin = np.min(self.m1.ix[self.m1.segment.isin(r.upsegs), elevation_columns].values)
            startmax = np.max(self.m1.ix[self.m1.segment == i, elevation_columns].values)
            cfelev = np.min([endsmin, startmax])
            confluences.loc[i, 'elev'] = cfelev

            # update Mat2
            self.m2.loc[r.upsegs, 'Min'] = cfelev
            self.m2.loc[i, 'Max'] = cfelev

        # check for Mins that are higher than maxes
        #diffs = np.ravel(np.diff(self.m2[['Min', 'Max']].values))

        # Go through again and make sure that Mins are all <= Maxes
        for i, r in confluences.iterrows():
            endsmin = np.min(self.m2.ix[self.m2.segment.isin(r.upsegs), ['Min', 'Max']].values)
            startmax = np.max(self.m2.ix[self.m2.segment==i, ['Min', 'Max']].values)
            cfelev = np.min([endsmin, startmax])
            confluences.loc[i, 'elev'] = cfelev

            self.m2.loc[r.upsegs, 'Min'] = cfelev
            self.m2.loc[i, 'Max'] = cfelev

        self.confluences = confluences
        print 'Done, see confluences attribute.'


class Widths(SFRdata):

    def __init__(self, sfrobject=None, Mat1=None, Mat2=None, sfr=None, to_km=0.0003048, Mat2_out=None):

        SFRdata.__init__(self, sfrobject=sfrobject, Mat1=Mat1, Mat2=Mat2, sfr=sfr, to_meters_mult=to_km*1000, Mat2_out=Mat2_out)

        self.to_km = to_km # multiplier from model units to km

    def widthcorrelation(self, arbolate):
        #estimate widths, equation from Feinstein and others (Lake
        #Michigan Basin model) width=0.1193*(x^0.5032)
        # x=arbolate sum of stream upstream of the COMID in meters
        #NHDPlus has arbolate sum in kilometers.
        #print a table with reachcode, order, estimated width, Fcode
        estwidth = 0.1193 * (1000 * arbolate) ** 0.5032
        return estwidth


    def map_upsegs(self):
        '''
        from Mat2, returns dataframe of all upstream segments (will not work with circular routing!)
        '''
        # make a list of adjacent upsegments keyed to outseg list in Mat2
        upsegs = dict([(o, self.m2.segment[self.m2.outseg == o].tolist()) for o in np.unique(self.m2.outseg)])
        self.outsegs = [k for k in upsegs.keys() if k > 0] # exclude 0, which is the outlet designator

        # for each outseg key, for each upseg, check for more upsegs, append until headwaters has been reached
        for outseg in self.outsegs:

            up = True
            upsegslist = upsegs[outseg]
            while up:
                added_upsegs = []
                for us in upsegslist:
                    if us in self.outsegs:
                        added_upsegs += upsegs[us]
                if len(added_upsegs) == 0:
                    up = False
                    break
                else:
                    upsegslist = added_upsegs
                    upsegs[outseg] += added_upsegs

        # the above algorithm is recursive, so lower order streams get duplicated many times
        # use a set to get unique upsegs
        self.upsegs = dict([(u, list(set(upsegs[u]))) for u in self.outsegs])


    def estimate_from_arbolate(self):

        print 'estimating stream widths...'

        self.m2.in_arbolate = self.m2.in_arbolate.fillna(0) # replace any nan values with zeros

        # map upsegments for each outlet segment
        self.map_upsegs()

        # compute starting arbolate sum values for all segments (sum lengths of all mapped upsegs)
        asum_start = {}
        for oseg in self.outsegs:
            asum_start[oseg] = 0
            for us in self.upsegs[oseg]:
                # add the sum of all lengths for the upsegment
                asum_start[oseg] += self.m1.ix[self.m1.segment == us, 'length_in_cell'].sum() * self.to_km

                # add on any starting arbolate sum values from outside the model (will add zero if nothing was entered)
                asum_start[oseg] += self.m2.in_arbolate[us]
            #print 'outseg: {}, starting arbolate sum: {}'.format(oseg, asum_start[oseg])

        # assign the starting arbolate sum values to Mat2
        asum_starts = [asum_start[i] if i in self.outsegs
                       else self.m2.ix[i, 'in_arbolate'] for i, r in self.m2.iterrows()]
        self.m2['starting_arbolate'] = asum_starts

        for s in self.segments:

            segdata = self.m1[self.m1.segment == s] # shouldn't need to sort, as that was done in __init__

            # compute arbolate sum at each reach, in km, including starting values from upstream segments
            asums = segdata.length_in_cell.cumsum() * self.to_km + self.m2.ix[s, 'starting_arbolate']

            # compute width, assign to Mat1
            self.m1.loc[self.m1.segment == s, 'width_in_cell'] = [self.widthcorrelation(asum) for asum in asums]

        self.m1.to_csv(self.Mat1_out, index=False)
        print 'Done, updated Mat1 saved to {}.'.format(self.Mat1_out)

        if self.Mat2_out:
            self.m2.to_csv(self.Mat2_out, index=False)
            print 'saved arbolate sum information to {}.'.format(self.Mat2_out)


class Outsegs(SFRdata):

    def plot_outsegs(self, outpdf='complete_profiles.pdf', units='ft', add_profiles={}):
        print "listing unique segment sequences..."
        # make a dataframe of all outsegs
        self.map_outsegs()

        # list of lists of connected segments
        seglists = [[p]+[s for s in self.outsegs.ix[p, :] if s > 0] for p in self.outsegs.index]

        # reduce seglists down to unique sets of segments
        unique_lists = []
        for l in seglists:
            # find all lists for which the current list is a subset
            overlap = np.where([set(l).issubset(set(ll)) for ll in seglists])
            # length 1 means that there are no other lists for which the current is a subset
            if len(overlap[0]) == 1:
                unique_lists.append(l)
            else:
                continue

        # add model top elevations to dictionary using node numbers
        if 'model_top' not in self.m1.columns:
            self.m1['model_top'] = [self.elevs_by_cellnum[c] for c in self.m1[self.node_column]]

        pdf = PdfPages(outpdf)
        for l in unique_lists:
            if len(l) > 1:
                cdist = []
                cdist_end = 0
                cdist_ends = [0]
                sbtop = []
                modeltop = []
                self.profile_elevs = {}

                for p in add_profiles.iterkeys():
                    self.profile_elevs[p] = []

                # make each profile from its segments
                for s in l:
                    # make a view of the Mat 1 dataframe that only includes the segment
                    df = self.m1[self.m1.segment == s].sort('reach')

                    # calculate cumulative distances at cell centers
                    cdist += (np.cumsum(df['length_in_cell']) - 0.5 * df['length_in_cell'] + cdist_end).tolist()

                    cdist_end += np.sum(df['length_in_cell']) # add length of segment to total
                    #print cdist_end
                    cdist_ends.append(cdist_end) # record distances at segment boundaries

                    # streambed tops
                    sbtop += df['top_streambed'].tolist()

                    # model tops
                    modeltop += df['model_top'].tolist()

                    for p in add_profiles.iterkeys():
                        #elevs = self.cells2vertices(df[add_profiles[p]].tolist())
                        self.profile_elevs[p] += df[add_profiles[p]].tolist()

                # plot the profile
                fig = plt.figure()
                ax = fig.add_subplot(111)
                plt.plot(cdist, modeltop, label='Model top', lw=0.5)
                plt.plot(cdist, sbtop, label='Streambed top')

                # add any additional profiles
                for p in add_profiles.iterkeys():
                    #print len(cdist)
                    #print len(self.elevs[p])
                    plt.plot(cdist, self.profile_elevs[p], label=p, lw=0.5)

                # add in lines marking segment starts/ends
                # add in annotations labeling each segment
                for i, x in enumerate(cdist_ends[1:]):
                    plt.axvline(x, c='0.25', lw=0.5)

                    tx = 0.5 * (cdist_ends[i] + cdist_ends[i-1])
                    ty = ax.get_ylim()[0] + 0.9 * (ax.get_ylim()[1] - ax.get_ylim()[0])
                    ax.text(tx, ty, str(l[i-1]), ha='center', fontsize=8)

                ax.set_xlabel('Distance along SFR path')
                ax.set_ylabel('Elevation, {}'.format(units))
                title = 'SFR profile for segments:\n'
                for s in l:
                    title += ' {:.0f}'.format(s)
                ax.set_title(title, fontsize=10, fontweight='bold', loc='left')
                ax.set_ylim(int(np.floor(ax.get_ylim()[0])), int(np.ceil(ax.get_ylim()[1])))
                plt.legend(loc='best', fontsize=10)
                pdf.savefig()
                plt.close()
        pdf.close()

    def plot_routing(self, outpdf='routing.pdf'):

        # make a dataframe of all outsegs
        self.map_outsegs()

        # make a matrix of the model grid, assign outlet values to SFR cells
        self.watersheds = np.empty((self.dis.nrow * self.dis.ncol))
        self.watersheds[:] = np.nan

        for s in self.segments:

            segdata = self.m1[self.m1.segment == s]
            cns = segdata[self.node_column].values
            self.watersheds[cns] = segdata.Outlet.values[0]

        self.watersheds = np.reshape(self.watersheds, (self.dis.nrow, self.dis.ncol))

        # now plot it up!
        self.watersheds[self.watersheds == 999999] = np.nan # screen out 999999 values
        self.routingfig = plt.figure()
        ax = self.routingfig.add_subplot(111)
        plt.imshow(self.watersheds)
        cb = plt.colorbar()
        cb.set_label('Outlet segment')
        plt.savefig(outpdf, dpi=300)

        self.m1.to_csv(self.Mat1, index=False)
        print 'Done, Mat1 updated with outlet information; saved to {}.'.format(self.Mat1)


class Streamflow(SFRdata):

    def __init__(self, sfrobject=None, Mat1=None, Mat2=None, sfr=None, streamflow_file=None):

        SFRdata.__init__(self, sfrobject=sfrobject, Mat1=Mat1, Mat2=Mat2, sfr=sfr, streamflow_file=streamflow_file)

    def read_streamflow_file(self, streamflow_file=None):

        if streamflow_file is not None:
            self.streamflow_file = streamflow_file

        h = header(self.streamflow_file)
        ofp = open(self.streamflow_file)
        for i in np.arange(h):
            ofp.readline()

        nreaches = len(self.m1)
        sfrresults = {}
        for i in np.arange(nreaches):
            line = ofp.readline().strip().split()
            l, r, c, s, reach = map(int, line[0:5])
            Qin, Qgw, Qout, Qovr, Qp, Qet, S, d, w, Cond, sb_slope = map(float, line[5:])

            Qstream = 0.5 * (Qin + Qout)

            if Qgw > 0:
                state = 'losing'
            elif Qgw < 0:
                state = 'gaining'
            elif Qstream == 0:
                state= 'dry'

            sfrresults[i] = {'layer': l,
                             'row': r,
                             'column': c,
                             'segment': s,
                             'reach': reach,
                             'Qin': Qin,
                             'Qgw': Qgw,
                             'Qstream': Qstream,
                             'Qout': Qout,
                             'Qovr': Qovr,
                             'Qp': Qp,
                             'Qet': Qet,
                             'stage': S,
                             'depth': d,
                             'width': w,
                             'cond': Cond,
                             'sb_slope': sb_slope,
                             'state': state}

        return pd.DataFrame.from_dict(sfrresults, orient='index')

    def write_streamflow_shp(self, lines_shapefile=None, node_col=None):

        if node_col is None:
            node_col = self.node_column

        streamflow_shp = self.streamflow_file[:-4] + '.shp'

        # read in the streamflow results from the MODFLOW run
        df = self.read_streamflow_file()

        # join in node column from mat1
        # segments and reaches in SFR results and Mat1 must be identical! (and sorted)
        if np.max(self.m1.reach - df.reach) != 0 or np.max(self.m1.segment - df.segment) != 0:
            raise IndexError('Mismatch in segment and reach ordering between Mat1 and {}!\
            \nCheck that segments and reaches in Mat1 are identical to those in the SFR output.')\
                .format(self.streamflow_file)

        df['node'] = self.m1[self.node_attribute]
        df['length'] = self.m1.length_in_cell

        # if a shapefile is provided, get the geometries from there (by node)
        if lines_shapefile is not None:
            df_lines = GISio.shp2df(lines_shapefile)
            prj = lines_shapefile[:-4] + '.prj'

            if df_lines[node_col].dtype.name == 'object':
                raise TypeError('Node number column of shapefile is not numeric!')

            if 'segment' not in df_lines.columns:
                print "segment and reach information not in linework shapefile!"
                df_lines = self.segment_reach2linework_shapefile(lines_shapefile)

            df_lines.sort(['segment', 'reach'], inplace=True)

            # join in node column from mat1
            # segments and reaches in SFR results and Mat1 must be identical! (and sorted)
            if np.max(df_lines.reach - df.reach) != 0 or np.max(df_lines.segment - df.segment) != 0:
                raise IndexError('Segments and reaches in {} are different than {}!'\
                                 .format(lines_shapefile, self.streamflow_file))

            df['geometry'] = df_lines.geometry
            '''
            # first assign geometries for model cells with only 1 SFR reach
            df_lines_geoms = df_lines.geometry.values
            geom_idx = [df_lines.index[df_lines.node == n] for n in df.node.values]
            df['geometry'] = [df_lines_geoms[g[0]] if len(g) == 1 else LineString() for g in geom_idx]

            # then assign geometries for model cells with multiple SFR reaches
            # use length to figure out which geometry goes with which reach
            shared_cells = np.unique(df.ix[df.node.duplicated(), self.node_attribute])
            for n in shared_cells:
                # make dataframe of all reaches within the model cell
                reaches = pd.DataFrame(df_lines.ix[df_lines[node_col] == n, 'geometry'].copy())
                reaches['length'] = [g.length for g in reaches.geometry]

                # streamflow results for that model cell
                dfs = df.ix[df[self.node_attribute] == n]
                # this inner loop may be somewhat inefficient,
                # but the number of collocated reaches is small enough for it not to matter
                for i, r in reaches.iterrows():
                    # index of streamflow results with length closest to reach
                    ind = np.argmin(np.abs(dfs.length - r.length))
                    # assign the reach geometry to the streamflow results at that index
                    df.loc[ind, 'geometry'] = r.geometry
            '''
        # otherwise get the geometries from the dis file and model origin
        # (right now this does not support rotated grids)
        # the geometries are for the model cells- collocated reaches will be represented by
        # one model cell polygon for each reach
        else:
            self.get_cell_geometries()
            df_lines = pd.DataFrame.from_dict({'geometry': self.cell_geometries}, orient='rows')
            prj = self.prj
            df['geometry'] = [df_lines.geometry[n] for n in df.node.values]

        GISio.df2shp(df, streamflow_shp, prj=prj)


class Segments(SFRdata):

    def __init__(self, sfrobject=None, Mat1=None, Mat2=None, sfr=None):

        SFRdata.__init__(self, sfrobject=sfrobject, Mat1=Mat1, Mat2=Mat2, sfr=sfr)

    def index_downstream_reaches(self):
        """Get index of downstream reach for each SFR reach in Mat
        """
        print "indexing downstream reaches..."
        '''
        # list outsegs for each reach in m1
        outsegs = [int(self.m2.outseg[self.m2.segment == s]) for s in self.m1.segment]

        # list index of outseg reach 1 for each reach in m1
        outseg_reach1_inds = [self.m1[(self.m1.segment == o) & (self.m1.reach == 1)].index[0]
                       if o != 0 and o != 999999
                       else -999999 for o in outsegs]

        # downstream reach is next lowest index if segment is same, otherwise reach 1 of outseg
        # (m1 is sorted by segment and then reach)
        self.m1['downreach_ind'] = [i-1 if i > 1 and self.m1.segment[i-1] == self.m1.segment[i]
                                    else outseg_reach1_inds[i] for i, r in enumerate(self.m1.reach)]
        '''
        m1index = self.m1.index.values
        m2segments = self.m2.segment.values
        m2outsegs = self.m2.outseg.values
        segments = self.m1.segment.values
        reaches = self.m1.reach.values
        downreach_inds = []
        for i, seg in enumerate(segments):

            # if downstream reach is in current seg and not on reach index 0
            if i > 1 and segments[i-1] == seg:
                downreach_inds.append(i - 1)
            # otherwise, find the index of the first reach of the outseg
            else:
                outseg = int(m2outsegs[m2segments == seg])
                if outseg != 0 and outseg != 999999:
                    outseg_reach1_ind = m1index[(segments == outseg) & (reaches == 1)][0]
                else:
                    outseg_reach1_ind = -999999

                downreach_inds.append(outseg_reach1_ind)

        self.m1['downreach_ind'] = downreach_inds

    def renumber_SFR_cells(self, sfr_cells_list):
        """Renumbers the segments and reaches for a list of SFR cells, or a list containing lists of SFR cells.
        Each list of SFR cells must be contiguous.

        Parameters
        ----------
        sfr_cells_list : list
            List of model cells (by cell or node number) containing SFR reaches to rename. Can be a list of lists.

        Returns
        -------
        Updates segment and reach information in Mat1 and Mat2 dataframes
        """

        if not isinstance(sfr_cells_list[0], list):
            sfr_cells_list = [sfr_cells_list]

        def renumber(seg, old_reaches):

            # correct any offset in first reach
            reach1 = np.min(old_reaches)
            old_reaches -= reach1 - 1

            new_reaches = old_reaches
            new_segments = np.array([seg])

            reach1 = 1
            gaps = np.diff(old_reaches)
            for i, diff in enumerate(gaps):
                if diff > 1:
                    reach1 = old_reaches[i+1]
                    seg +=1

                new_reaches[i+1] -= reach1 - 1
                new_segments = np.append(new_segments, seg)

            return new_segments, new_reaches

        self.index_downstream_reaches()

        # for each intersect feature, create new SFR segments and renumber their reaches
        self.m1['old_segment'] = self.m1.segment
        self.m1['old_reach'] = self.m1.reach

        # create local variables to speed up loops
        m1segments = self.m1.segment.values
        m1reaches = self.m1.reach.values
        m1nodes = self.m1.node
        downreach_inds = self.m1.downreach_ind.values
        m2 = self.m2

        print 'Adding new SFR segments and renumbering reaches... (may be slow for large datasets)'
        nfeatures = len(sfr_cells_list)
        for i, new_seg_nodes in enumerate(sfr_cells_list):
            print '\r{:.0f}%'.format(100 * i/nfeatures),
            # skip entries that are empty (e.g. waterbodies that don't intersect SFR)
            if len(new_seg_nodes) == 0:
                continue

            old_segments = np.unique(m1segments[np.in1d(m1nodes, new_seg_nodes)])
            #old_segments = np.unique(self.m1.ix[self.m1.node.isin(new_seg_nodes), 'segment'])

            # iterate through each segment intersecting the feature
            # This code is slow; could probably be greatly sped up by "removing the dots"
            for seg in old_segments:

                #print '{}\t-->\t'.format(seg),
                # reaches that do not intersect the polygon feature
                inds1 = (~np.in1d(m1nodes, new_seg_nodes) & (m1segments == seg))
                #inds1 = ((~self.m1['node'].isin(new_seg_nodes)) & (m1segments == seg)).values

                # reaches that intersect the polygon feature
                inds2 = (np.in1d(m1nodes, new_seg_nodes) & (m1segments == seg))
                #inds2 = (self.m1['node'].isin(new_seg_nodes) & (m1segments == seg)).values

                # renumber the reaches
                newsegs = []
                for i, inds in enumerate([inds1, inds2]):

                    firstnewseg = np.max(m1segments) + 1
                    old_reaches = m1reaches[inds]

                    # in case the whole segment is in the waterbody, continue
                    # (no reaches are in inds1)
                    if len(old_reaches) == 0:
                        continue

                    # renumber these
                    # the segments start at nseg, and then nseg + 1
                    new_segments_inds, new_reaches = renumber(firstnewseg + 0, old_reaches)

                    #print old_reaches
                    m1segments[inds] = new_segments_inds
                    m1reaches[inds] = new_reaches
                    #self.m1.loc[inds, 'segment'] = new_segments_inds
                    #self.m1.loc[inds, 'reach'] = new_reaches
                    newsegs = list(set(newsegs + list(new_segments_inds)))

                # renumber the last new segment to the previous segment number
                # (can't have gaps in segment numbering; kind of clunky, but allows
                # the renumbering algorithm to be general
                lastnewseg = np.max(newsegs)
                m1segments[m1segments == lastnewseg] = seg
                #self.m1.loc[self.m1.segment == lastnewseg, 'segment'] = seg
                newsegs.remove(lastnewseg)

                # update routing
                for newseg in newsegs:
                    #print '{} '.format(newseg),

                    # get the outseg via the index of the next downstream reach 1
                    lastreach = np.max(m1reaches[m1segments == newseg])
                    outseg_lastreach_ind = downreach_inds[(m1segments == newseg) & (m1reaches == lastreach)][0]
                    #lastreach = np.max(m1reaches[m1segments == newseg])
                    #downreach_ind = downreach_inds[(m1segments == newseg) & (m1reaches == outseg_last_reach)][0]
                    #downreach_ind = downreach_inds[(m1segments == newseg) & (m1reaches == lastreach)][0]
                    #downreach_ind = self.m1.ix[(m1segments == newseg) & (m1reaches == lastreach), 'downreach_ind'].values[0]
                    if outseg_lastreach_ind == -999999:
                        outseg = 0
                    else:
                        outseg = m1segments[outseg_lastreach_ind]
                        #outseg = int(self.m1.ix[downreach_ind, 'segment'])

                    m2.loc[newseg, :] = m2.ix[seg, :] # copy all Mat2 info from old seg
                    m2.loc[newseg, ['outseg', 'segment']] = outseg, newseg # update to new segment and outsegment numbers

        self.m1['segment'] = m1segments
        self.m1['reach'] = m1reaches
        self.m2 = m2

        print '\nDone'


class Spatial(SFRdata):

    def __init__(self, sfrobject=None, Mat1=None, Mat2=None, sfr=None, xll=0.0, yll=0.0,
                 GIS_mult=0.3048, prj=None, proj4=None, epsg=None):

        SFRdata.__init__(self, sfrobject=sfrobject, Mat1=Mat1, Mat2=Mat2, sfr=sfr, xll=xll, yll=yll,
                         GIS_mult=GIS_mult, prj=prj, proj4=proj4, epsg=epsg)

    def intersect_with_SFR_cells(self, intersect_df=None, intersect_shapefile=None, intersect_prj=None,
                                 sfr_shapefile=None,
                                 node_attribute=None, GIS_mult=None):

        if GIS_mult is not None:
            self.GIS_mult = GIS_mult

        # get geometries for SFR cells
        if sfr_shapefile is None:

            self.get_cell_geometries()

            self.m1['geometry'] = [self.cell_geometries[n] for n in self.m1.node]
            if self.proj4 is None:
                print 'No coordinate projection supplied for SFR cells.'

        else:
            df = GISio.shp2df(sfr_shapefile, geometry=True)
            df.sort(['segment', 'reach'], inplace=True)

            # check to make sure that SFR shapefile indexing is consistent with Mat1
            diff = np.max(self.m1.reach.values - df.reach.values)

            if diff == 0:
                self.m1['geometry'] = df.geometry.tolist()
            else:
                raise IndexError('Segments and reaches in SFR shapefile are not consistent with Mat1')

            self.proj4 = GISio.get_proj4(sfr_shapefile)

        # read in intersect feature(s)
        if intersect_df is None and intersect_shapefile is not None:
            dfi = GISio.shp2df(intersect_shapefile, geometry=True)
            intersect_name = intersect_shapefile
        elif isinstance(intersect_df, pd.DataFrame):
            dfi = intersect_df
            intersect_name = 'intersect_df'
        else:
            raise ValueError("No intersect feature supplied!")
        print 'Intersecting features in {} with SFR cells'.format(intersect_name)

        # reproject intersect features into coordinate system of model grid (SFR cells)
        if self.proj4 is not None:
            if intersect_shapefile is not None:
                intersect_proj4 = GISio.get_proj4(intersect_shapefile)
            elif intersect_prj is not None:
                intersect_proj4 = GISio.get_proj4(intersect_prj)

            print 'Reprojecting from:\n{}\nto:\n{}\n...'.format(intersect_proj4, self.proj4)
            dfi['geometry'] = GISops.projectdf(dfi, intersect_proj4, self.proj4)
        else:
            print 'SFR cells are assumed to be in same coordinate system as {}.'.format(intersect_name)

        # limit intersect features to those within bounding box for SFR network
        # note that because this uses cell centroids,
        # it may miss intersect features that only nick the outside of the grid
        print 'Discarding features outside of SFR network bounding box...'
        geom = self.m1.geometry.tolist()
        xmin, ymin = np.min([p.centroid.xy[0] for p in geom]), np.min([p.centroid.xy[1] for p in geom])
        xmax, ymax = np.max([p.centroid.xy[0] for p in geom]), np.max([p.centroid.xy[1] for p in geom])
        bbox = Polygon([(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin)])
        dfi = dfi.ix[[f.intersects(bbox) for f in dfi.geometry], :].copy()

        if len(dfi) == 0:
            print 'No intersections between SFR cells and intersect features! Check coordinate projections ' \
                  'and that model origin is in GIS coordinate units.'

        print 'Listing SFR cell intersections for remaining features...'
        sfr_geom = self.m1.geometry.tolist()
        sfr_nodes = self.m1.node.tolist()
        poly_geom = dfi.geometry.tolist()
        try:
            import rtree
            intersections = GISops.intersect_rtree(sfr_geom, poly_geom)
        except:
            print 'Rtree not found. Features will be intersected without spatial indexing ' \
                  '(may take an hour or more for very large problems).'
            intersections = GISops.intersect_brute_force(sfr_geom, poly_geom)

        dfi['sfr_cells'] = [[sfr_nodes[i] for i in p] for p in intersections]
        dfi.to_csv('waterbodies_intersected.csv')
        return dfi
