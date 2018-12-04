import os
import numpy as np
import numpy.lib.recfunctions as np_recfn
import json
import copy
from sqlalchemy import text
from lsst.utils import getPackageDir
from lsst.sims.utils import HalfSpace, levelFromHtmid
from lsst.sims.utils import halfSpaceFromRaDec
from lsst.sims.utils import halfSpaceFromPoints
from lsst.sims.utils import intersectHalfSpaces
from lsst.sims.utils import cartesianFromSpherical, sphericalFromCartesian
from lsst.sims.catalogs.db import ChunkIterator
from lsst.sims.catUtils.baseCatalogModels import GalaxyTileObj

__all__ = ["FatboyTiles"]


class Tile(object):

    def __init__(self, box_corners):
        self._trixel_bounds = None
        self._trixel_bound_level = None
        self._hs_list = []
        if len(box_corners) == 0:
            return
        self._init_from_corners(box_corners)

    def _init_from_corners(self, box_corners):
        ra_range = [c[0] for c in box_corners]
        ra_min = min(ra_range)
        ra_max = max(ra_range)
        dec_range = [c[1] for c in box_corners]
        dec_min = min(dec_range)
        dec_max = max(dec_range)
        #print(tile_id,dec_min,dec_max)
        tol = 1.0e-10
        for i_c1 in range(len(box_corners)):
            c1 = box_corners[i_c1]
            pt1 = cartesianFromSpherical(np.degrees(c1[0]), np.degrees(c1[1]))
            for i_c2 in range(i_c1+1, len(box_corners), 1):
                hs = None
                c2 = box_corners[i_c2]
                pt2 = cartesianFromSpherical(np.degrees(c2[0]), np.degrees(c2[1]))
                if np.abs(1.0-np.dot(pt1, pt2))<tol:
                    continue

                dra = np.abs(c1[0]-c2[0])
                ddec = np.abs(c1[1]-c2[1])

                if dra<tol and ddec>tol:
                    # The RAs of the two corners is identical, but the Decs are
                    # different; this Half Space is defined by a Great Circle
                    if np.abs(c1[0]-ra_min)<tol:
                        inner_pt = (ra_min+0.001, dec_min+0.001)
                    else:
                        inner_pt = (ra_max-0.001, dec_min+0.001)
                    hs = halfSpaceFromPoints(c1, c2, inner_pt)
                elif ddec<tol and dra>tol:
                    # The Decs of the two corners is identical, bu the RAs are
                    # different; this Half Space is defined by a line of constant
                    # Dec and should be centered at one of the poles
                    if np.abs(c1[1]-dec_min)<tol:
                        hs = halfSpaceFromRaDec(0.0, 90.0, 90.0-dec_min)
                    else:
                        hs = halfSpaceFromRaDec(0.0, -90.0, 90.0+dec_max)
                else:
                    continue

                if hs is None:
                    raise RuntimeError("Somehow Half Space == None")
                self._hs_list.append(hs)

    def contains_many_pts(self, pts):
        result = None
        for hs in self.half_space_list:
            valid = hs.contains_many_pts(pts)
            if result is None:
                result = valid
            else:
                result &= valid
        return result

    @property
    def half_space_list(self):
        return self._hs_list

    def rotate(self, matrix):
        new_tile = Tile([])
        for hs in self.half_space_list:
            vv = np.dot(matrix, hs.vector)
            new_hs = HalfSpace(vv, hs.dd)
            new_tile._hs_list.append(new_hs)
        return new_tile

    def intersects_circle(self, center_pt, radius_rad):
        gross_is_contained = True
        for hs in self.half_space_list:
            if not hs.intersects_circle(center_pt, radius_rad):
                gross_is_contained = False
                break
        if not gross_is_contained:
            return False

        hs_interest = HalfSpace(center_pt, np.cos(radius_rad))
        for i_h1 in range(len(self.half_space_list)):
            hs1 = self.half_space_list[i_h1]
            roots = intersectHalfSpaces(hs1, hs_interest)
            if len(roots) == 0:
                continue

            for i_h2 in range(len(self.half_space_list)):
                if i_h1 == i_h2:
                    continue
                hs2 = self.half_space_list[i_h2]
                local_contained = False
                for rr in roots:
                    if hs2.contains_pt(rr):
                        local_contained = True
                        break
                if not local_contained:
                    return False
        return True

    def _generate_all_trixels(self, level):
        output = None
        for hs in self.half_space_list:
            local_limits = hs.findAllTrixels(level)
            if output is None:
                output = local_limits
            else:
                output = HalfSpace.join_trixel_bound_sets(output, local_limits)
        self._trixel_bounds = output
        self._trixel_bound_level = level
        return None

    @property
    def trixel_bound_level(self):
        return self._trixel_bound_level

    @property
    def trixel_bounds(self):
        return self._trixel_bounds

    def find_all_trixels(self, level):
        if self._trixel_bounds is None or self.trixel_bound_level != level:
            self._generate_all_trixels(level)
        return self._trixel_bounds


class FatboyTiles(object):
    """
    A class to store the fatboy galaxy tiles as a series of Half Spaces
    """

    def __init__(self):
        data_dir = os.path.join(getPackageDir('sims_catUtils'), 'data')
        data_file = os.path.join(data_dir, 'tile_data.txt')
        dtype = np.dtype([('id', int), ('ra', float), ('dec', float),
                          ('box', str, 500)])
        tile_data = np.genfromtxt(data_file, dtype=dtype, delimiter=';')

        self._tile_id = tile_data['id']
        self._tile_ra = {}
        self._tile_dec = {}
        self._rotation_matrix_dict = {}
        for ii, rr, dd, in zip(tile_data['id'], tile_data['ra'], tile_data['dec']):
            self._tile_ra[ii] = rr
            self._tile_dec[ii] = dd
            ra_rad = np.radians(rr)
            dec_rad = np.radians(dd)

            ra_mat = np.array([[np.cos(ra_rad), np.sin(ra_rad), 0.0],
                               [-np.sin(ra_rad), np.cos(ra_rad), 0.0],
                               [0.0, 0.0, 1.0]])

            dec_mat = np.array([[np.cos(dec_rad), 0.0, np.sin(dec_rad)],
                                [0.0, 1.0, 0.0],
                                [-np.sin(dec_rad), 0.0, np.cos(dec_rad)]])

            full_mat = np.dot(dec_mat, ra_mat)
            self._rotation_matrix_dict[ii] = full_mat

        self._tile_dict = {}
        for tile_id, box in zip(tile_data['id'], tile_data['box']):
            box_corners = json.loads(box)
            self._tile_dict[tile_id] = Tile(box_corners)


    def tile_ra(self, tile_idx):
        return self._tile_ra[tile_idx]

    def tile_dec(self, tile_idx):
        return self._tile_dec[tile_idx]

    def rotation_matrix(self, tile_idx):
        return self._rotation_matrix_dict[tile_idx]

    def tile(self, tile_idx):
        return self._tile_dict[tile_idx]

    def find_all_tiles(self, ra, dec, radius):
        """
        ra, dec, radius are all in degrees

        returns a numpy array of tile IDs that intersect the circle
        """
        valid_id = []
        radius_rad = np.radians(radius)
        center_pt = cartesianFromSpherical(np.radians(ra), np.radians(dec))
        for tile_id in self._tile_dict:
            tile = self._tile_dict[tile_id]
            is_contained = tile.intersects_circle(center_pt, radius_rad)
            if is_contained:
                valid_id.append(tile_id)

        return np.array(valid_id)


class LocalGalaxyChunkIterator(ChunkIterator):

    def __init__(self, dbobj, colnames, obs_metadata, chunk_size, constraint):
        """
        Parameters
        ----------
        dbobj -- a CatalogDBObject connected to the 'galaxies' table on fatboy

        colnames -- a list of the columns to query

        chunk_size -- size of chunks to return

        constraint -- a string specifying a SQL 'WHERE' clause
        """
        self.arbitrarySQL = False
        self.dbobj = dbobj
        if 'ra' not in colnames:
            query_colnames = ['htmid', 'galid', 'ra', 'dec'] + colnames
        else:
            query_colnames = ['htmid', 'galid'] + colnames
        self._column_query = dbobj._get_column_query(query_colnames)
        self.chunk_size = chunk_size
        tile_idx_list = np.sort(self._find_tiles(obs_metadata))
        self._trixel_search_level = 9
        self.obs_metadata = obs_metadata
        total_trixel_bounds = []
        self._00_bounds = []
        self._rotate_to_sky = []
        self._sky_tile = []
        self._tile_idx = []

        # construct a HalfSpace based on obs_metadata
        self.obs_hs = halfSpaceFromRaDec(obs_metadata.pointingRA,
                                         obs_metadata.pointingDec,
                                         obs_metadata.boundLength)

        for tile_idx in tile_idx_list:
            rotate_to_00 = self.fatboy_tiles.rotation_matrix(tile_idx)

            # find the bounds for trixels contained by the field of view
            # when rotated from the current tile to RA=Dec=0
            new_vv = np.dot(rotate_to_00, self.obs_hs.vector)
            obs_hs_00 = HalfSpace(new_vv, self.obs_hs.dd)
            obs_hs_00_trixels = obs_hs_00.findAllTrixels(self._trixel_search_level)

            # find the trixels in the current tile when it is rotated
            # to RA=Dec=0
            sky_tile = self.fatboy_tiles.tile(tile_idx)
            single_tile = sky_tile.rotate(rotate_to_00)
            local_bounds = single_tile.find_all_trixels(self._trixel_search_level)
            local_bounds = HalfSpace.join_trixel_bound_sets(local_bounds, obs_hs_00_trixels)

            total_trixel_bounds += local_bounds

            self._sky_tile.append(sky_tile)
            self._00_bounds.append(local_bounds)
            self._rotate_to_sky.append(np.linalg.inv(rotate_to_00))
            self._tile_idx.append(tile_idx)

        print('last pass on trixel_bounds')
        total_trixel_bounds = HalfSpace.merge_trixel_bounds(total_trixel_bounds)
        print('time to write where clause')

        where_clause = "("
        for bound in total_trixel_bounds:
            if where_clause != "(":
                where_clause += " OR "
            htmid_min = bound[0] << 2*(21-self._trixel_search_level)
            htmid_max = (bound[1]+1) << 2*(21-self._trixel_search_level)
            assert levelFromHtmid(htmid_min) == 21
            assert levelFromHtmid(htmid_max) == 21
            assert htmid_min<htmid_max
            where_clause += "(htmid>=%d AND htmid<=%d)" % (htmid_min, htmid_max)
        where_clause += ")"
        print('got where clause')

        if constraint is not None:
            where_clause += "AND (%s)" % constraint

        query = self._column_query
        query = query.filter(text(where_clause))

        print(query)

        self._galaxy_query = dbobj.connection.session.execute(query)
        self._tile_to_do = 0

        self._has_J2000 = False
        if 'raJ2000' in colnames:
            self._has_J2000 = True


    def __next__(self):
        #print('running on tile %d of %d' % (self._tile_to_do, len(self._rotate_to_sky)))
        if self._tile_to_do == 0:
            if self.chunk_size is None and not self._galaxy_query.closed:
                results = self._galaxy_query.fetchall()
            elif self.chunk_size is not None:
                results = self._galaxy_query.fetchmany(self.chunk_size)
            else:
                raise StopIteration
            self._galaxy_cache = self.dbobj._convert_results_to_numpy_recarray_dbobj(results)
            if len(self._galaxy_cache) == 0:
                raise StopIteration

        #print("galaxy_cache is ",type(self._galaxy_cache),self._galaxy_cache['htmid'].min())
        current_chunk = copy.deepcopy(self._galaxy_cache)
        rot_mat = self._rotate_to_sky[self._tile_to_do]
        bounds = self._00_bounds[self._tile_to_do]
        sky_tile = self._sky_tile[self._tile_to_do]
        tile_idx = self._tile_idx[self._tile_to_do]

        make_the_cut = None
        for bb in bounds:
            htmid_min = bb[0] << 2*(21-self._trixel_search_level)
            htmid_max = (bb[1]+1) << 2*(21-self._trixel_search_level)
            valid = ((current_chunk['htmid']>=htmid_min) & (current_chunk['htmid']<=htmid_max))
            if make_the_cut is None:
                make_the_cut = valid
            else:
                make_the_cut |= valid

        good_dexes = np.where(make_the_cut)[0]
        if len(good_dexes) < len(current_chunk):
            current_chunk = current_chunk[good_dexes]

        self._tile_to_do += 1
        if self._tile_to_do >= len(self._rotate_to_sky):
            self._tile_to_do = 0

        if len(current_chunk) == 0:
            return self.__next__()

        #print(current_chunk)
        xyz = cartesianFromSpherical(np.radians(current_chunk['ra']),
                                     np.radians(current_chunk['dec']))

        xyz_sky = np.dot(rot_mat, xyz.transpose()).transpose()

        final_cut = sky_tile.contains_many_pts(xyz_sky)
        final_cut &= self.obs_hs.contains_many_pts(xyz_sky)
        final_cut = np.where(final_cut)

        xyz_sky = xyz_sky[final_cut]
        current_chunk = current_chunk[final_cut]
        if len(current_chunk) == 0:
            return self.__next__()

        ra_dec_sky = sphericalFromCartesian(xyz_sky)
        current_chunk['ra'] = np.degrees(ra_dec_sky[0]) % 360.0
        current_chunk['dec'] = np.degrees(ra_dec_sky[1]) % 360.0
        current_chunk['dec'] = np.where(current_chunk['dec']<270.0,
                                        current_chunk['dec'],
                                        current_chunk['dec']-360.0)
        current_chunk['dec'] = np.where(np.abs(current_chunk['dec'])<=90.0,
                                        current_chunk['dec'],
                                        180.0-current_chunk['dec'])
        if self._has_J2000:
            current_chunk['raJ2000'] = ra_dec_sky[0] % (2.0*np.pi)
            _dec = ra_dec_sky[1] % (2.0*np.pi)
            current_chunk['decJ2000'] = np.where(_dec<1.5*np.pi,
                                                 _dec,
                                                 _dec-2.0*np.pi)
            current_chunk['decJ2000'] = np.where(np.abs(current_chunk['decJ2000'])<=0.5*np.pi,
                                                 current_chunk['decJ2000'],
                                                 np.pi-current_chunk['decJ2000'])


        #print('current_chunk is ',type(current_chunk))

        #>>> r2 = recfunc.append_fields(r,['d','e'],d,dtypes=[float, int], usemask=False, asrecarray=True)

        galtileid = tile_idx*100000000+current_chunk['id']
        current_chunk = np_recfn.append_fields(current_chunk, ['galtileid'], [galtileid],
                                               dtypes=[int], usemask=False, asrecarray=True)

        return self._postprocess_results(current_chunk)

    @property
    def fatboy_tiles(self):
        if not hasattr(self, '_fatboy_tiles'):
            self._fatboy_tiles = FatboyTiles()
        return self._fatboy_tiles

    def _find_tiles(self, obs_metadata):

        if obs_metadata.boundType != 'circle':
            raise RuntimeError("Cannot use ObservationMetaData with "
                               "boundType == %s in LocalGalaxyTileObj" % obs_metadata.boundType)

        return self.fatboy_tiles.find_all_tiles(obs_metadata.pointingRA,
                                                obs_metadata.pointingDec,
                                                obs_metadata.boundLength)



class LocalGalaxyTileObj(GalaxyTileObj):

    def query_columns(self, colnames=None, chunk_size=None, obs_metadata=None, constraint=None,
                      limit=None):
        """Execute a query

        **Parameters**

            * colnames : list or None
              a list of valid column names, corresponding to entries in the
              `columns` class attribute.  If not specified, all columns are
              queried.
            * chunksize : int (optional)
              if specified, then return an iterator object to query the database,
              each time returning the next `chunksize` elements.  If not
              specified, all matching results will be returned.
            * obs_metadata : object (optional)
              object containing information on the observation including the region of the sky
              to query and time of the observation.
            * constraint : string (optional)
              if specified, the predicate is added to the query verbatim using AND
            * limit:
              This kwarg is not actually used.  It exists to preserve the same interface
              as other definitions of query_columns elsewhere in CatSim.  If not None,
              a warning will be emitted, pointing out to the user that 'limit' is not used.

        **Returns**

            * result : structured array or iterator
              If chunksize is not specified, then result is a structured array of all
              items which match the specified query with columns named by the column
              names in the columns class attribute.  If chunksize is specified,
              then result is an iterator over structured arrays of the given size.

        """

        if colnames is None:
            colnames = [k for k in self.columnMap.keys()]

        # We know that galtileid comes back with the query, but we don't want
        # to add it to the query since it's generated on the fly.
        #
        # 25 August 2015
        # The code below has been modified to remove all column names
        # that contain 'galtileid.'  This is to accommodate the
        # CompoundInstanceCatalog and CompoundDBObject classes, which
        # mangle column names such that they include the objid of the
        # specific CatalogDBObject that is asking for them.
        for name in colnames:
            if 'galtileid' in name:
                colnames.remove(name)

        mappedcolnames = ["%s as %s"%(self.columnMap[x], x) for x in colnames]
        mappedcolnames = ",".join(mappedcolnames)

        tile_id_list = self._find_tiles(obs_metadata)

        #if constraint is not None:
        #    query += ", @WhereClause = '%s'"%(constraint)

        if limit is not None:
            warnings.warn("You specified a row number limit in your query of a GalaxyTileObj "
                          "daughter class.  Because of the way GalaxyTileObj is searched, row "
                          "number limits are not possible.  If you really want to limit the number "
                          "of rows returned by you query, consider using GalaxyObj (note that you "
                          "will have to you limit your search to -2.5<RA<2.5 -2.25<Dec<2.25 -- both in "
                          "degrees -- as this is the only region where galaxies exist in GalaxyObj).")

        # should probably write a new ChunkIterator that will do the query once
        # and then selectively munge the outputs per relevant tile
        return LocalGalaxyChunkIterator(self, query, chunk_size, constraint)
