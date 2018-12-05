"""
This file will contain unit tests that compare the behavior of
LocalGalaxyTileObj and the supporting classes against the stored procedure
on fatboy.

This module does not include the usual LSST memory tests.  This is because
we are testing the connection to fatboy with a global object that will not
be cleaned up before the tests are finished running.  This is ultimately fine
because the dependence of this test on fatboy connectivity means that it will
never actually run on Jenkins.
"""
import unittest
import numpy as np
import time

import lsst.sims.utils as sims_utils
from lsst.sims.catUtils.baseCatalogModels import GalaxyTileObj, GalaxyObj
import lsst.sims.catUtils.baseCatalogModels.LocalGalaxyModels as LocalGalaxy

_can_access_fatboy = True
try:
    _fatboy_tileobj = GalaxyTileObj(database='LSSTCATSIM',
                                    host='fatboy.phys.washington.edu',
                                    port=1433, driver='mssql+pymssql')

    _fatboy_galaxies = GalaxyObj(database='LSSTCATSIM',
                                 host='fatboy.phys.washington.edu',
                                 port=1433, driver='mssql+pymssql')

except:
    _can_access_fatboy = False
    raise


@unittest.skipIf(not _can_access_fatboy, 'No connection to fatboy')
class ChunkIteratorTestCase(unittest.TestCase):

    def setUp(self):
        self.tiles = LocalGalaxy.FatboyTiles()

    def _validate_fatboy_v_local(self, local_galtileid, local_ra, local_dec,
                                 fatboy_galtileid, fatboy_ra, fatboy_dec,
                                 obs):
            # Because of floating point errors, there are going to be some
            # galaxies that fatboy returns which the local classes do not
            # return.  Verify that these are a small fraction of the total
            # galaxies and that they are just barely outside of the relevant
            # sky tiles, according to the local implementation.

            fatboy_tileidx = fatboy_galtileid//100000000
            local_tileidx = local_galtileid//100000000

            local_in_fatboy = np.in1d(local_galtileid,
                                      fatboy_galtileid,
                                      assume_unique=True)

            self.assertTrue(local_in_fatboy.all(),
                            msg="Not all galaxies found by the local "
                                 "implementation were found by fatboy")

            overlap = np.in1d(fatboy_galtileid, local_galtileid,
                              assume_unique=True)
            offenders = np.where(np.logical_not(overlap))
            print('\noffenders %e out of %e' % (len(offenders[0]), len(fatboy_galtileid)))

            for oo in offenders[0]:
                tile_idx = fatboy_tileidx[oo]
                xyz = sims_utils.cartesianFromSpherical(np.radians(fatboy_ra[oo]),
                                                        np.radians(fatboy_dec[oo]))
                sky_tile = self.tiles.tile(tile_idx)
                in_tile = True
                delta_tile_dotproducts = []
                for hs in sky_tile.half_space_list:
                    dotproduct = np.dot(hs.vector, xyz)
                    delta_tile_dotproducts.append(dotproduct-hs.dd)
                    if dotproduct < hs.dd:
                        in_tile = False
                        break

                # now check whether or not the galaxy is in the tile
                # when the tile is rotated to RA=Dec=0
                rot_mat = self.tiles.rotation_matrix(tile_idx)
                xyz_origin = np.dot(rot_mat, xyz)
                origin_tile = sky_tile.rotate(rot_mat)
                for hs in origin_tile.half_space_list:
                    dotproduct = np.dot(hs.vector, xyz_origin)
                    delta_tile_dotproducts.append(dotproduct-hs.dd)
                    if dotproduct < hs.dd:
                        in_tile = False
                        break

                # the np.min(delta_tile_dotproducts) condition has us ignore
                # points that the local implementation thinks are just outside
                # of the tile, likely a result of floating point roundoff errors
                if in_tile and np.min(delta_tile_dotproducts)>1.0e-10:                        
                    # If in_tile is True, check that the angular separation between
                    # the boresite and the galaxy is within the requested bounds;
                    # fatboy returns things that are just barely outside the requested
                    # field of view
                    dd = sims_utils.angularSeparation(obs.pointingRA,
                                                      obs.pointingDec,
                                                      fatboy_ra[oo],
                                                      fatboy_dec[oo])
                    delta_str = '\n'
                    for delta in delta_tile_dotproducts:
                        delta_str += '%e\n' % delta
                    self.assertGreater(dd, obs.boundLength,
                                       msg='One of the galaxies returned by fatboy '
                                      'but not by the local implementation '
                                      '*was* in the expected tile\n' +
                                      delta_str)
                elif not in_tile:
                    # If it was outside the tile, check that it was only just barely
                    # so (i.e. if fatboy and the local implementation *really* disagree
                    # about what is and is not in the requested sky area, we probably
                    # have a problem)
                    self.assertLess(np.abs(np.min(delta_tile_dotproducts)), 2.0e-8,
                                    msg="Not 'just barely' outside of tile")

            overlap = np.where(overlap)
            self.assertLess(len(overlap[0]), len(fatboy_galtileid))
            self.assertLess(len(offenders[0]), len(fatboy_galtileid)/1000)
            np.testing.assert_array_equal(fatboy_galtileid[overlap],
                                          local_galtileid)

            dd = sims_utils.angularSeparation(local_ra, local_dec,
                                              fatboy_ra[overlap],
                                              fatboy_dec[overlap])

            print('\nmax dd %e arcsec' % (3600.0*dd.max()))
            max_dex = np.argmax(dd)
            lra = local_ra[max_dex]
            ldec = local_dec[max_dex]
            fra = fatboy_ra[overlap][max_dex]
            fdec = fatboy_dec[overlap][max_dex]
            print("%e %e" % (lra, ldec))
            print("%e %e" % (fra, fdec))
            print("delta %e %e" % (lra-fra, ldec-fdec))
            l_xyz = sims_utils.cartesianFromSpherical(np.radians(lra), np.radians(ldec))
            f_xyz = sims_utils.cartesianFromSpherical(np.radians(fra), np.radians(fdec))
            print("1.0-dot %e" % (1.0-np.dot(l_xyz,f_xyz)))

            # around the poles there are points that return dd between 1 and 2 milliarcsec;
            # however, if you convert the fatboy and local RA, Dec of those points into
            # Cartesian coordinates and take the dot product of those vectors, you get
            # exactly 1.0
            self.assertLess(dd.max(), 0.002/3600.0)

    def test_ra_dec_galtileid_chunkIterator(self):
        """
        Test that the LocalGalaxyChunkIterator selects the correct galaxies
        and assigns them the correct RA, Dec
        """
        ra_list = [112.0, 154.2, 54.2]
        dec_list = [31.0, -89.7, 89.8]
        radius = 0.4
        chunk_size = 100000
        for ra, dec in zip(ra_list, dec_list):
            obs = sims_utils.ObservationMetaData(pointingRA=ra,
                                                 pointingDec=dec,
                                                 boundType='circle',
                                                 boundLength=radius)


            local_ra = []
            local_dec = []
            local_galtileid = []
            t_start = time.time()
            local_iter = LocalGalaxy.LocalGalaxyChunkIterator(_fatboy_galaxies,
                                                  ['raJ2000', 'decJ2000',
                                                   'ra', 'dec'], obs,
                                                  chunk_size, None)

            ct_chunk = 0
            for chunk in local_iter:
                ct_chunk += 1
                dd = sims_utils.angularSeparation(chunk['ra'], chunk['dec'],
                                                  np.degrees(chunk['raJ2000']),
                                                  np.degrees(chunk['decJ2000']))
                self.assertLess(dd.max(), 0.0005/3600.0)

                local_ra.append(chunk['ra'])
                local_dec.append(chunk['dec'])
                local_galtileid.append(chunk['galtileid'])
            t_local = time.time()-t_start

            self.assertGreater(ct_chunk, 2)
            local_ra = np.concatenate(local_ra)
            local_dec = np.concatenate(local_dec)
            local_galtileid = np.concatenate(local_galtileid)
            sorted_dex = np.argsort(local_galtileid)
            local_galtileid = local_galtileid[sorted_dex]
            local_ra = local_ra[sorted_dex]
            local_dec = local_dec[sorted_dex]

            fatboy_ra = []
            fatboy_dec = []
            fatboy_galtileid = []
            t_start = time.time()
            fatboy_iter = _fatboy_tileobj.query_columns(colnames=['raJ2000',
                                                                 'decJ2000'],
                                                        obs_metadata=obs,
                                                        chunk_size=chunk_size)

            for chunk in fatboy_iter:
                fatboy_ra.append(np.degrees(chunk['raJ2000']))
                fatboy_dec.append(np.degrees(chunk['decJ2000']))
                fatboy_galtileid.append(chunk['galtileid'])
            t_fatboy = time.time()-t_start
            print("t_local %e t_fatboy %e" % (t_local, t_fatboy))

            fatboy_ra = np.concatenate(fatboy_ra)
            fatboy_dec = np.concatenate(fatboy_dec)
            fatboy_galtileid = np.concatenate(fatboy_galtileid)
            sorted_dex = np.argsort(fatboy_galtileid)

            fatboy_ra = fatboy_ra[sorted_dex]
            fatboy_dec = fatboy_dec[sorted_dex]
            fatboy_galtileid = fatboy_galtileid[sorted_dex]

            self._validate_fatboy_v_local(local_galtileid,
                                          local_ra,
                                          local_dec,
                                          fatboy_galtileid,
                                          fatboy_ra,
                                          fatboy_dec,
                                          obs)

    def test_local_galaxy_tile_obj(self):
        """
        Test that LocalGalaxyTileObj returns correct quantities.
        """
        pass


if __name__ == "__main__":
    unittest.main()
