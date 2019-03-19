import unittest
import numpy as np
import sqlite3
import shutil
import tempfile
import os

import lsst.utils.tests
from lsst.utils import getPackageDir

import lsst.sims.utils.htmModule as htm
from lsst.sims.utils import xyz_from_ra_dec
from lsst.sims.utils import ra_dec_from_xyz
from lsst.sims.utils import angularSeparation
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import LocalGalaxyModels as LocGal

from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catalogs.decorators import cached


ROOT = os.path.abspath(os.path.dirname(__file__))


def setup_module(module):
    lsst.utils.tests.init()


class GalaxyTileObjTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._tile_radius = angularSeparation(0.0, 0.0, 2.0, 2.0)

        dir_name = os.path.join(getPackageDir('sims_catUtils'), 'tests')
        cls._tmpdir = tempfile.mkdtemp(dir=dir_name)

        cls._temp_gal_db = os.path.join(cls._tmpdir, 'test_galaxies.db')

        ra_min = -2.25
        cls._ra_min = ra_min
        ra_max = 2.251
        cls._d_ra = 0.05
        ra_grid = np.arange(ra_min, ra_max, cls._d_ra)
        dec_grid = np.arange(ra_min, ra_max, cls._d_ra)
        print('raw grid %d' % len(ra_grid))
        ra_dec_mesh = np.meshgrid(ra_grid, dec_grid)
        ra_grid = ra_dec_mesh[0].flatten()
        dec_grid = ra_dec_mesh[1].flatten()

        # add a very small offset so that numerical precision
        # does not foul things up on the tile boundaries
        rng = np.random.RandomState(7163)
        ra_grid += 1.0e-5*(rng.random_sample(len(ra_grid))-0.5)
        dec_grid += 1.0e-5*(rng.random_sample(len(dec_grid))-0.5)


        galtag = (100*np.round(45 + ra_grid/cls._d_ra) +
                  np.round(45+dec_grid/cls._d_ra)).astype(int)
        assert len(galtag) == len(np.unique(galtag))
        htmid_grid = htm.findHtmid(ra_grid, dec_grid, 21)
        print('got htmid %d' % len(htmid_grid))
        print(htm.levelFromHtmid(htmid_grid.min()))
        print(htm.levelFromHtmid(htmid_grid.max()))
        assert htm.levelFromHtmid(htmid_grid.min())==21
        assert htm.levelFromHtmid(htmid_grid.max())==21
        gid = np.arange(len(ra_grid), dtype=int)
        assert len(galtag) == len(np.unique(galtag))
        print(ra_grid.max(),ra_grid.min())

        with sqlite3.connect(cls._temp_gal_db) as conn:
            c = conn.cursor()
            query = '''CREATE TABLE galaxy(htmid int, id int,
                       galid text, ra real, dec real, galtag int)'''
            c.execute(query).fetchall()
            values = ((int(hh), int(ii), str(ii), r, d, int(g)) for hh, ii, r, d, g in
                      zip(htmid_grid, gid, ra_grid, dec_grid, galtag))
            c.executemany('INSERT INTO galaxy VALUES (?,?,?,?,?,?)', values)

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(cls._tmpdir):
            shutil.rmtree(cls._tmpdir)


    def test_equatorial_tile(self):
        """
        Query an observation at the center of an equatorial tile; make sure the right galaxies
        get returned
        """

        # select all of the galaxies that should be in the tile
        with sqlite3.connect(self._temp_gal_db) as db_conn:
            c = db_conn.cursor()
            results = c.execute("SELECT galtag FROM galaxy WHERE ra>=-2.0 AND ra<=2.0 "
                                "AND dec>=-2.0 AND dec<=2.0")

            expected_galtag = set(int(r[0]) for r in results)

        dbobj = LocGal.LocalGalaxyTileObj(database=self._temp_gal_db, driver='sqlite')
        obs = ObservationMetaData(pointingRA=34.0, pointingDec=0.0,
                                  boundType='circle',
                                  boundLength=self._tile_radius)

        data_iter = dbobj.query_columns(['galtileid', 'ra', 'dec', 'galtag'],
                                        obs_metadata=obs)

        n_data = 0
        n_valid = 0
        for chunk in data_iter:
            n_data += len(chunk)

            # select only those galaxies that are inside the tile on which we
            # are centered
            valid = np.where(np.logical_and(chunk['ra']>=obs.pointingRA-2.0,
                             np.logical_and(chunk['ra']<=obs.pointingRA+2.0,
                             np.logical_and(chunk['dec']>=obs.pointingDec-2.0,
                                            chunk['dec']<=obs.pointingDec+2.0))))
            n_valid += len(valid[0])
            if len(valid[0])==0:
                continue

            # validate that the correct galaxies were rotated to the correct position
            # when the galaxy database was rotated into the field of view
            valid_data = chunk[valid]
            tag_shld = 100*np.round(45+(valid_data['ra']-obs.pointingRA)/self._d_ra).astype(int)
            tag_shld += np.round(45+(valid_data['dec']-obs.pointingDec)/self._d_ra).astype(int)
            np.testing.assert_array_equal(tag_shld, valid_data['galtag'])
            for tag in valid_data['galtag']:
                self.assertIn(tag, expected_galtag)
            for tag in expected_galtag:
                self.assertIn(tag, valid_data['galtag'])

        self.assertGreater(n_data, 10)
        self.assertGreater(n_valid, 1000)

    def test_two_equatorial_tiles(self):
        """
        Query a field of view directly between two equatorial tiles and make sure that
        all expected galaxies are returned where we think they ought to be
        """
        radius = 1.0
        expected_galaxies_left = []  # galaxies that come from left edge of original tile
        expected_galaxies_right = []  # galaxies that come from right edge of original tile
        with sqlite3.connect(self._temp_gal_db) as db_conn:
            c = db_conn.cursor()
            query = "SELECT ra, dec, galtag FROM galaxy "
            query += "WHERE ra>=-2.0 AND ra<=2.0 "
            query += "AND dec>=-2.0 AND dec<=2.0"
            results = c.execute(query).fetchall()
            for gal in results:
                dd_left = angularSeparation(-2.0, 0.0, gal[0], gal[1])
                if dd_left <= radius:
                    expected_galaxies_left.append(gal)
                dd_right = angularSeparation(2.0, 0.0, gal[0], gal[1])
                if dd_right <= radius:
                    expected_galaxies_right.append(gal)

        # create sets of the expected galtag values
        gal_tag_left = set(g[2] for g in expected_galaxies_left)
        gal_tag_right = set(g[2] for g in expected_galaxies_right)
        self.assertGreater(len(gal_tag_left), 100)
        self.assertGreater(len(gal_tag_right), 100)

        # construct field of view on the border between two equatorial tiles
        ra1 = 34.0
        ra2 = 38.0
        obs = ObservationMetaData(pointingRA=0.5*(ra1+ra2), pointingDec=0.0,
                                  boundType='circle', boundLength=radius)

        dbobj = LocGal.LocalGalaxyTileObj(database=self._temp_gal_db, driver='sqlite')
        data_iter = dbobj.query_columns(['galtileid', 'ra', 'dec', 'galtag'],
                                        obs_metadata=obs)

        found_left = set()
        found_right = set()
        for chunk in data_iter:
            valid_left = np.where(chunk['ra']>obs.pointingRA)
            valid_right = np.where(chunk['ra']<obs.pointingRA)
            self.assertEqual(len(valid_left[0])+len(valid_right[0]), len(chunk))
            data_left = chunk[valid_left]
            tag_left = 100*np.round(45+(data_left['ra']-ra2)/self._d_ra)
            tag_left += np.round(45+(data_left['dec']-obs.pointingDec)/self._d_ra)
            tag_left = tag_left.astype(int)
            for tag in tag_left:
                self.assertIn(tag, gal_tag_left)
                found_left.add(tag)

            data_right = chunk[valid_right]
            tag_right = 100*np.round(45+(data_right['ra']-ra1)/self._d_ra)
            tag_right += np.round(45+(data_right['dec']-obs.pointingDec)/self._d_ra)
            tag_right = tag_right.astype(int)
            for tag in tag_right:
                self.assertIn(tag, gal_tag_right)
                found_right.add(tag)

        # make sure that the galaxies returned by the query
        # match what we expected
        for tag in found_left:
            self.assertIn(tag, gal_tag_left)
        for tag in gal_tag_left:
            self.assertIn(tag, found_left)
        for tag in found_right:
            self.assertIn(tag, gal_tag_right)
        for tag in gal_tag_right:
            self.assertIn(tag, found_right)

    def test_four_corners(self):
        """
        Test field of view centered at a corner between four tiles
        (not on the equator)
        """
        ra_obs = 44.0
        dec_obs = 34.0
        radius = 2.8  # this should be enough to get galaxies to appear
                      # in more than one tile

        # construct bounding half spaces for the four quadrants
        quadrant_half_spaces = []
        quadrant_centers = []  # prime tile coords of tile corners
        for ra_q, dec_q in zip([46.0, 42.0, 42.0, 46.0],
                                    [36.0, 36.0, 32.0, 32.0]):

            cos_ra = np.cos(np.radians(ra_q))
            sin_ra = np.sin(np.radians(ra_q))
            cos_dec = np.cos(np.radians(dec_q))
            sin_dec = np.sin(np.radians(dec_q))

            m_ra = np.array([[cos_ra, sin_ra, 0.0],
                             [-sin_ra, cos_ra, 0.0],
                             [0.0, 0.0, 1.0]])

            m_dec = np.array([[cos_dec, 0.0, sin_dec],
                              [0.0, 1.0, 0.0],
                              [-sin_dec, 0.0, cos_dec]])

            upper_hs_vv = np.array([0.0, 0.0, -1.0])
            upper_hs_vv = np.dot(m_dec,
                                 np.dot(m_ra, upper_hs_vv))
            rr, dd = ra_dec_from_xyz(upper_hs_vv[0],
                                     upper_hs_vv[1],
                                     upper_hs_vv[2])
            upper_hs = htm.halfSpaceFromRaDec(rr, dd, 90.0+dec_q+2.0)

            lower_hs_vv = np.array([0.0, 0.0, 1.0])
            lower_hs_vv = np.dot(m_dec,
                                 np.dot(m_ra, lower_hs_vv))
            rr, dd = ra_dec_from_xyz(lower_hs_vv[0],
                                     lower_hs_vv[1],
                                     lower_hs_vv[2])

            lower_hs = htm.halfSpaceFromRaDec(rr, dd, 92.0-dec_q)

            left_hs_vv = xyz_from_ra_dec(ra_q+88.0, 0.0)
            left_hs_vv = np.dot(m_dec,
                                np.dot(m_ra, left_hs_vv))
            rr, dd = ra_dec_from_xyz(left_hs_vv[0],
                                     left_hs_vv[1],
                                     left_hs_vv[2])
            left_hs = htm.halfSpaceFromRaDec(rr, dd, 90.0)

            right_hs_vv = xyz_from_ra_dec(ra_q-88.0, 0.0)
            right_hs_vv = np.dot(m_dec,
                                 np.dot(m_ra, right_hs_vv))
            rr, dd = ra_dec_from_xyz(right_hs_vv[0],
                                     right_hs_vv[1],
                                     right_hs_vv[2])
            right_hs = htm.halfSpaceFromRaDec(rr, dd, 90.0)
            quadrant_half_spaces.append((upper_hs, lower_hs,
                                         left_hs, right_hs))

            corner = xyz_from_ra_dec(44.0, 34.0)
            rot_corner = np.dot(m_dec,
                                np.dot(m_ra, corner))
            quadrant_centers.append(ra_dec_from_xyz(rot_corner[0],
                                                    rot_corner[1],
                                                    rot_corner[2]))


        gal_tag_1st_quad = set()
        gal_tag_2nd_quad = set()
        gal_tag_3rd_quad = set()
        gal_tag_4th_quad = set()
        n_multiple_quad = 0
        with sqlite3.connect(self._temp_gal_db) as db_conn:
            c = db_conn.cursor()
            query = "SELECT ra, dec, galtag FROM galaxy "
            results = c.execute(query).fetchall()
            for gal in results:
                n_quad = 0
                vv = xyz_from_ra_dec(gal[0], gal[1])
                dd = list([angularSeparation(c[0], c[1], gal[0], gal[1])
                           for c in quadrant_centers])

                if gal[2]== 6785:
                    for hs in quadrant_half_spaces[3]:
                        print('6785 contained ',hs.contains_pt(vv))

                if (dd[0] <= radius and
                    quadrant_half_spaces[0][0].contains_pt(vv) and
                    quadrant_half_spaces[0][1].contains_pt(vv) and
                    quadrant_half_spaces[0][2].contains_pt(vv) and
                    quadrant_half_spaces[0][3].contains_pt(vv)):
                    n_quad += 1
                    gal_tag_1st_quad.add(gal[2])

                if (dd[1] <= radius and
                    quadrant_half_spaces[1][0].contains_pt(vv) and
                    quadrant_half_spaces[1][1].contains_pt(vv) and
                    quadrant_half_spaces[1][2].contains_pt(vv) and
                    quadrant_half_spaces[1][3].contains_pt(vv)):
                    n_quad += 1
                    gal_tag_2nd_quad.add(gal[2])

                if (dd[2] <= radius and
                    quadrant_half_spaces[2][0].contains_pt(vv) and
                    quadrant_half_spaces[2][1].contains_pt(vv) and
                    quadrant_half_spaces[2][2].contains_pt(vv) and
                    quadrant_half_spaces[2][3].contains_pt(vv)):
                    n_quad += 1
                    gal_tag_3rd_quad.add(gal[2])

                if (dd[3] <= radius and
                    quadrant_half_spaces[3][0].contains_pt(vv) and
                    quadrant_half_spaces[3][1].contains_pt(vv) and
                    quadrant_half_spaces[3][2].contains_pt(vv) and
                    quadrant_half_spaces[3][3].contains_pt(vv)):
                    n_quad += 1
                    gal_tag_4th_quad.add(gal[2])

                if n_quad>1:
                    n_multiple_quad += 1

        self.assertGreater(n_multiple_quad, 100)

        gal_found_1st_quad = set()
        gal_found_2nd_quad = set()
        gal_found_3rd_quad = set()
        ra_dec_3 = {}
        gal_found_4th_quad = set()
        ra_dec_4 = {}

        obs = ObservationMetaData(pointingRA=ra_obs, pointingDec=dec_obs,
                                  boundType='circle', boundLength=radius)

        dbobj = LocGal.LocalGalaxyTileObj(database=self._temp_gal_db, driver='sqlite')
        data_iter = dbobj.query_columns(['galtileid', 'ra', 'dec', 'galtag'],
                                        obs_metadata=obs)

        for chunk in data_iter:
            for gal in chunk:
                if gal['ra'] >= ra_obs and gal['dec'] >= dec_obs:
                    quad_set = gal_found_1st_quad
                elif gal['ra'] < ra_obs and gal['dec'] >= dec_obs:
                    quad_set = gal_found_2nd_quad
                elif gal['ra'] < ra_obs and gal['dec'] < dec_obs:
                    quad_set = gal_found_3rd_quad
                    ra_dec_3[gal['galtag']] = (gal['ra'], gal['dec'])
                elif gal['ra'] >= ra_obs and gal['dec'] < dec_obs:
                    quad_set = gal_found_4th_quad
                    ra_dec_4[gal['galtag']] = (gal['ra'], gal['dec'])
                else:
                    raise RuntimeError("Unsure what quadrant galaxy belongs in")
                assert gal['galtag'] not in quad_set
                quad_set.add(gal['galtag'])

        test_sum = 0
        control_sum = 0
        for test, control, ra_dec in zip([gal_found_1st_quad, gal_found_2nd_quad,
                                  gal_found_3rd_quad, gal_found_4th_quad],
                                 [gal_tag_1st_quad, gal_tag_2nd_quad,
                                  gal_tag_3rd_quad, gal_tag_4th_quad],
                                 [None, None, ra_dec_3, ra_dec_4]):

            n_erroneous_test = 0
            for tag in test:
                if tag not in control:
                    n_erroneous_test += 1
                    dd = angularSeparation(ra_obs, dec_obs,
                            ra_dec[tag][0], ra_dec[tag][1])
                    print('    bad test %d %e -- %.4f %.4f' % (tag, dd,
                    ra_dec[tag][0], ra_dec[tag][1]))
            n_erroneous_control = 0
            for tag in control:
                if tag not in test:
                    n_erroneous_control += 1
            print('n_test %d bad %d' % (len(test), n_erroneous_test))
            print('n_control %d bad %d\n' % (len(control), n_erroneous_control))
            test_sum += len(test)
            control_sum += len(control)
        print('test_sum %d' % test_sum)
        print('control_sum %d' % control_sum)

        self.assertEqual(len(gal_found_1st_quad), len(gal_tag_1st_quad))
        self.assertEqual(len(gal_found_2nd_quad), len(gal_tag_2nd_quad))
        self.assertEqual(len(gal_found_3rd_quad), len(gal_tag_3rd_quad))
        self.assertEqual(len(gal_found_4th_quad), len(gal_tag_4th_quad))

        for tag in gal_found_1st_quad:
            self.assertIn(tag, gal_tag_1st_quad)
        for tag in gal_found_2nd_quad:
            self.assertIn(tag, gal_tag_2nd_quad)
        for tag in gal_found_3rd_quad:
            self.assertIn(tag, gal_tag_3rd_quad)
        for tag in gal_found_4th_quad:
            self.assertIn(tag, gal_tag_4th_quad)

    def test_local_instance_catalog(self):
        """
        Test that the local galaxy obj can interface with the
        InstanceCatalog as expected
        """

        class LocalGalDummyICat(InstanceCatalog):
            column_outputs=['raJ2000', 'decJ2000',
                            'ra', 'dec', 'galtag',
                            'someFloat']

            delimiter = ' '
            default_formats = {'f': '%.9g'}

            @cached
            def get_someFloat(self):
                t = self.column_by_name('galtag')
                r = self.column_by_name('ra')
                d = self.column_by_name('dec')
                return r**2+d+0.5*t

        obs = ObservationMetaData(pointingRA=34.0, pointingDec=44.0,
                                  boundType='circle', boundLength=3.0)

        dbobj = LocGal.LocalGalaxyTileObj(database=self._temp_gal_db, driver='sqlite')

        cat = LocalGalDummyICat(dbobj, obs_metadata=obs)
        scratch_dir = tempfile.mkdtemp(dir=ROOT, prefix='local_gal')
        cat_name = tempfile.mktemp(dir=scratch_dir,
                                   prefix='local_gal',
                                   suffix='.txt')

        cat.write_catalog(cat_name)
        dt = np.dtype([('raJ2000', float), ('decJ2000', float),
                       ('ra', float), ('dec', float), ('galtag', int),
                       ('someFloat', float)])

        data = np.genfromtxt(cat_name, dtype=dt)
        dd = angularSeparation(data['ra'], data['dec'],
                               np.degrees(data['raJ2000']),
                               np.degrees(data['decJ2000']))
        self.assertLess(dd.max(), 0.001/3600.0)
        ff = data['ra']**2+data['dec']+0.5*data['galtag']
        np.testing.assert_array_almost_equal(ff, data['someFloat'], decimal=5)

        if os.path.exists(cat_name):
            os.unlink(cat_name)
        shutil.rmtree(scratch_dir)

class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
