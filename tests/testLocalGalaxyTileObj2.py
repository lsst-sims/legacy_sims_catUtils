import unittest
import numpy as np
import sqlite3
import shutil
import tempfile
import os

import lsst.utils.tests
from lsst.utils import getPackageDir

import lsst.sims.utils.htmModule as htm
from lsst.sims.utils import angularSeparation
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import LocalGalaxyModels as LocGal


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


        galtag = (100*np.round(45 + ra_grid/0.05) + np.round(45+dec_grid/0.05)).astype(int)
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

        self.assertGreater(n_data, 10)
        self.assertGreater(n_valid, 1000)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
