import unittest
import numpy as np
import sqlite3
import shutil
import tempfile
import os

import lsst.utils.tests
from lsst.utils import getPackageDir

import lsst.sims.utils.htmModule as htm
from lsst.sims.catUtils.baseCatalogModels import LocalGalaxyModels as LocGal


def setup_module(module):
    lsst.utils.tests.init()


class GalaxyTileObjTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        dir_name = os.path.join(getPackageDir('sims_catUtils'), 'tests')
        cls._tmpdir = tempfile.mkdtemp(dir=dir_name)

        cls._temp_gal_db = os.path.join(cls._tmpdir, 'test_galaxies.db')

        d_pos = 0.05
        ra_grid = np.arange(-2.25, 2.251, d_pos)
        dec_grid = np.arange(-2.25, 2.251, d_pos)
        print('raw grid %d' % len(ra_grid))
        ra_dec_mesh = np.meshgrid(ra_grid, dec_grid)
        ra_grid = ra_dec_mesh[0].flatten()
        dec_grid = ra_dec_mesh[1].flatten()
        galtag = 100*(45 + ra_grid/0.05) + (45+dec_grid/0.05)
        htmid_grid = htm.findHtmid(ra_grid, dec_grid, 21)
        print('got htmid %d' % len(htmid_grid))
        gid = np.arange(len(ra_grid), dtype=int)
        assert len(galtag) == len(np.unique(galtag))
        print(ra_grid.max(),ra_grid.min())

        with sqlite3.connect(cls._temp_gal_db) as conn:
            c = conn.cursor()
            query = '''CREATE TABLE galaxy(htmid int, galid int,
                       ra real, dec real, galtag int)'''
            c.execute(query).fetchall()
            values = ((hh, ii, r, d, g) for hh, ii, r, d, g in
                      zip(htmid_grid, gid, ra_grid, dec_grid, galtag))
            c.executemany('INSERT INTO galaxy VALUES (?,?,?,?,?)', values)

    @classmethod
    def tearDownClass(cls):
        if os.path.isdir(cls._tmpdir):
            shutil.rmtree(cls._tmpdir)

    def test_method(self):
        pass


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
