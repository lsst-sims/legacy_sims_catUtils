from __future__ import with_statement
import unittest
import os
import numpy as np
import sqlite3
import lsst.utils.tests
from lsst.utils import getPackageDir
from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.catUtils.mixins import VariabilityStars
from lsst.sims.catUtils.mixins import MLTflaringMixin
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catUtils.mixins import PhotometryStars
from lsst.sims.utils import ObservationMetaData

def setup_module(module):
    lsst.utils.tests.init()


class MLT_test_DB(CatalogDBObject):
    objid = 'mlt_test'
    tableid = 'mlt_test'
    idColKey = 'id'
    raColName = 'ra'
    decColName = 'decl'
    objectTypeId = 66

    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('parallax', 'parallax*PI()/648000000.'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', str, 40)]

class MLT_flare_test_case(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.scratch_dir = os.path.join(getPackageDir('sims_catUtils'),
                                       'tests', 'scratchSpace')
        # Create a dummy mlt light curve file
        cls.mlt_lc_name = os.path.join(cls.scratch_dir,
                                       'test_mlt_lc_file.npz')

        lc_files = {}
        amp = 1.0e32
        lc_files['lc_1_time'] = np.arange(0.0, 3652.5, 0.7)
        lc_files['lc_1_u'] = amp*np.power(np.sin(lc_files['lc_1_time']/100.0), 2)
        lc_files['lc_1_g'] = amp*np.power(np.cos(lc_files['lc_1_time']/100.0), 2)

        amp = 2.0e31
        lc_files['lc_2_time'] = np.arange(0.0, 365.25, 1.3)
        lc_files['lc_2_u'] = amp*np.power(np.sin(lc_files['lc_2_time']/50.0), 2)
        lc_files['lc_2_g'] = amp*np.power(np.cos(lc_files['lc_2_time']/50.0), 2)

        with open(cls.mlt_lc_name, 'wb') as file_handle:
            np.savez(file_handle, **lc_files)

        # Create a database of stars using these light curves
        cls.db_name = os.path.join(cls.scratch_dir, 'test_mlt_db.db')

        conn = sqlite3.connect(cls.db_name)
        cursor = conn.cursor()
        cursor.execute('''CREATE TABLE mlt_test
                       (simobjid int, ra real, decl real, sedfilename text,
                        varParamStr text, parallax real, galacticAv real,
                        magNorm real)''')
        conn.commit()

        cursor.execute('''INSERT INTO mlt_test VALUES( 1, 25.0, 31.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz',
                       '{"varMethodName": "applyMLTflaring",
                         "pars": {"lc": "lc_1.txt", "t0": 456.2}}',
                       0.25, 2.4, 17.1)''')

        cursor.execute('''INSERT INTO mlt_test VALUES( 2, 25.2, 32.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz',
                       '{"varMethodName": "applyMLTflaring",
                         "pars": {"lc": "lc_1.txt", "t0": 41006.2}}',
                       0.15, 1.8, 17.2)''')

        cursor.execute('''INSERT INTO mlt_test VALUES( 3, 25.3, 10.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz',
                       '{"varMethodName": "applyMLTflaring",
                         "pars": {"lc": "lc_2.txt", "t0": 117.2}}',
                       0.3, 2.6, 17.3)''')

        cursor.execute('''INSERT INTO mlt_test VALUES( 4, 25.4, 11.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz',
                       '{"varMethodName": "applyMLTflaring",
                         "pars": {"lc": "lc_2.txt", "t0": 10456.2}}',
                       0.22, 2.3, 17.4)''')
        conn.commit()
        conn.close()


    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.mlt_lc_name):
            os.unlink(cls.mlt_lc_name)

        if os.path.exists(cls.db_name):
            os.unlink(cls.db_name)

    def test_flare_magnitudes(self):
        """
        Test that we get the expectedmagnitudes out
        """
        db = MLT_test_DB(database=self.db_name, driver='sqlite')

        obs = ObservationMetaData(mjd=60000.0)

        class QuiescentCatalog(PhotometryStars, InstanceCatalog):
            column_outputs = ['id', 'lsst_u', 'lsst_g']

        quiet_cat = QuiescentCatalog(db, obs_metadata=obs)
        quiet_cat_name = os.path.join(self.scratch_dir, 'mlt_quiet_cat.txt')
        quiet_cat.write_catalog(quiet_cat_name)

        class FlaringCatalog(PhotometryStars, VariabilityStars,
                             MLTflaringMixin, InstanceCatalog):
            column_outputs = ['id', 'lsst_u', 'lsst_g']

        db.show_db_columns()

        flare_cat = FlaringCatalog(db, obs_metadata=obs)
        flare_cat._mlt_lc_file = self.mlt_lc_name
        flare_cat_name = os.path.join(self.scratch_dir, 'mlt_flaring_cat.txt')
        flare_cat.write_catalog(flare_cat_name)

class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
