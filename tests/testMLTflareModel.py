from __future__ import with_statement
import unittest
import os
import numpy as np
import sqlite3
import lsst.utils.tests
from lsst.utils import getPackageDir
from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.catUtils.mixins import VariabilityStars
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catUtils.mixins import PhotometryStars

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
               ('magNorm', '(-2.5*log(flux_scale)/log(10.)) - 18.402732642'),
               ('parallax', 'parallax*PI()/648000000.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
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
        lc_files['lc_1_u'] = amp*np.sin(lc_files['lc_1_time']/100.0)
        lc_files['lc_1_g'] = amp*np.cos(lc_files['lc_1_time']/100.0)

        amp = 2.0e31
        lc_files['lc_2_time'] = np.arange(0.0, 365.25, 1.3)
        lc_files['lc_2_u'] = amp*np.sin(lc_files['lc_2_time']/50.0)
        lc_files['lc_2_g'] = amp*np.cos(lc_files['lc_2_time']/50.0)

        with open(cls.mlt_lc_name, 'wb') as file_handle:
            np.savez(file_handle, lc_files)

        # Create a database of stars using these light curves
        cls.db_name = os.path.join(cls.scratch_dir, 'test_mlt_db.db')

        conn = sqlite3.connect(cls.db_name)
        cursor = conn.cursor()
        cursor.execute('''CREATE TABLE mlt_test
                       (simobjid int, ra real, decl real, sedfilename text,
                        varParamStr text, parallax real, ebv,
                        flux_scale)''')
        conn.commit()

        cursor.execute('''INSERT INTO mlt_test VALUES( 1, 25.0, 31.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz',
                       '{"varMethodName": "applyMLTflaring",
                         "pars": {"lc": "lc_1.txt", "t0": 456.2}}',
                       0.25, 0.698, 1.2e-18)''')

        cursor.execute('''INSERT INTO mlt_test VALUES( 2, 25.2, 32.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz',
                       '{"varMethodName": "applyMLTflaring",
                         "pars": {"lc": "lc_1.txt", "t0": 41006.2}}',
                       0.15, 0.714, 1.3e-18)''')

        cursor.execute('''INSERT INTO mlt_test VALUES( 3, 25.3, 10.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz',
                       '{"varMethodName": "applyMLTflaring",
                         "pars": {"lc": "lc_2.txt", "t0": 117.2}}',
                       0.3, 0.816, 2.2e-18)''')

        cursor.execute('''INSERT INTO mlt_test VALUES( 4, 25.4, 11.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz',
                       '{"varMethodName": "applyMLTflaring",
                         "pars": {"lc": "lc_2.txt", "t0": 10456.2}}',
                       0.22, 0.798, 1.5e-18)''')
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
        pass


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
