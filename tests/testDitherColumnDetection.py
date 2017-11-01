from builtins import zip
from builtins import range
import unittest
import os
import sqlite3
import numpy as np
import tempfile
import shutil

import lsst.utils.tests
from lsst.utils import getPackageDir
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.utils.CodeUtilities import sims_clean_up

ROOT = os.path.abspath(os.path.dirname(__file__))


def setup_module(module):
    lsst.utils.tests.init()


class ObsMetaDataGenDitherTestClass(unittest.TestCase):
    """
    This TestCase will verify that the ObservationMetaDataGenerator
    puts summary columns which are not hardcoded into its interface
    into the OpsimMetaData of the ObservationMetaData it creates.
    """

    @classmethod
    def tearDownClass(cls):
        sims_clean_up()
        if os.path.exists(cls.fake_db_name):
            os.unlink(cls.fake_db_name)
        if os.path.exists(cls.scratch_space):
            shutil.rmtree(cls.scratch_space)

    @classmethod
    def setUpClass(cls):
        cls.scratch_space = tempfile.mkdtemp(dir=ROOT, prefix='ObsMetaDataGenDitherTestClass-')

        cls.fake_db_name = os.path.join(cls.scratch_space,
                                        'dither_test_fake_opsim_sqlite.db')

        conn = sqlite3.connect(cls.fake_db_name)
        curs = conn.cursor()
        curs.execute('''CREATE TABLE Summary (fieldRA real, fieldDec real,
                      obsHistID int, rotSkyPos real, m5 real,
                      raTestDithering real, decTestDithering real, expMJD real,
                      filter text)''')

        conn.commit()

        n_ptngs = 10
        rng = np.random.RandomState(18341)
        ra_list = rng.random_sample(n_ptngs)*2.0*np.pi
        dec_list = rng.random_sample(n_ptngs)*np.pi-0.5*np.pi
        rotSkyPos_list = rng.random_sample(n_ptngs)*np.pi
        m5_list = rng.random_sample(n_ptngs)*20.0
        ra_dither_list = ra_list + rng.random_sample(n_ptngs)*0.1+0.1
        dec_dither_list = dec_list + rng.random_sample(n_ptngs)*0.1+0.1
        expMJD_list = rng.random_sample(n_ptngs)*1000.0

        cls.db_control = []

        for obsid, (ra, dec, rot, m5, raDith, decDith, mjd) in \
        enumerate(zip(ra_list, dec_list, rotSkyPos_list, m5_list,
                      ra_dither_list, dec_dither_list, expMJD_list)):

            cls.db_control.append({'ra': ra, 'dec': dec, 'rot': rot, 'm5': m5,
                                   'raDith': raDith, 'decDith': decDith, 'mjd': mjd})
            curs.execute('''INSERT INTO Summary VALUES
                         (%.12f, %.12f, %d, %.12f, %.12f, %.12f, %.12f, %.12f, '%s')''' %
                         (ra, dec, obsid, rot, m5, raDith, decDith, mjd, 'g'))

        conn.commit()
        conn.close()

    def test_query(self):
        """
        Use ObservationMetaData to query an OpSim-like database that contains
        dithering columns.  Make sure that the dithering columns get carried
        over into the OpsimMetaData of the resulting ObservationMetaData.
        """

        gen = ObservationMetaDataGenerator(database=self.fake_db_name, driver='sqlite')
        obs_list = gen.getObservationMetaData(fieldRA=(0.0, 180.0))
        self.assertGreater(len(obs_list), 0)
        found_list = []
        for obs in obs_list:
            obsid = obs.OpsimMetaData['obsHistID']
            control_dict = self.db_control[obsid]
            self.assertAlmostEqual(obs._pointingRA, control_dict['ra'], 11)
            self.assertAlmostEqual(obs._pointingDec, control_dict['dec'], 11)
            self.assertAlmostEqual(obs._rotSkyPos, control_dict['rot'], 11)
            self.assertAlmostEqual(obs.OpsimMetaData['m5'], control_dict['m5'], 11)
            self.assertAlmostEqual(obs.OpsimMetaData['raTestDithering'], control_dict['raDith'], 11)
            self.assertAlmostEqual(obs.OpsimMetaData['decTestDithering'], control_dict['decDith'], 11)
            self.assertAlmostEqual(obs.mjd.TAI, control_dict['mjd'], 11)
            self.assertEqual(obs.bandpass, 'g')
            self.assertGreaterEqual(obs.pointingRA, 0.0)
            self.assertLessEqual(obs.pointingRA, 180.0)
            found_list.append(obs.OpsimMetaData['obsHistID'])

        # check that the entries not returned do, in fact, violate the query
        for ix in range(len(self.db_control)):
            if ix not in found_list:
                self.assertGreater(self.db_control[ix]['ra'], np.radians(180.0))


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
