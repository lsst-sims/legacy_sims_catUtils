import unittest
import os
import sqlite3
import lsst.utils.tests
from lsst.utils import getPackageDir
from lsst.sims.catUtils.utils import ObservationMetaDataGeneratorV4


def setup_module(module):
    lsst.utils.tests.init()


class ObsGenV4TestCase(unittest.TestCase):

    def test_spatial_query(self):
        """
        Test that spatial queries work
        """
        db_dir = os.path.join(getPackageDir('sims_data'), 'OpSimData')
        assert os.path.isdir(db_dir)
        db_file = os.path.join(db_dir, 'astro-lsst-01_2014.db')
        obs_gen = ObservationMetaDataGeneratorV4(db_file)
        obs_list = obs_gen.getObservationMetaData(fieldRA=(20.0, 40.0),
                                                  fieldDec=(-30.0, -10.0))
        self.assertGreater(len(obs_list), 10)
        with sqlite3.connect(db_file) as conn:
            cursor = conn.cursor()
            query = '''SELECT observationId, fieldRA, fieldDec,
                       observationStartMJD, filter
                       FROM SummaryAllProps WHERE
                       fieldRA BETWEEN 20.0 AND 40.0 AND
                       fieldDec BETWEEN -30.0 AND -10.0
                       ORDER BY observationId'''
            control = cursor.execute(query).fetchall()
        self.assertEqual(len(control), len(obs_list))
        for ii in range(len(obs_list)):
            self.assertEqual(obs_list[ii].OpsimMetaData['observationId'],
                             int(control[ii][0]))

            self.assertAlmostEqual(obs_list[ii].pointingRA,
                                   float(control[ii][1]), 10)
            self.assertAlmostEqual(obs_list[ii].pointingDec,
                                   float(control[ii][2]), 10)
            self.assertAlmostEqual(obs_list[ii].mjd.TAI,
                                   float(control[ii][3]), 7)
            self.assertEqual(obs_list[ii].bandpass,
                             str(control[ii][4]))

            self.assertGreaterEqual(obs_list[ii].pointingRA, 20.0)
            self.assertLessEqual(obs_list[ii].pointingRA, 40.0)
            self.assertGreaterEqual(obs_list[ii].pointingDec, -30.0)
            self.assertLessEqual(obs_list[ii].pointingDec, -10.0)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
