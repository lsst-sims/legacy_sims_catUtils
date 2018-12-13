import unittest
import os
import sqlite3
import lsst.utils.tests
from lsst.utils import getPackageDir
import lsst.sims.utils.htmModule as htm
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator


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
        obs_gen = ObservationMetaDataGenerator(db_file)
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

    def test_htmid_query(self):
        """
        Test that queries on HTMID work
        """
        db_dir = os.path.join(getPackageDir('sims_data'), 'OpSimData')
        db_file = os.path.join(db_dir, 'astro-lsst-01_2014.db')
        obs_gen = ObservationMetaDataGenerator(db_file)

        boundLength = 0.7
        htmid  = 36277

        tx = htm.trixelFromHtmid(htmid)
        radius = tx.get_radius()
        ra_c, dec_c = tx.get_center()

        # construct a set of rectangular bounds
        # that will encompass the entire trixel
        # plus extra pointings.  This will be our
        # "control" list
        d_ra_dec = 3.50+2.0*(radius+boundLength)
        ra_min = ra_c - d_ra_dec
        ra_max = ra_c + d_ra_dec
        dec_min = dec_c - d_ra_dec
        dec_max = dec_c + d_ra_dec

        obs_control_list = obs_gen.getObservationMetaData(fieldRA=(ra_min, ra_max),
                                                          fieldDec=(dec_min, dec_max),
                                                          boundLength=boundLength,
                                                          boundType='circle')

        control_id_list = [obs.OpsimMetaData['observationId']
                           for obs in obs_control_list]
        obs_control_dict = {}
        for obs in obs_control_list:
            obs_control_dict[obs.OpsimMetaData['observationId']] = obs

        obs_htmid_list = obs_gen.obsFromHtmid(htmid,
                                              boundLength=boundLength,
                                              transform_to_degrees=None)

        htmid_id_list = [obs.OpsimMetaData['observationId']
                         for obs in obs_htmid_list]


        obs_htmid_dict = {}
        for obs in obs_htmid_list:
            obs_htmid_dict[obs.OpsimMetaData['observationId']] = obs

        self.assertGreater(len(htmid_id_list), 5)
        self.assertGreater(len(control_id_list), len(htmid_id_list))

        # make sure that every pointing selected on htmid is
        # in the control set
        control_id_set = set(control_id_list)
        for hh in htmid_id_list:
            self.assertIn(hh, control_id_set)

        htmid_id_set = set(htmid_id_list)
        for cc in control_id_list:
            obs = obs_control_dict[cc]
            hs = htm.halfSpaceFromRaDec(obs.pointingRA, obs.pointingDec, obs.boundLength)
            if cc not in htmid_id_set:
                # if a control pointing is not in the set selected on htmid,
                # make sure it does not overlap the trixel
                self.assertEqual(hs.contains_trixel(tx), 'outside')
            else:
                # make sure that pointings in both sets overlap the trixel
                # and are identical
                self.assertNotEqual(hs.contains_trixel(tx), 'outside')
                self.assertEqual(obs_control_dict[cc], obs_htmid_dict[cc])


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
