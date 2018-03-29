import unittest
import numpy as np

import lsst.utils.tests
from lsst.utils import getPackageDir
from lsst.sims.utils.CodeUtilities import sims_clean_up
from lsst.sims.photUtils import EBVbase


def setup_module(module):
    lsst.utils.tests.init()


class EbvTestCase(unittest.TestCase):

    @classmethod
    def tearDownClass(cls):
        sims_clean_up()

    def test_cache(self):
        """
        Test that EBVbase() only loads each dust map once
        """
        sims_clean_up()
        self.assertEqual(len(EBVbase._ebv_map_cache), 0)

        ebv1 = EBVbase()
        ebv1.load_ebvMapNorth()
        ebv1.load_ebvMapSouth()
        self.assertEqual(len(EBVbase._ebv_map_cache), 2)

        ebv2 = EBVbase()
        ebv2.load_ebvMapNorth()
        ebv2.load_ebvMapSouth()
        self.assertEqual(len(EBVbase._ebv_map_cache), 2)

        rng = np.random.RandomState(881)
        ra_list = rng.random_sample(10)*2.0*np.pi
        dec_list = rng.random_sample(10)*np.pi - 0.5*np.pi

        ebv1_vals = ebv1.calculateEbv(equatorialCoordinates=np.array([ra_list, dec_list]))
        ebv2_vals = ebv2.calculateEbv(equatorialCoordinates=np.array([ra_list, dec_list]))

        self.assertEqual(len(EBVbase._ebv_map_cache), 2)

        np.testing.assert_array_equal(ebv1_vals, ebv2_vals)

    def testEBV(self):

        ebvObject = EBVbase()
        ra = []
        dec = []
        gLat = []
        gLon = []
        for i in range(10):
            ra.append(i*2.0*np.pi/10.0)
            dec.append(i*np.pi/10.0)

            gLat.append(-0.5*np.pi+i*np.pi/10.0)
            gLon.append(i*2.0*np.pi/10.0)

            equatorialCoordinates = np.array([ra, dec])
            galacticCoordinates = np.array([gLon, gLat])

        ebvOutput = ebvObject.calculateEbv(equatorialCoordinates=equatorialCoordinates)
        self.assertEqual(len(ebvOutput), len(ra))

        ebvOutput = ebvObject.calculateEbv(galacticCoordinates=galacticCoordinates)
        self.assertEqual(len(ebvOutput), len(gLon))
        self.assertGreater(len(ebvOutput), 0)

        self.assertRaises(RuntimeError, ebvObject.calculateEbv, equatorialCoordinates=equatorialCoordinates,
                          galacticCoordinates=galacticCoordinates)
        self.assertRaises(RuntimeError, ebvObject.calculateEbv,
                          equatorialCoordinates=None, galacticCoordinates=None)
        self.assertRaises(RuntimeError, ebvObject.calculateEbv)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
