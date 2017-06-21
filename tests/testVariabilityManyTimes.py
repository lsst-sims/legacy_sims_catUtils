import unittest
import lsst.utils.tests
import numpy as np

from lsst.sims.catUtils.mixins import StellarVariabilityModels


def setup_module(module):
    lsst.utils.tests.init()

class Variability_at_many_times_case(unittest.TestCase):
    """
    This test case will verify that all of the variability
    models in CatSim deal correctly with receiving a vector
    of time values.
    """

    def setUp(self):
        self.star_var = StellarVariabilityModels()
        self.star_var.initializeVariability()

    def test_RRLy_many(self):
        rng = np.random.RandomState(42)
        params = {}
        params['filename'] = ['rrly_lc/RRc/959802_per.txt',
                              'rrly_lc/RRc/1078860_per.txt',
                              'rrly_lc/RRab/98874_per.txt',
                              'rrly_lc/RRab/3879827_per.txt']

        n_obj = len(params['filename'])

        params['tStartMjd'] = rng.random_sample(n_obj)*1000.0+40000.0

        mjd_arr = rng.random_sample(10)*3653.3+59580.0

        n_time = len(mjd_arr)

        dmag_vector = self.star_var.applyRRly([np.arange(n_obj, dtype=int)],
                                              params,
                                              mjd_arr)

        self.assertEqual(dmag_vector.shape, (6, n_obj, n_time))

        for i_time, mjd in enumerate(mjd_arr):
            dmag_test = self.star_var.applyRRly([np.arange(n_obj,dtype=int)],
                                                params,
                                                mjd)

            for i_star in range(n_obj):
                for i_band in range(6):
                    self.assertEqual(dmag_vector[i_band][i_star][i_time],
                                     dmag_test[i_band][i_star])

    def test_Cepeheid_may(self):
        rng = np.random.RandomState(8123)
        params = {}
        params['lcfile'] = ['cepheid_lc/classical_longPer_specfile',
                            'cepheid_lc/classical_medPer_specfile',
                            'cepheid_lc/classical_shortPer_specfile',
                            'cepheid_lc/classical_shortPer_specfile',
                            'cepheid_lc/popII_longPer_specfile',
                            'cepheid_lc/popII_shortPer_specfile']

        n_obj = len(params['lcfile'])

        params['t0'] = rng.random_sample(n_obj)*1000.0+40000.0

        mjd_arr = rng.random_sample(10)*3653.3 + 59580.0

        n_time = len(mjd_arr)

        valid_dexes = [np.arange(n_obj, dtype=int)]

        dmag_vector = self.star_var.applyCepheid(valid_dexes, params, mjd_arr)
        self.assertEqual(dmag_vector.shape, (6, n_obj, n_time))

        for i_time, mjd in enumerate(mjd_arr):
            dmag_test = self.star_var.applyCepheid(valid_dexes, params, mjd)
            for i_band in range(6):
                for i_obj in range(n_obj):
                    self.assertEqual(dmag_test[i_band][i_obj],
                                     dmag_vector[i_band][i_obj][i_time])

    def test_Eb_many(self):
        rng = np.random.RandomState(814512)
        params = {}
        params['lcfile'] = ['eb_lc/EB.2294.inp',
                            'eb_lc/EB.1540.inp',
                            'eb_lc/EB.2801.inp']

        n_obj = len(params['lcfile'])
        params['t0'] = rng.random_sample(n_obj)*20000.0+30000.0
        valid_dexes = [np.arange(n_obj, dtype=int)]
        mjd_arr = rng.random_sample(19)*3653.3+59580.0
        n_time = len(mjd_arr)
        dmag_vector = self.star_var.applyEb(valid_dexes, params, mjd_arr)
        self.assertEqual(dmag_vector.shape, (6, n_obj, n_time))

        for i_time, mjd in enumerate(mjd_arr):
            dmag_test = self.star_var.applyEb(valid_dexes, params, mjd)
            for i_band in range(6):
                for i_obj in range(n_obj):
                    self.assertEqual(dmag_test[i_band][i_obj],
                                     dmag_vector[i_band][i_obj][i_time])

    def test_MicroLens_many(self):
        rng = np.random.RandomState(65123)
        n_obj = 5
        that = rng.random_sample(n_obj)*40.0+40.0
        umin = rng.random_sample(n_obj)
        mjDisplacement = rng.random_sample(n_obj)*50.0
        params = {}
        params['that'] = that
        params['umin'] = umin
        params['t0'] = mjDisplacement+51580.0

        mjd_arr = rng.random_sample(17)*3653.3+59580.0
        n_time = len(mjd_arr)
        valid_dexes = [np.arange(n_obj, dtype=int)]
        dmag_vector = self.star_var.applyMicrolens(valid_dexes, params, mjd_arr)
        self.assertEqual(dmag_vector.shape, (6, n_obj, n_time))

        for i_time, mjd in enumerate(mjd_arr):
            dmag_test = self.star_var.applyMicrolens(valid_dexes, params, mjd)
            self.assertEqual(dmag_test.shape, (6, n_obj))
            for i_band in range(6):
                for i_obj in range(n_obj):
                    self.assertEqual(dmag_test[i_band][i_obj],
                                     dmag_vector[i_band][i_obj][i_time])


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
