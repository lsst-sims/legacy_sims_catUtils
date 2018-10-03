import numpy as np
import unittest
import lsst.utils.tests

from lsst.sims.catUtils.mixins import VariabilityAGN


def setup_module(module):
    lsst.utils.tests.init()


class AgnModelTestCase(unittest.TestCase):

    longMessage = True

    def test_agn_structure_function(self):
        """
        Test that the structure_function of the delta_magnitudes
        resulting from our AGN variability model is consistent with
        equation 3 of MacLeod et al. 2010 (ApJ 721, 1014)
        """

        agn_obj = VariabilityAGN()
        n_obj = 10
        d_mjd = 1.0
        rng = np.random.RandomState(88)
        mjd_grid = np.arange(61000.0, 237000.0, d_mjd)
        agn_params = {}
        agn_params['seed'] = rng.randint(10, high=1000, size=n_obj)
        agn_params['agn_tau'] = rng.random_sample(n_obj)*25.0+75.0
        for bp in ('u', 'g', 'r', 'i', 'z', 'y'):
            agn_params['agn_sf%s' % bp] = rng.random_sample(n_obj)*100.0+5.0

        redshift = np.zeros(n_obj, dtype=float)
        dmag_arr = agn_obj.applyAgn([range(n_obj)], agn_params, mjd_grid,
                                    redshift=redshift)

        self.assertEqual(dmag_arr.shape, (6, n_obj, len(mjd_grid)))

        max_dev = -1.0
        for i_obj in range(n_obj):
            tau = agn_params['agn_tau'][i_obj]
            seed = agn_params['seed'][i_obj]
            for i_bp, bp in enumerate(('u', 'g', 'r', 'i', 'z', 'y')):
                sf_inf = agn_params['agn_sf%s' % bp][i_obj]

                # loop over different time lags, calculating the structure
                # function of the light curves and comparing to the
                # expected value of the structure function
                for delta_i_t in range(5,len(mjd_grid)//2, len(mjd_grid)//20):

                    delta_t = d_mjd*delta_i_t

                    dmag_0 = dmag_arr[i_bp][i_obj][:-delta_i_t]
                    dmag_1 = dmag_arr[i_bp][i_obj][delta_i_t:]
                    self.assertEqual(len(dmag_0), len(dmag_1))

                    # expected structure funciton value taken from
                    # equation 3 of MacLeod et al
                    sf_th = sf_inf*np.sqrt(1.0-np.exp(-delta_t/tau))

                    # use definition of structure function from
                    # section 2.1 of Hughes et al. 1992
                    # (ApJ 396, 469)
                    sf_test = np.sqrt(np.mean((dmag_1-dmag_0)**2))

                    # verify that the structure function is within 10%
                    # of the expected value
                    self.assertLess(np.abs(1.0-sf_test/sf_th), 0.1)

    def test_agn_mean(self):
        """
        Test that the mean of time lagged AGN light curves approaches
        delta_magnitude = 0 as the time lag gets larger than
        the AGN variability time scale
        """

        agn_obj = VariabilityAGN()
        n_obj = 10
        d_mjd = 1.0
        rng = np.random.RandomState(11273)
        mjd_grid = np.arange(61000.0, 237000.0, d_mjd)
        agn_params = {}
        agn_params['seed'] = rng.randint(10, high=1000, size=n_obj)
        agn_params['agn_tau'] = rng.random_sample(n_obj)*25.0+75.0
        for bp in ('u', 'g', 'r', 'i', 'z', 'y'):
            agn_params['agn_sf%s' % bp] = rng.random_sample(n_obj)*100.0+5.0

        redshift = np.zeros(n_obj, dtype=float)
        dmag_arr = agn_obj.applyAgn([range(n_obj)], agn_params, mjd_grid,
                                    redshift=redshift)

        self.assertEqual(dmag_arr.shape, (6, n_obj, len(mjd_grid)))

        max_dev = -1.0
        for i_obj in range(n_obj):
            tau = agn_params['agn_tau'][i_obj]
            seed = agn_params['seed'][i_obj]
            for i_bp, bp in enumerate(('u', 'g', 'r', 'i', 'z', 'y')):
                sf_inf = agn_params['agn_sf%s' % bp][i_obj]

                # loop over different time lags, calculating the mean
                # of the light curve taken at those time lags; make
                # sure delta_mag is within 1-sigma of zero
                #
                # only consider lags that are greater than 5*tau
                delta_i_t_min = int(np.round(tau/d_mjd))
                self.assertLess(5*delta_i_t_min, len(mjd_grid)//20)

                ct_lags = 0
                for delta_i_t in range(5*delta_i_t_min, len(mjd_grid)//20, 100):
                    ct_lags += 1
                    t_dexes = range(delta_i_t+rng.randint(0,high=10), len(mjd_grid), delta_i_t)
                    dmag_subset = dmag_arr[i_bp][i_obj][t_dexes]
                    self.assertGreater(len(dmag_subset), 19)
                    dmag_mean = np.mean(dmag_subset)
                    dmag_stdev = np.std(dmag_subset)
                    msg = 'failed with %d samples' % len(dmag_subset)
                    self.assertLess(np.abs(dmag_mean)/dmag_stdev, 1.0,
                                    msg=msg)
                self.assertGreater(ct_lags, 10)

    def test_agn_structure_function_with_redshift(self):
        """
        Test that the structure_function of the delta_magnitudes
        resulting from our AGN variability model is consistent with
        equation 3 of MacLeod et al. 2010 (ApJ 721, 1014)

        This test is done for the case of non-zero redshift
        """

        agn_obj = VariabilityAGN()
        n_obj = 10
        d_mjd = 1.0
        rng = np.random.RandomState(443)
        mjd_grid = np.arange(61000.0, 237000.0, d_mjd)
        agn_params = {}
        agn_params['seed'] = rng.randint(10, high=1000, size=n_obj)
        agn_params['agn_tau'] = rng.random_sample(n_obj)*25.0+75.0
        for bp in ('u', 'g', 'r', 'i', 'z', 'y'):
            agn_params['agn_sf%s' % bp] = rng.random_sample(n_obj)*100.0+5.0

        redshift = rng.random_sample(n_obj)*2.0+0.1
        dmag_arr = agn_obj.applyAgn([range(n_obj)], agn_params, mjd_grid,
                                    redshift=redshift)

        self.assertEqual(dmag_arr.shape, (6, n_obj, len(mjd_grid)))

        max_dev = -1.0
        for i_obj in range(n_obj):
            tau = agn_params['agn_tau'][i_obj]
            seed = agn_params['seed'][i_obj]
            time_dilation = 1.0+redshift[i_obj]
            for i_bp, bp in enumerate(('u', 'g', 'r', 'i', 'z', 'y')):
                sf_inf = agn_params['agn_sf%s' % bp][i_obj]

                # loop over different time lags, calculating the structure
                # function of the light curves and comparing to the
                # expected value of the structure function
                for delta_i_t in range(5,len(mjd_grid)//2, len(mjd_grid)//20):

                    delta_t = d_mjd*delta_i_t/time_dilation

                    dmag_0 = dmag_arr[i_bp][i_obj][:-delta_i_t]
                    dmag_1 = dmag_arr[i_bp][i_obj][delta_i_t:]
                    self.assertEqual(len(dmag_0), len(dmag_1))

                    # expected structure funciton value taken from
                    # equation 3 of MacLeod et al
                    sf_th = sf_inf*np.sqrt(1.0-np.exp(-delta_t/tau))

                    # use definition of structure function from
                    # section 2.1 of Hughes et al. 1992
                    # (ApJ 396, 469)
                    sf_test = np.sqrt(np.mean((dmag_1-dmag_0)**2))

                    # verify that the structure function is within 10%
                    # of the expected value
                    self.assertLess(np.abs(1.0-sf_test/sf_th), 0.1)

    def test_agn_mean_with_redshift(self):
        """
        Test that the mean of time lagged AGN light curves approaches
        delta_magnitude = 0 as the time lag gets larger than
        the AGN variability time scale

        This test is done in the case of non-zero redshift
        """

        agn_obj = VariabilityAGN()
        n_obj = 10
        d_mjd = 1.0
        rng = np.random.RandomState(2273)
        mjd_grid = np.arange(61000.0, 237000.0, d_mjd)
        agn_params = {}
        agn_params['seed'] = rng.randint(10, high=1000, size=n_obj)
        agn_params['agn_tau'] = rng.random_sample(n_obj)*25.0+75.0
        for bp in ('u', 'g', 'r', 'i', 'z', 'y'):
            agn_params['agn_sf%s' % bp] = rng.random_sample(n_obj)*100.0+5.0

        redshift = rng.random_sample(n_obj)*2.0+0.1
        dmag_arr = agn_obj.applyAgn([range(n_obj)], agn_params, mjd_grid,
                                    redshift=redshift)

        self.assertEqual(dmag_arr.shape, (6, n_obj, len(mjd_grid)))

        max_dev = -1.0
        for i_obj in range(n_obj):
            tau = agn_params['agn_tau'][i_obj]
            seed = agn_params['seed'][i_obj]
            time_dilation = 1.0+redshift[i_obj]
            for i_bp, bp in enumerate(('u', 'g', 'r', 'i', 'z', 'y')):
                sf_inf = agn_params['agn_sf%s' % bp][i_obj]

                # loop over different time lags, calculating the mean
                # of the light curve taken at those time lags; make
                # sure delta_mag is within 1-sigma of zero
                #
                # only consider lags that are greater than 5*tau
                delta_i_t_min = int(np.round(tau/(time_dilation*d_mjd)))
                self.assertLess(5*delta_i_t_min, len(mjd_grid)//20)

                ct_lags = 0
                for delta_i_t in range(5*delta_i_t_min, len(mjd_grid)//20, 100):
                    ct_lags += 1
                    t_dexes = range(delta_i_t+rng.randint(0,high=10), len(mjd_grid), delta_i_t)
                    dmag_subset = dmag_arr[i_bp][i_obj][t_dexes]
                    self.assertGreater(len(dmag_subset), 19)
                    dmag_mean = np.mean(dmag_subset)
                    dmag_stdev = np.std(dmag_subset)
                    msg = 'failed with %d samples' % len(dmag_subset)
                    self.assertLess(np.abs(dmag_mean)/dmag_stdev, 1.0,
                                    msg=msg)

                self.assertGreater(ct_lags, 10)

    def test_threading(self):
        """
        Test that running applyAgn with multithreading does not change the answers
        """
        agn_obj = VariabilityAGN()
        agn_obj_2 = VariabilityAGN()
        agn_obj_2._agn_threads = 4
        self.assertEqual(agn_obj._agn_threads, 1)
        self.assertNotEqual(agn_obj._agn_threads, agn_obj_2._agn_threads)

        rng = np.random.RandomState(61422)
        n_agn = 11
        redshift = rng.random_sample(n_agn)*2.0+1.1
        agn_params = {}
        agn_params['agn_tau'] = rng.random_sample(n_agn)*10.0+1.0
        for bp in 'ugrizy':
            agn_params['agn_sf%s' % bp] = rng.random_sample(n_agn)*2.0+0.1
        agn_params['seed'] = rng.randint(2, high=100, size=n_agn)
        mjd = 61923.5
        dmag_control = agn_obj.applyAgn([np.arange(n_agn, dtype=int)],
                                        agn_params, mjd, redshift=redshift)

        self.assertEqual(dmag_control.shape, (6,n_agn))
        n_zero = np.where(np.abs(dmag_control.flatten())<1.0e-10)
        self.assertEqual(len(n_zero[0]), 0)

        dmag_threaded = agn_obj_2.applyAgn([np.arange(n_agn, dtype=int)],
                                           agn_params, mjd, redshift=redshift)

        np.testing.assert_array_equal(dmag_control, dmag_threaded)

        ##### now test it on a numpy array of mjd

        mjd = 59580.0+rng.random_sample(13)*2000.0
        dmag_control = agn_obj.applyAgn([np.arange(n_agn, dtype=int)],
                                        agn_params, mjd, redshift=redshift)

        self.assertEqual(dmag_control.shape, (6,n_agn, len(mjd)))
        n_zero = np.where(np.abs(dmag_control.flatten())<1.0e-10)
        self.assertEqual(len(n_zero[0]), 0)

        dmag_threaded = agn_obj_2.applyAgn([np.arange(n_agn, dtype=int)],
                                           agn_params, mjd, redshift=redshift)

        np.testing.assert_array_equal(dmag_control, dmag_threaded)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass



if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
