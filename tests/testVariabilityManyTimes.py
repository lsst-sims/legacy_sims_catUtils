import unittest
import lsst.utils.tests
import numpy as np
import copy
import numbers

from lsst.sims.catUtils.mixins import StellarVariabilityModels
from lsst.sims.catUtils.mixins import ExtraGalacticVariabilityModels


def applyAmcvn_original(valid_dexes, params, expmjd_in):
    """
    Copied from VariabilityMixin.py before attempt was made
    to make applyAmcvn handle vectors of MJD.  We will use
    this method to verify that we have not broken the logic
    inside of applyAmcvn.
    """

    if not isinstance(expmjd_in, numbers.Number):
        raise RuntimeError("cannot pass multiple MJD to "
                           "applyAmcvn_original")

    if len(params) == 0:
        return np.array([[],[],[],[],[],[]])

    n_obj = len(params[list(params.keys())[0]])

    maxyears = 10.
    dMag = np.zeros((6, n_obj))
    epoch = expmjd_in

    amplitude = params['amplitude'].astype(float)[valid_dexes]
    t0 = params['t0'].astype(float)[valid_dexes]
    period = params['period'].astype(float)[valid_dexes]
    burst_freq = params['burst_freq'].astype(float)[valid_dexes]
    burst_scale = params['burst_scale'].astype(float)[valid_dexes]
    amp_burst = params['amp_burst'].astype(float)[valid_dexes]
    color_excess = params['color_excess_during_burst'].astype(float)[valid_dexes]
    does_burst = params['does_burst'][valid_dexes]

    # get the light curve of the typical variability
    uLc   = amplitude*np.cos((epoch - t0)/period)
    gLc   = copy.deepcopy(uLc)
    rLc   = copy.deepcopy(uLc)
    iLc   = copy.deepcopy(uLc)
    zLc   = copy.deepcopy(uLc)
    yLc   = copy.deepcopy(uLc)

    # add in the flux from any bursting
    local_bursting_dexes = np.where(does_burst==1)
    for i_burst in local_bursting_dexes[0]:
        adds = 0.0
        for o in np.linspace(t0[i_burst] + burst_freq[i_burst],\
                                t0[i_burst] + maxyears*365.25, \
                                np.ceil(maxyears*365.25/burst_freq[i_burst])):
            tmp = np.exp( -1*(epoch - o)/burst_scale[i_burst])/np.exp(-1.)
            adds -= amp_burst[i_burst]*tmp*(tmp < 1.0)  ## kill the contribution
        ## add some blue excess during the outburst
        uLc[i_burst] += adds +  2.0*color_excess[i_burst]
        gLc[i_burst] += adds + color_excess[i_burst]
        rLc[i_burst] += adds + 0.5*color_excess[i_burst]
        iLc[i_burst] += adds
        zLc[i_burst] += adds
        yLc[i_burst] += adds

    dMag[0][valid_dexes] += uLc
    dMag[1][valid_dexes] += gLc
    dMag[2][valid_dexes] += rLc
    dMag[3][valid_dexes] += iLc
    dMag[4][valid_dexes] += zLc
    dMag[5][valid_dexes] += yLc
    return dMag



def setup_module(module):
    lsst.utils.tests.init()

class StellarVariability_at_many_times_case(unittest.TestCase):
    """
    This test case will verify that all of the stellar variability
    models in CatSim deal correctly with receiving a vector
    of time values.
    """

    longMessage = True

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

    def test_RRLy_many_some_invalid(self):
        """
        Test that the correct thing happens when some of the objects
        do not use the variability model.
        """
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
        valid_dexes = [np.array([1,3])]

        dmag_vector = self.star_var.applyRRly(valid_dexes,
                                              params,
                                              mjd_arr)

        self.assertEqual(dmag_vector.shape, (6, n_obj, n_time))

        for i_time, mjd in enumerate(mjd_arr):
            dmag_test = self.star_var.applyRRly(valid_dexes,
                                                params,
                                                mjd)

            for i_star in range(n_obj):
                if i_star in (0, 2):
                    for i_band in range(6):
                        self.assertEqual(dmag_test[i_band][i_star], 0.0,
                                         msg='failed on obj %d; band %d; time %d' % (i_star, i_band, i_time))
                for i_band in range(6):
                    self.assertEqual(dmag_vector[i_band][i_star][i_time],
                                     dmag_test[i_band][i_star],
                                     msg='failed on obj %d; band %d; time %d' % (i_star, i_band, i_time))


    def test_Cepeheid_many(self):
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

    def test_Cepeheid_many_some_invalid(self):
        """
        Test that the correct thing happens when some of the objects
        do not use the variability model.
        """
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

        valid_dexes = [np.array([1,3])]

        dmag_vector = self.star_var.applyCepheid(valid_dexes, params, mjd_arr)
        self.assertEqual(dmag_vector.shape, (6, n_obj, n_time))

        for i_time, mjd in enumerate(mjd_arr):
            dmag_test = self.star_var.applyCepheid(valid_dexes, params, mjd)
            for i_band in range(6):
                for i_obj in range(n_obj):
                    if i_obj not in (1, 3):
                        self.assertEqual(dmag_test[i_band][i_obj], 0.0)
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

    def test_Eb_many_some_invalid(self):
        """
        Test that the correct thing happens when some of the objects
        do not use the variability model.
        """
        rng = np.random.RandomState(814512)
        params = {}
        params['lcfile'] = ['eb_lc/EB.2294.inp',
                            'eb_lc/EB.1540.inp',
                            'eb_lc/EB.2801.inp']

        n_obj = len(params['lcfile'])
        params['t0'] = rng.random_sample(n_obj)*20000.0+30000.0
        valid_dexes = [np.array([0,2])]
        mjd_arr = rng.random_sample(19)*3653.3+59580.0
        n_time = len(mjd_arr)
        dmag_vector = self.star_var.applyEb(valid_dexes, params, mjd_arr)
        self.assertEqual(dmag_vector.shape, (6, n_obj, n_time))

        for i_time, mjd in enumerate(mjd_arr):
            dmag_test = self.star_var.applyEb(valid_dexes, params, mjd)
            for i_band in range(6):
                for i_obj in range(n_obj):
                    if i_obj not in (0, 2):
                        self.assertEqual(dmag_test[i_band][i_obj], 0.0)
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

    def test_MicroLens_many_some_invalid(self):
        """
        Test that the correct thing happens when some of the objects
        do not use the variability model.
        """
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
        valid_dexes = [np.array([1,3])]
        dmag_vector = self.star_var.applyMicrolens(valid_dexes, params, mjd_arr)
        self.assertEqual(dmag_vector.shape, (6, n_obj, n_time))

        for i_time, mjd in enumerate(mjd_arr):
            dmag_test = self.star_var.applyMicrolens(valid_dexes, params, mjd)
            self.assertEqual(dmag_test.shape, (6, n_obj))
            for i_band in range(6):
                for i_obj in range(n_obj):
                    if i_obj not in (1,3):
                        self.assertEqual(dmag_test[i_band][i_obj], 0.0)
                    self.assertEqual(dmag_test[i_band][i_obj],
                                     dmag_vector[i_band][i_obj][i_time])


    def test_Amcvn_many(self):
        rng = np.random.RandomState(71242)
        n_obj = 20
        doesBurst = rng.randint(0, 2, size=n_obj)
        burst_freq = rng.randint(10, 150, size=n_obj)
        burst_scale = np.array([1500.0]*n_obj)
        amp_burst = rng.random_sample(n_obj)*8.0
        color_excess_during_burst = rng.random_sample(n_obj)*0.2-0.4
        amplitude = rng.random_sample(n_obj)*0.2
        period = rng.random_sample(n_obj)*200.0
        mjDisplacement = rng.random_sample(n_obj)*500.0
        params = {}
        params['does_burst'] = doesBurst
        params['burst_freq'] = burst_freq
        params['burst_scale'] = burst_scale
        params['amp_burst'] = amp_burst
        params['color_excess_during_burst'] = color_excess_during_burst
        params['amplitude'] = amplitude
        params['period'] = period
        params['t0'] = 51500.0-mjDisplacement

        mjd_arr = rng.random_sample(100)*3653.3+59580.0
        n_time = len(mjd_arr)

        valid_dexes = [np.arange(n_obj, dtype=int)]

        dmag_vector = self.star_var.applyAmcvn(valid_dexes, params, mjd_arr)

        self.assertEqual(dmag_vector.shape, (6, n_obj, n_time))
        for i_time, mjd in enumerate(mjd_arr):
            dmag_test = self.star_var.applyAmcvn(valid_dexes, params, mjd)
            self.assertEqual(dmag_test.shape, (6, n_obj))
            dmag_old = applyAmcvn_original(valid_dexes, params,mjd)
            for i_obj in range(n_obj):
                for i_band in range(6):
                    self.assertEqual(dmag_test[i_band][i_obj],
                                     dmag_old[i_band][i_obj])
                    self.assertEqual(dmag_test[i_band][i_obj],
                                     dmag_vector[i_band][i_obj][i_time])

    def test_Amcvn_many_some_invalid(self):
        """
        Test that the correct thing happens when some of the objects
        do not use the variability model.
        """
        rng = np.random.RandomState(71242)
        n_obj = 20
        doesBurst = rng.randint(0, 2, size=n_obj)
        burst_freq = rng.randint(10, 150, size=n_obj)
        burst_scale = np.array([1500.0]*n_obj)
        amp_burst = rng.random_sample(n_obj)*8.0
        color_excess_during_burst = rng.random_sample(n_obj)*0.2-0.4
        amplitude = rng.random_sample(n_obj)*0.2
        period = rng.random_sample(n_obj)*200.0
        mjDisplacement = rng.random_sample(n_obj)*500.0
        params = {}
        params['does_burst'] = doesBurst
        params['burst_freq'] = burst_freq
        params['burst_scale'] = burst_scale
        params['amp_burst'] = amp_burst
        params['color_excess_during_burst'] = color_excess_during_burst
        params['amplitude'] = amplitude
        params['period'] = period
        params['t0'] = 51500.0-mjDisplacement

        mjd_arr = rng.random_sample(100)*3653.3+59580.0
        n_time = len(mjd_arr)

        valid_dexes = [np.array([1,5,6])]

        dmag_vector = self.star_var.applyAmcvn(valid_dexes, params, mjd_arr)

        self.assertEqual(dmag_vector.shape, (6, n_obj, n_time))
        for i_time, mjd in enumerate(mjd_arr):
            dmag_test = self.star_var.applyAmcvn(valid_dexes, params, mjd)
            self.assertEqual(dmag_test.shape, (6, n_obj))
            dmag_old = applyAmcvn_original(valid_dexes, params,mjd)
            for i_obj in range(n_obj):
                for i_band in range(6):
                    if i_obj not in(1,5,6):
                        self.assertEqual(dmag_test[i_band][i_obj], 0.0)
                    self.assertEqual(dmag_test[i_band][i_obj],
                                     dmag_old[i_band][i_obj])
                    self.assertEqual(dmag_test[i_band][i_obj],
                                     dmag_vector[i_band][i_obj][i_time])


    def test_BHMicrolens_many(self):
        rng = np.random.RandomState(5132)
        params = {}
        params['filename'] = ['microlens/bh_binary_source/lc_14_25_75_8000_0_0.05_316',
                              'microlens/bh_binary_source/lc_14_25_4000_8000_0_phi1.09_0.005_100',
                              'microlens/bh_binary_source/lc_14_25_75_8000_0_tets2.09_0.005_316']

        n_obj = len(params['filename'])
        params['t0'] = rng.random_sample(n_obj)*10.0+59580.0

        mjd_arr = rng.random_sample(12)*4.0+59590.0
        n_time = len(mjd_arr)
        valid_dexes = [np.arange(n_obj, dtype=int)]

        dmag_vector = self.star_var.applyBHMicrolens(valid_dexes, params, mjd_arr)
        self.assertEqual(dmag_vector.shape, (6, n_obj, n_time))

        for i_time, mjd in enumerate(mjd_arr):
            dmag_test = self.star_var.applyBHMicrolens(valid_dexes, params, mjd)
            self.assertEqual(dmag_test.shape, (6, n_obj))
            for i_band in range(6):
                for i_obj in range(n_obj):
                    self.assertTrue(isinstance(dmag_test[i_band][i_obj], numbers.Number))
                    self.assertEqual(dmag_test[i_band][i_obj],
                                     dmag_vector[i_band][i_obj][i_time])

    def test_BHMicrolens_many_some_invalid(self):
        """
        Test that the correct thing happens when some of the objects
        do not use the variability model.
        """
        rng = np.random.RandomState(5132)
        params = {}
        params['filename'] = ['microlens/bh_binary_source/lc_14_25_75_8000_0_0.05_316',
                              'microlens/bh_binary_source/lc_14_25_4000_8000_0_phi1.09_0.005_100',
                              'microlens/bh_binary_source/lc_14_25_75_8000_0_tets2.09_0.005_316']

        n_obj = len(params['filename'])
        params['t0'] = rng.random_sample(n_obj)*10.0+59580.0

        mjd_arr = rng.random_sample(12)*4.0+59590.0
        n_time = len(mjd_arr)
        valid_dexes = [np.array([1])]

        dmag_vector = self.star_var.applyBHMicrolens(valid_dexes, params, mjd_arr)
        self.assertEqual(dmag_vector.shape, (6, n_obj, n_time))

        for i_time, mjd in enumerate(mjd_arr):
            dmag_test = self.star_var.applyBHMicrolens(valid_dexes, params, mjd)
            self.assertEqual(dmag_test.shape, (6, n_obj))
            for i_band in range(6):
                for i_obj in range(n_obj):
                    if i_obj != 1:
                        self.assertEqual(dmag_test[i_band][i_obj], 0.0)
                    self.assertTrue(isinstance(dmag_test[i_band][i_obj], numbers.Number))
                    self.assertEqual(dmag_test[i_band][i_obj],
                                     dmag_vector[i_band][i_obj][i_time])


class AgnVariability_at_many_times_case(unittest.TestCase):
    """
    This test case will verify that the AGN variability
    model in CatSim deals correctly with receiving a vector
    of time values.
    """

    longMessage = True

    def test_agn_many(self):
        agn_model = ExtraGalacticVariabilityModels()
        rng = np.random.RandomState(7153)
        params  = {}
        n_obj = 5
        params['agn_tau'] = rng.random_sample(n_obj)*100.0+100.0
        params['agn_sfu'] = rng.random_sample(n_obj)*2.0
        params['agn_sfg'] = rng.random_sample(n_obj)*2.0
        params['agn_sfr'] = rng.random_sample(n_obj)*2.0
        params['agn_sfi'] = rng.random_sample(n_obj)*2.0
        params['agn_sfz'] = rng.random_sample(n_obj)*2.0
        params['agn_sfy'] = rng.random_sample(n_obj)*2.0
        params['t0_mjd'] = 48000.0+rng.random_sample(n_obj)*5.0
        params['seed'] = rng.randint(0, 20000, size=n_obj)
        redshift_arr = rng.random_sample(n_obj)*5.0

        mjd_arr = np.sort(rng.random_sample(6)*3653.3+59580.0)
        n_time = len(mjd_arr)

        valid_dexes = [np.arange(n_obj, dtype=int)]
        dmag_vector = agn_model.applyAgn(valid_dexes, params, mjd_arr,
                                         redshift=redshift_arr)
        self.assertEqual(dmag_vector.shape, (6, n_obj,n_time))

        for i_time, mjd in enumerate(mjd_arr):
            dmag_test = agn_model.applyAgn(valid_dexes, params, mjd,
                                           redshift=redshift_arr)
            self.assertEqual(dmag_test.shape, (6, n_obj))
            for i_band in range(6):
                for i_obj in range(n_obj):
                    self.assertAlmostEqual(dmag_vector[i_band][i_obj][i_time],
                                           dmag_test[i_band][i_obj], 6,
                                           msg='failed on band %d obj %d time %d' % (i_band, i_obj, i_time))

    def test_agn_many_some_invalid(self):
        """
        Test that the correct thing happens when some of the objects
        do not use the variability model.
        """
        agn_model = ExtraGalacticVariabilityModels()
        rng = np.random.RandomState(7153)
        params  = {}
        n_obj = 5
        params['agn_tau'] = rng.random_sample(n_obj)*100.0+100.0
        params['agn_sfu'] = rng.random_sample(n_obj)*2.0
        params['agn_sfg'] = rng.random_sample(n_obj)*2.0
        params['agn_sfr'] = rng.random_sample(n_obj)*2.0
        params['agn_sfi'] = rng.random_sample(n_obj)*2.0
        params['agn_sfz'] = rng.random_sample(n_obj)*2.0
        params['agn_sfy'] = rng.random_sample(n_obj)*2.0
        params['t0_mjd'] = 48000.0+rng.random_sample(n_obj)*5.0
        params['seed'] = rng.randint(0, 20000, size=n_obj)
        redshift_arr = rng.random_sample(n_obj)*5.0

        mjd_arr = np.sort(rng.random_sample(6)*3653.3+59580.0)
        n_time = len(mjd_arr)

        valid_dexes = [np.array([1,2,4])]
        dmag_vector = agn_model.applyAgn(valid_dexes, params, mjd_arr,
                                         redshift=redshift_arr)
        self.assertEqual(dmag_vector.shape, (6, n_obj,n_time))

        for i_time, mjd in enumerate(mjd_arr):
            dmag_test = agn_model.applyAgn(valid_dexes, params, mjd,
                                           redshift=redshift_arr)
            self.assertEqual(dmag_test.shape, (6, n_obj))
            for i_band in range(6):
                for i_obj in range(n_obj):
                    if i_obj not in (1,2,4):
                        self.assertEqual(dmag_test[i_band][i_obj], 0.0)
                    self.assertAlmostEqual(dmag_vector[i_band][i_obj][i_time],
                                           dmag_test[i_band][i_obj], 6,
                                           msg='failed on band %d obj %d time %d' % (i_band, i_obj, i_time))


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
