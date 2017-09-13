import unittest
import tempfile
import gzip
import os
import numpy as np

import lsst.utils.tests

from lsst.sims.catUtils.mixins import KeplerLightCurveMixin
from lsst.sims.utils.CodeUtilities import sims_clean_up


def setup_module(module):
    lsst.utils.tests.init()


class KeplerLightCurve_testCase(unittest.TestCase):

    def test_calc_dflux(self):
        """
        Test the method that calculates the flux of
        Kepler light curves by generating a fake light
        curve library with known parameters, calculating
        the fluxes, and comparing to the expected results.
        """
        lc_temp_file_name = tempfile.mktemp(prefix='test_calc_dflux_lc',
                                            suffix='.gz')

        rng = np.random.RandomState(7124)
        n_c_1 = 10
        a1_list = rng.random_sample(n_c_1)*5.0
        b1_list = (rng.random_sample(n_c_1)-0.5)*2.0
        c1_list = (rng.random_sample(n_c_1)-0.5)*0.1
        omega1_list = rng.random_sample(n_c_1)*20.0
        tau1_list = rng.random_sample(n_c_1)*100.0
        median1 = 100.0

        n_c_2 = 15
        a2_list = rng.random_sample(n_c_2)*5.0
        b2_list = (rng.random_sample(n_c_2)-0.5)*2.0
        c2_list = (rng.random_sample(n_c_2)-0.5)*0.1
        omega2_list = rng.random_sample(n_c_2)*20.0
        tau2_list = rng.random_sample(n_c_2)*100.0
        median2 = 200.0

        with gzip.open(lc_temp_file_name, 'w') as out_file:
            out_file.write('# a header\n')
            out_file.write('kplr990000000_lc.txt 100 1.0e+02 %d ' % n_c_1)
            for i_c in range(n_c_1):
                out_file.write('%e ' % (1.0/(i_c+1)))
            out_file.write('%e ' % median1)
            for i_c in range(n_c_1):
                out_file.write('%.15e %.15e %.15e %.15e %.15e ' %
                               (a1_list[i_c], b1_list[i_c], c1_list[i_c],
                                omega1_list[i_c], tau1_list[i_c]))
            out_file.write('\n')

            out_file.write('kplr990000001_lc.txt 100 1.0e+02 %d ' % n_c_2)
            for i_c in range(n_c_2):
                out_file.write('%e ' % (1.0/(i_c+1)))
            out_file.write('%e ' % median2)
            for i_c in range(n_c_2):
                out_file.write('%.15e %.15e %.15e %.15e %.15e ' %
                               (a2_list[i_c], b2_list[i_c], c2_list[i_c],
                                omega2_list[i_c], tau2_list[i_c]))
            out_file.write('\n')

        expmjd = rng.random_sample(100)*200.0
        kp = KeplerLightCurveMixin()
        kp.load_kepler_light_curves(lc_temp_file_name)
        
        q_flux, d_flux = kp._calc_dflux(990000000, expmjd)
        self.assertAlmostEqual(q_flux, median1+c1_list.sum(), 10)
        
        true_flux = np.zeros(len(expmjd))
        for i_c in range(n_c_1):
            arg = omega1_list[i_c]*(expmjd-tau1_list[i_c])
            true_flux += a1_list[i_c]*np.cos(arg)
            true_flux += b1_list[i_c]*np.sin(arg)
        self.assertEqual(len(d_flux), len(true_flux))
        np.testing.assert_allclose(d_flux, true_flux, rtol=0.0, atol=1.0e-10)

        q_flux, d_flux = kp._calc_dflux(990000001, expmjd)
        self.assertAlmostEqual(q_flux, median2+c2_list.sum(), 10)
        
        true_flux = np.zeros(len(expmjd))
        for i_c in range(n_c_2):
            arg = omega2_list[i_c]*(expmjd-tau2_list[i_c])
            true_flux += a2_list[i_c]*np.cos(arg)
            true_flux += b2_list[i_c]*np.sin(arg)
        self.assertEqual(len(d_flux), len(true_flux))
        np.testing.assert_allclose(d_flux, true_flux, rtol=0.0, atol=1.0e-10)

        sims_clean_up()
        if os.path.exists(lc_temp_file_name):
            os.unlink(lc_temp_file_name)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
