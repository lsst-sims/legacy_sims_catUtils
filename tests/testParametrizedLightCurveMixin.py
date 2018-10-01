import unittest
import tempfile
import gzip
import os
import numpy as np

import lsst.utils.tests

from lsst.sims.catUtils.mixins import ParametrizedLightCurveMixin
from lsst.sims.catUtils.mixins import VariabilityStars
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catalogs.db import fileDBObject
from lsst.sims.utils import ObservationMetaData
from lsst.sims.utils.CodeUtilities import sims_clean_up


def setup_module(module):
    lsst.utils.tests.init()


class ParametrizedLightCurve_testCase(unittest.TestCase):

    def test_calc_dflux(self):
        """
        Test the method that calculates the flux of
        parametrized light curves by generating a fake light
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
            out_file.write(b'# a header\n')
            out_file.write(b'kplr990000000_lc.txt 100 1.0e+02 %d ' % n_c_1)
            for i_c in range(n_c_1):
                out_file.write(b'%e ' % (1.0/(i_c+1)))
            out_file.write(b'%e ' % median1)
            for i_c in range(n_c_1):
                out_file.write(b'%.15e %.15e %.15e %.15e %.15e ' %
                               (a1_list[i_c], b1_list[i_c], c1_list[i_c],
                                omega1_list[i_c], tau1_list[i_c]))
            out_file.write(b'\n')

            out_file.write(b'kplr990000001_lc.txt 100 1.0e+02 %d ' % n_c_2)
            for i_c in range(n_c_2):
                out_file.write(b'%e ' % (1.0/(i_c+1)))
            out_file.write(b'%e ' % median2)
            for i_c in range(n_c_2):
                out_file.write(b'%.15e %.15e %.15e %.15e %.15e ' %
                               (a2_list[i_c], b2_list[i_c], c2_list[i_c],
                                omega2_list[i_c], tau2_list[i_c]))
            out_file.write(b'\n')

        expmjd = rng.random_sample(100)*200.0
        kp = ParametrizedLightCurveMixin()
        kp.load_parametrized_light_curves(lc_temp_file_name)

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

    def test_applyParametrizedLightCurve_singleExpmjd(self):
        """
        test applyParametrizedLightCurve on a single expmjd value
        by creating a dummy light curve file with known
        parameters, generating magnitudes, and comparing to
        the expected outputs.

        We will use _calc_dflux() to calculate the known truth,
        since that method was tested in test_calc_dflux()
        """

        lc_temp_file_name = tempfile.mktemp(prefix='test_applyParametrizedLightCurve_singleexpmjd',
                                            suffix='.gz')

        rng = np.random.RandomState(5245)
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
            out_file.write(b'# a header\n')
            out_file.write(b'kplr999000000_lc.txt 100 1.0e+02 %d ' % n_c_1)
            for i_c in range(n_c_1):
                out_file.write(b'%e ' % (1.0/(i_c+1)))
            out_file.write(b'%e ' % median1)
            for i_c in range(n_c_1):
                out_file.write(b'%.15e %.15e %.15e %.15e %.15e ' %
                               (a1_list[i_c], b1_list[i_c], c1_list[i_c],
                                omega1_list[i_c], tau1_list[i_c]))
            out_file.write(b'\n')

            out_file.write(b'kplr999000001_lc.txt 100 1.0e+02 %d ' % n_c_2)
            for i_c in range(n_c_2):
                out_file.write(b'%e ' % (1.0/(i_c+1)))
            out_file.write(b'%e ' % median2)
            for i_c in range(n_c_2):
                out_file.write(b'%.15e %.15e %.15e %.15e %.15e ' %
                               (a2_list[i_c], b2_list[i_c], c2_list[i_c],
                                omega2_list[i_c], tau2_list[i_c]))
            out_file.write(b'\n')

        params = {}
        params['lc'] = np.array([999000001, 999000000, None, 999000001])
        params['t0'] = np.array([223.1, 1781.45, None, 32.0])

        kp = ParametrizedLightCurveMixin()
        kp.load_parametrized_light_curves(lc_temp_file_name)

        # first test that passing in an empty set of params
        # results in an empty numpy array (so that the 'dry
        # run' of catalog generation does not fail)
        d_mag_out = kp.applyParametrizedLightCurve([],{},1.0)
        np.testing.assert_array_equal(d_mag_out,
                                      np.array([[],[],[],[],[],[]]))

        expmjd = 59580.0
        d_mag_out = kp.applyParametrizedLightCurve([], params, expmjd)
        self.assertEqual(d_mag_out.shape, (6, 4))

        for i_obj in range(4):
            if i_obj == 2:
                for i_filter in range(6):
                    self.assertEqual(d_mag_out[i_filter][i_obj], 0.0)
            else:
                q_flux, d_flux = kp._calc_dflux(params['lc'][i_obj],
                                                expmjd-params['t0'][i_obj])

                d_mag_truth = -2.5*np.log10(1.0+d_flux/q_flux)
                self.assertFalse(np.isnan(d_mag_truth))
                for i_filter in range(6):
                    self.assertAlmostEqual(d_mag_out[i_filter][i_obj]/d_mag_truth, 1.0, 12)

        sims_clean_up()
        if os.path.exists(lc_temp_file_name):
            os.unlink(lc_temp_file_name)

    def test_applyParametrizedLightCurve_singleExpmjd_as_array(self):
        """
        test applyParametrizedLightCurve on an aray of expmjd values
        that only contains one value by creating a dummy light curve
        file with known parameters, generating magnitudes, and comparing
        to the expected outputs.

        We will use _calc_dflux() to calculate the known truth,
        since that method was tested in test_calc_dflux()
        """

        lc_temp_file_name = tempfile.mktemp(prefix='test_applyParametrizedLightCurve_singleexpmjd',
                                            suffix='.gz')

        rng = np.random.RandomState(5245)
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
            out_file.write(b'# a header\n')
            out_file.write(b'kplr999000000_lc.txt 100 1.0e+02 %d ' % n_c_1)
            for i_c in range(n_c_1):
                out_file.write(b'%e ' % (1.0/(i_c+1)))
            out_file.write(b'%e ' % median1)
            for i_c in range(n_c_1):
                out_file.write(b'%.15e %.15e %.15e %.15e %.15e ' %
                               (a1_list[i_c], b1_list[i_c], c1_list[i_c],
                                omega1_list[i_c], tau1_list[i_c]))
            out_file.write(b'\n')

            out_file.write(b'kplr999000001_lc.txt 100 1.0e+02 %d ' % n_c_2)
            for i_c in range(n_c_2):
                out_file.write(b'%e ' % (1.0/(i_c+1)))
            out_file.write(b'%e ' % median2)
            for i_c in range(n_c_2):
                out_file.write(b'%.15e %.15e %.15e %.15e %.15e ' %
                               (a2_list[i_c], b2_list[i_c], c2_list[i_c],
                                omega2_list[i_c], tau2_list[i_c]))
            out_file.write(b'\n')

        params = {}
        params['lc'] = np.array([999000001, 999000000, None, 999000001])
        params['t0'] = np.array([223.1, 1781.45, None, 32.0])

        kp = ParametrizedLightCurveMixin()
        kp.load_parametrized_light_curves(lc_temp_file_name)

        # first test that passing in an empty set of params
        # results in an empty numpy array (so that the 'dry
        # run' of catalog generation does not fail)
        d_mag_out = kp.applyParametrizedLightCurve([],{},1.0)
        np.testing.assert_array_equal(d_mag_out,
                                      np.array([[],[],[],[],[],[]]))

        expmjd = np.array([59580.0])
        d_mag_out = kp.applyParametrizedLightCurve([], params, expmjd)
        self.assertEqual(d_mag_out.shape, (6, 4, 1))

        for i_obj in range(4):
            if i_obj == 2:
                for i_filter in range(6):
                    self.assertEqual(d_mag_out[i_filter][i_obj], 0.0)
            else:
                q_flux, d_flux = kp._calc_dflux(params['lc'][i_obj],
                                                expmjd-params['t0'][i_obj])

                d_mag_truth = -2.5*np.log10(1.0+d_flux/q_flux)
                self.assertFalse(np.isnan(d_mag_truth))
                for i_filter in range(6):
                    self.assertEqual(len(d_mag_out[i_filter][i_obj]), 1)
                    self.assertAlmostEqual(d_mag_out[i_filter][i_obj][0]/d_mag_truth, 1.0, 12)

        sims_clean_up()
        if os.path.exists(lc_temp_file_name):
            os.unlink(lc_temp_file_name)


    def test_applyParametrizedLightCurve_manyExpmjd(self):
        """
        test applyParametrizedLightCurve on an array of expmjd values
        by creating a dummy light curve file with known
        parameters, generating magnitudes, and comparing to
        the expected outputs.

        We will use _calc_dflux() to calculate the known truth,
        since that method was tested in test_calc_dflux()
        """

        lc_temp_file_name = tempfile.mktemp(prefix='test_applyParametrizedLightCurve_manyexpmjd',
                                            suffix='.gz')

        rng = np.random.RandomState(13291)
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
            out_file.write(b'# a header\n')
            out_file.write(b'kplr999900000_lc.txt 100 1.0e+02 %d ' % n_c_1)
            for i_c in range(n_c_1):
                out_file.write(b'%e ' % (1.0/(i_c+1)))
            out_file.write(b'%e ' % median1)
            for i_c in range(n_c_1):
                out_file.write(b'%.15e %.15e %.15e %.15e %.15e ' %
                               (a1_list[i_c], b1_list[i_c], c1_list[i_c],
                                omega1_list[i_c], tau1_list[i_c]))
            out_file.write(b'\n')

            out_file.write(b'kplr999900001_lc.txt 100 1.0e+02 %d ' % n_c_2)
            for i_c in range(n_c_2):
                out_file.write(b'%e ' % (1.0/(i_c+1)))
            out_file.write(b'%e ' % median2)
            for i_c in range(n_c_2):
                out_file.write(b'%.15e %.15e %.15e %.15e %.15e ' %
                               (a2_list[i_c], b2_list[i_c], c2_list[i_c],
                                omega2_list[i_c], tau2_list[i_c]))
            out_file.write(b'\n')

        params = {}
        params['lc'] = np.array([999900001, 999900000, None, 999900001])
        params['t0'] = np.array([223.1, 1781.45, None, 32.0])

        kp = ParametrizedLightCurveMixin()
        kp.load_parametrized_light_curves(lc_temp_file_name)

        # first test that passing in an empty set of params
        # results in an empty numpy array (so that the 'dry
        # run' of catalog generation does not fail)
        d_mag_out = kp.applyParametrizedLightCurve([],{},1.0)
        np.testing.assert_array_equal(d_mag_out,
                                      np.array([[],[],[],[],[],[]]))

        expmjd = rng.random_sample(10)*10000.0 + 59580.0
        d_mag_out = kp.applyParametrizedLightCurve([], params, expmjd)
        self.assertEqual(d_mag_out.shape, (6, 4, 10))

        for i_obj in range(4):
            if i_obj == 2:
                for i_filter in range(6):
                    np.testing.assert_array_equal(d_mag_out[i_filter][i_obj],
                                                  np.zeros(10))
            else:
                q_flux, d_flux = kp._calc_dflux(params['lc'][i_obj],
                                                expmjd-params['t0'][i_obj])

                d_mag_truth = -2.5*np.log10(1.0+d_flux/q_flux)
                nan_vals = np.where(np.isnan(d_mag_truth))
                self.assertEqual(len(nan_vals[0]), 0)
                for i_filter in range(6):
                    np.testing.assert_array_equal(d_mag_out[i_filter][i_obj], d_mag_truth)

        sims_clean_up()
        if os.path.exists(lc_temp_file_name):
            os.unlink(lc_temp_file_name)

    def test_ParametrizedLightCurve_in_catalog(self):
        """
        Test the performance of applyParametrizedLightCurve()
        in the context of an InstanceCatalog
        """

        # Create dummy light curve parameters
        lc_temp_file_name = tempfile.mktemp(prefix='test_ParametrizedLightCurve_in_catalog',
                                            suffix='.gz')

        rng = np.random.RandomState(1621145)
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
            out_file.write(b'# a header\n')
            out_file.write(b'kplr999990000_lc.txt 100 1.0e+02 %d ' % n_c_1)
            for i_c in range(n_c_1):
                out_file.write(b'%e ' % (1.0/(i_c+1)))
            out_file.write(b'%e ' % median1)
            for i_c in range(n_c_1):
                out_file.write(b'%.15e %.15e %.15e %.15e %.15e ' %
                               (a1_list[i_c], b1_list[i_c], c1_list[i_c],
                                omega1_list[i_c], tau1_list[i_c]))
            out_file.write(b'\n')

            out_file.write(b'kplr999990001_lc.txt 100 1.0e+02 %d ' % n_c_2)
            for i_c in range(n_c_2):
                out_file.write(b'%e ' % (1.0/(i_c+1)))
            out_file.write(b'%e ' % median2)
            for i_c in range(n_c_2):
                out_file.write(b'%.15e %.15e %.15e %.15e %.15e ' %
                               (a2_list[i_c], b2_list[i_c], c2_list[i_c],
                                omega2_list[i_c], tau2_list[i_c]))
            out_file.write(b'\n')

        # Create dummy database of astrophysical sources
        db_temp_file_name = tempfile.mktemp(prefix='test_ParametrizedLightCurve_in_catalog_db',
                                            suffix='.txt')

        lc_list = [999990001, None, 999990001, 999990000]
        t0_list = [1729.1, None, 2345.1, 10.9]

        with open(db_temp_file_name, 'w') as out_file:
            out_file.write('# a header\n')
            for i_obj in range(len(lc_list)):
                if lc_list[i_obj] is not None:
                    paramStr = '{"m":"kplr", "p":{"lc":%d, "t0":%.3f}}' % (lc_list[i_obj], t0_list[i_obj])
                else:
                    paramStr = None
                out_file.write('%d;10.0;20.0;0.01;0.01;%s\n' % (i_obj, paramStr))

        dtype = np.dtype([('simobjid', int), ('ra', float), ('dec', float),
                          ('ebv', float), ('parallax', float), ('varParamStr', str, 100)])
        db = fileDBObject(db_temp_file_name, runtable='test', dtype=dtype, delimiter=';',
                          idColKey='simobjid')

        class ParametrizedVarParamStrCat(InstanceCatalog, VariabilityStars):
            column_outputs = ['simobjid', 'delta_lsst_u', 'delta_lsst_g', 'delta_lsst_r',
                              'delta_lsst_i', 'delta_lsst_z', 'delta_lsst_y']
            default_formats = {'f':'%.15g'}


        obs = ObservationMetaData(mjd=59580.0)
        cat = ParametrizedVarParamStrCat(db, obs_metadata=obs)
        cat.load_parametrized_light_curves(lc_temp_file_name)
        cat_out_name = tempfile.mktemp(prefix='test_ParametrizedLightCurve_in_cat_out',
                                       suffix='.txt')

        cat.write_catalog(cat_out_name)

        kp = ParametrizedLightCurveMixin()
        cat_dtype = np.dtype([('simobjid', int), ('du', float), ('dg', float),
                              ('dr', float), ('di', float), ('dz', float),
                              ('dy', float)])

        cat_data = np.genfromtxt(cat_out_name, dtype=cat_dtype, delimiter=', ')

        for i_obj in range(len(cat_data)):
            obj_id = cat_data['simobjid'][i_obj]
            if lc_list[obj_id] is None:
                self.assertEqual(cat_data['du'][i_obj], 0.0)
                self.assertEqual(cat_data['dg'][i_obj], 0.0)
                self.assertEqual(cat_data['dr'][i_obj], 0.0)
                self.assertEqual(cat_data['di'][i_obj], 0.0)
                self.assertEqual(cat_data['dz'][i_obj], 0.0)
                self.assertEqual(cat_data['dy'][i_obj], 0.0)
            else:
                q_flux, d_flux = kp._calc_dflux(lc_list[obj_id], obs.mjd.TAI-t0_list[obj_id])
                d_mag_true = -2.5*np.log10(1.0+d_flux/q_flux)
                self.assertGreater(np.abs(d_mag_true), 0.0001)
                self.assertAlmostEqual(cat_data['du'][i_obj], d_mag_true, 15)
                self.assertAlmostEqual(cat_data['dg'][i_obj], d_mag_true, 15)
                self.assertAlmostEqual(cat_data['dr'][i_obj], d_mag_true, 15)
                self.assertAlmostEqual(cat_data['di'][i_obj], d_mag_true, 15)
                self.assertAlmostEqual(cat_data['dz'][i_obj], d_mag_true, 15)
                self.assertAlmostEqual(cat_data['dy'][i_obj], d_mag_true, 15)

        if os.path.exists(cat_out_name):
            os.unlink(cat_out_name)

        if os.path.exists(db_temp_file_name):
            os.unlink(db_temp_file_name)

        sims_clean_up()
        if os.path.exists(lc_temp_file_name):
            os.unlink(lc_temp_file_name)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
