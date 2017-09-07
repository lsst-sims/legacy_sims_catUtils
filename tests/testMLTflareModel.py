from __future__ import with_statement
import unittest
import os
import numpy as np
import sqlite3
import json
import tempfile
import shutil
from astropy.analytic_functions import blackbody_lambda
import lsst.utils.tests
from lsst.utils import getPackageDir
from lsst.utils.tests import getTempFilePath
from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.catUtils.mixins import VariabilityStars
from lsst.sims.catUtils.mixins import MLTflaringMixin
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catUtils.mixins import PhotometryStars
from lsst.sims.utils import ObservationMetaData
from lsst.sims.photUtils import SedList, BandpassDict, Sed
from lsst.sims.utils import radiansFromArcsec
from lsst.sims.photUtils import PhotometricParameters
from lsst.sims.utils.CodeUtilities import sims_clean_up

ROOT = os.path.abspath(os.path.dirname(__file__))


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
               ('parallax', 'parallax*PI()/648000000.'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', str, 40),
               ('ebv', 'galacticAv/3.1')]


class QuiescentCatalog(PhotometryStars, InstanceCatalog):
    column_outputs = ['id', 'lsst_u', 'lsst_g']

class FlaringCatalog(PhotometryStars, VariabilityStars,
                     MLTflaringMixin, InstanceCatalog):
    column_outputs = ['id', 'lsst_u', 'lsst_g']


class MLT_flare_test_case(unittest.TestCase):

    longMessage = True

    @classmethod
    def setUpClass(cls):
        cls.scratch_dir = tempfile.mkdtemp(dir=ROOT, prefix="MLT_flare_test_case-")

        # Create a dummy mlt light curve file
        cls.mlt_lc_name = os.path.join(cls.scratch_dir,
                                       'test_mlt_lc_file.npz')

        lc_files = {}
        amp = 1.0e36
        lc_files['lc_1_time'] = np.arange(0.0, 3652.51, 0.1)
        lc_files['lc_1_u'] = amp*(1.0+np.power(np.sin(lc_files['lc_1_time']/100.0), 2))
        lc_files['lc_1_g'] = amp*(1.0+np.power(np.cos(lc_files['lc_1_time']/100.0), 2))

        amp = 2.0e35
        lc_files['lc_2_time'] = np.arange(0.0, 365.251, 0.01)
        lc_files['lc_2_u'] = amp*(1.0+np.power(np.sin(lc_files['lc_2_time']/50.0), 2))
        lc_files['lc_2_g'] = amp*(1.0+np.power(np.cos(lc_files['lc_2_time']/50.0), 2))

        with open(cls.mlt_lc_name, 'wb') as file_handle:
            np.savez(file_handle, **lc_files)

        # Create a database of stars using these light curves
        cls.db_name = os.path.join(cls.scratch_dir, 'test_mlt_db.db')

        conn = sqlite3.connect(cls.db_name)
        cursor = conn.cursor()
        cursor.execute('''CREATE TABLE mlt_test
                       (simobjid int, ra real, decl real, sedfilename text,
                        varParamStr text, parallax real, galacticAv real,
                        magNorm real)''')
        conn.commit()

        cursor.execute('''INSERT INTO mlt_test VALUES( 0, 25.0, 31.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz',
                       '{"m": "MLT",
                         "p": {"lc": "lc_1.txt", "t0": 456.2}}',
                       0.25, 2.432, 17.1)''')

        cursor.execute('''INSERT INTO mlt_test VALUES( 1, 25.2, 32.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz',
                       '{"m": "MLT",
                         "p": {"lc": "lc_1.txt", "t0": 41006.2}}',
                       0.15, 1.876, 17.2)''')

        cursor.execute('''INSERT INTO mlt_test VALUES( 2, 25.3, 10.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz',
                       '{"m": "MLT",
                         "p": {"lc": "lc_2.txt", "t0": 117.2}}',
                       0.3, 2.654, 17.3)''')

        cursor.execute('''INSERT INTO mlt_test VALUES( 3, 25.4, 11.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz',
                       '{"m": "MLT",
                         "p": {"lc": "lc_2.txt", "t0": 10456.2}}',
                       0.22, 2.364, 17.4)''')
        conn.commit()
        conn.close()


    @classmethod
    def tearDownClass(cls):
        sims_clean_up()
        if os.path.exists(cls.mlt_lc_name):
            os.unlink(cls.mlt_lc_name)

        if os.path.exists(cls.db_name):
            os.unlink(cls.db_name)

        if os.path.exists(cls.scratch_dir):
            shutil.rmtree(cls.scratch_dir)

    def test_flare_lc_failure(self):
        """
        Test that the correct exception is thrown when you try
        to interpolate from an MLT light curve cache that does not
        exist
        """
        with getTempFilePath('.txt') as cat_name:
            db = MLT_test_DB(database=self.db_name, driver='sqlite')
            obs = ObservationMetaData(mjd=67432.1)
            flare_cat = FlaringCatalog(db, obs_metadata=obs)
            flare_cat._mlt_lc_file = 'a_nonexistent_cache'
            with self.assertRaises(RuntimeError) as context:
                flare_cat.write_catalog(cat_name)
            self.assertIn('get_mdwarf_flares.sh',
                          context.exception.args[0])

    def test_flare_magnitudes(self):
        """
        Test that we get the expected magnitudes out
        """
        db = MLT_test_DB(database=self.db_name, driver='sqlite')

        # load the quiescent SEDs of the objects in our catalog
        sed_list = SedList(['lte028-5.0+0.5a+0.0.BT-Settl.spec.gz']*4,
                           [17.1, 17.2, 17.3, 17.4],
                           galacticAvList = [2.432, 1.876, 2.654, 2.364])

        bp_dict = BandpassDict.loadTotalBandpassesFromFiles()

        # calculate the quiescent fluxes of the objects in our catalog
        baseline_fluxes = bp_dict.fluxListForSedList(sed_list)

        bb_wavelen = np.arange(100.0, 1600.0, 0.1)
        bb_flambda = blackbody_lambda(bb_wavelen*10.0, 9000.0)

        # this data is taken from the setUpClass() classmethod above
        t0_list = [456.2, 41006.2, 117.2, 10456.2]
        av_list = [2.432, 1.876, 2.654, 2.364]
        parallax_list = np.array([0.25, 0.15, 0.3, 0.22])
        distance_list = 1.0/(206265.0*radiansFromArcsec(0.001*parallax_list))
        distance_list *= 3.0857e18  # convert to cm

        dtype = np.dtype([('id', int), ('u', float), ('g', float)])

        photParams = PhotometricParameters()

        ss = Sed()

        quiet_cat_name = os.path.join(self.scratch_dir, 'mlt_quiet_cat.txt')
        flare_cat_name = os.path.join(self.scratch_dir, 'mlt_flaring_cat.txt')

        # loop over several MJDs and verify that, to within a
        # milli-mag, our flaring model gives us the magnitudes
        # expected, given the light curves specified in
        # setUpClass()
        for mjd in (59580.0, 60000.0, 70000.0, 80000.0):

            obs = ObservationMetaData(mjd=mjd)

            quiet_cat = QuiescentCatalog(db, obs_metadata=obs)
            quiet_cat.write_catalog(quiet_cat_name)

            flare_cat = FlaringCatalog(db, obs_metadata=obs)
            flare_cat._mlt_lc_file = self.mlt_lc_name
            flare_cat.write_catalog(flare_cat_name)

            quiescent_data = np.genfromtxt(quiet_cat_name, dtype=dtype, delimiter=',')
            flaring_data = np.genfromtxt(flare_cat_name, dtype=dtype, delimiter=',')

            for ix in range(len(flaring_data)):
                obj_id = flaring_data['id'][ix]
                self.assertEqual(obj_id, ix)

                # the models below are as specified in the
                # setUpClass() method
                if obj_id == 0 or obj_id == 1:
                    amp = 1.0e36
                    dt = 3652.5
                    t_min = flare_cat._survey_start - t0_list[obj_id]

                    tt = mjd - t_min
                    while tt > dt:
                        tt -= dt

                    u_flux = amp*(1.0+np.power(np.sin(tt/100.0), 2))
                    g_flux = amp*(1.0+np.power(np.cos(tt/100.0), 2))
                else:
                    amp = 2.0e35
                    dt = 365.25
                    t_min = flare_cat._survey_start - t0_list[obj_id]

                    tt = mjd - t_min
                    while tt > dt:
                        tt -= dt
                    u_flux = amp*(1.0+np.power(np.sin(tt/50.0), 2))
                    g_flux = amp*(1.0+np.power(np.cos(tt/50.0), 2))

                # calculate the multiplicative effect of dust on a 9000K
                # black body
                bb_sed = Sed(wavelen=bb_wavelen, flambda=bb_flambda)
                u_bb_flux = bb_sed.calcFlux(bp_dict['u'])
                g_bb_flux = bb_sed.calcFlux(bp_dict['g'])
                a_x, b_x = bb_sed.setupCCMab()
                bb_sed.addCCMDust(a_x, b_x, A_v=av_list[obj_id])
                u_bb_dusty_flux = bb_sed.calcFlux(bp_dict['u'])
                g_bb_dusty_flux = bb_sed.calcFlux(bp_dict['g'])

                dust_u = u_bb_dusty_flux/u_bb_flux
                dust_g = g_bb_dusty_flux/g_bb_flux

                area = 4.0*np.pi*np.power(distance_list[obj_id], 2)
                tot_u_flux = baseline_fluxes[obj_id][0] + u_flux*dust_u*photParams.effarea/area
                tot_g_flux = baseline_fluxes[obj_id][1] + g_flux*dust_g*photParams.effarea/area

                msg = ('failed on object %d; mjd %.2f\n u_quiet %e u_flare %e\n g_quiet %e g_flare %e' %
                       (obj_id, mjd, quiescent_data['u'][obj_id], flaring_data['u'][obj_id],
                        quiescent_data['g'][obj_id], flaring_data['g'][obj_id]))

                self.assertEqual(quiescent_data['id'][obj_id], flaring_data['id'][obj_id], msg=msg)
                self.assertAlmostEqual(ss.magFromFlux(baseline_fluxes[obj_id][0]),
                                       quiescent_data['u'][obj_id], 3, msg=msg)
                self.assertAlmostEqual(ss.magFromFlux(baseline_fluxes[obj_id][1]),
                                       quiescent_data['g'][obj_id], 3, msg=msg)
                self.assertAlmostEqual(ss.magFromFlux(tot_u_flux), flaring_data['u'][obj_id],
                                       3, msg=msg)
                self.assertAlmostEqual(ss.magFromFlux(tot_g_flux), flaring_data['g'][obj_id],
                                       3, msg=msg)
                self.assertGreater(np.abs(flaring_data['g'][obj_id]-quiescent_data['g'][obj_id]),
                                   0.001, msg=msg)
                self.assertGreater(np.abs(flaring_data['u'][obj_id]-quiescent_data['u'][obj_id]),
                                   0.001, msg=msg)

        if os.path.exists(quiet_cat_name):
            os.unlink(quiet_cat_name)
        if os.path.exists(flare_cat_name):
            os.unlink(flare_cat_name)

    def test_MLT_many_mjd(self):
        """
        This test will verify that applyMLTflaring responds properly when given
        an array/vector of MJD values.
        """
        db = MLT_test_DB(database=self.db_name, driver='sqlite')
        rng = np.random.RandomState(16)
        mjd_arr = rng.random_sample(17)*3653.3+59580.0
        objid = []
        parallax = []
        ebv = []
        quiescent_u = []
        quiescent_g = []
        delta_u = []
        delta_g = []
        varparams = []
        for mjd in mjd_arr:
            obs = ObservationMetaData(mjd=mjd)
            cat = FlaringCatalog(db, obs_metadata=obs,
                                 column_outputs=['parallax', 'ebv',
                                                 'quiescent_lsst_u',
                                                 'quiescent_lsst_g',
                                                 'varParamStr',
                                                 'delta_lsst_u',
                                                 'delta_lsst_g'])

            cat._mlt_lc_file = self.mlt_lc_name
            n_obj = 0
            for line in cat.iter_catalog():
                n_obj += 1
                objid.append(line[0])
                parallax.append(line[3])
                ebv.append(line[4])
                quiescent_u.append(line[5])
                quiescent_g.append(line[6])
                varparams.append(line[7])
                delta_u.append(line[8])
                delta_g.append(line[9])

        objid = np.array(objid)
        parallax = np.array(parallax)
        ebv = np.array(ebv)
        quiescent_u = np.array(quiescent_u)
        quiescent_g = np.array(quiescent_g)
        delta_u = np.array(delta_u)
        delta_g = np.array(delta_g)

        self.assertEqual(len(parallax), n_obj*len(mjd_arr))
        np.testing.assert_array_equal(objid, np.array([0,1,2,3]*len(mjd_arr)))

        quiescent_mags = {}
        quiescent_mags['u'] = quiescent_u
        quiescent_mags['g'] = quiescent_g

        params = {}
        params['lc'] = []
        params['t0'] = []
        for ix in range(4):
            local_dict = json.loads(varparams[ix])
            params['lc'].append(local_dict['p']['lc'])
            params['t0'].append(local_dict['p']['t0'])
        params['lc'] = np.array(params['lc'])
        params['t0'] = np.array(params['t0'])

        mlt_obj = MLTflaringMixin()
        mlt_obj.photParams = PhotometricParameters()
        mlt_obj.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()
        mlt_obj._mlt_lc_file = cat._mlt_lc_file
        mlt_obj._actually_calculated_columns = ['delta_lsst_u', 'delta_lsst_g']
        valid_dex = [np.arange(4, dtype=int)]
        delta_mag_vector = mlt_obj.applyMLTflaring(valid_dex, params, mjd_arr,
                                                   parallax=parallax,
                                                   ebv=ebv,
                                                   quiescent_mags=quiescent_mags)
        n_time = len(mjd_arr)
        n_obj = 4
        self.assertEqual(delta_mag_vector.shape, (6, n_obj, n_time))

        for i_time, mjd in enumerate(mjd_arr):
            for i_obj in range(n_obj):
                for i_band in range(6):
                    if i_band==0:
                        self.assertEqual(delta_mag_vector[i_band][i_obj][i_time],
                                         delta_u[i_time*n_obj + i_obj])
                    elif i_band==1:
                        self.assertEqual(delta_mag_vector[i_band][i_obj][i_time],
                                         delta_g[i_time*n_obj + i_obj])
                    else:
                        self.assertEqual(delta_mag_vector[i_band][i_obj][i_time], 0.0)

    def test_MLT_many_mjd_some_invalid(self):
        """
        This test will verify that applyMLTflaring responds properly when given
        an array/vector of MJD values in the case where some stars are not
        marked as valid MLT stars.
        """
        db = MLT_test_DB(database=self.db_name, driver='sqlite')
        rng = np.random.RandomState(16)
        mjd_arr = rng.random_sample(17)*3653.3+59580.0
        objid = []
        parallax = []
        ebv = []
        quiescent_u = []
        quiescent_g = []
        delta_u = []
        delta_g = []
        varparams = []
        for mjd in mjd_arr:
            obs = ObservationMetaData(mjd=mjd)
            cat = FlaringCatalog(db, obs_metadata=obs,
                                 column_outputs=['parallax', 'ebv',
                                                 'quiescent_lsst_u',
                                                 'quiescent_lsst_g',
                                                 'varParamStr',
                                                 'delta_lsst_u',
                                                 'delta_lsst_g'])

            cat._mlt_lc_file = self.mlt_lc_name
            n_obj = 0
            for line in cat.iter_catalog():
                n_obj += 1
                objid.append(line[0])
                parallax.append(line[3])
                ebv.append(line[4])
                quiescent_u.append(line[5])
                quiescent_g.append(line[6])
                varparams.append(line[7])
                delta_u.append(line[8])
                delta_g.append(line[9])

        objid = np.array(objid)
        parallax = np.array(parallax)
        ebv = np.array(ebv)
        quiescent_u = np.array(quiescent_u)
        quiescent_g = np.array(quiescent_g)
        delta_u = np.array(delta_u)
        delta_g = np.array(delta_g)

        self.assertEqual(len(parallax), n_obj*len(mjd_arr))
        np.testing.assert_array_equal(objid, np.array([0,1,2,3]*len(mjd_arr)))

        quiescent_mags = {}
        quiescent_mags['u'] = quiescent_u
        quiescent_mags['g'] = quiescent_g

        params = {}
        params['lc'] = []
        params['t0'] = []
        for ix in range(n_obj):
            local_dict = json.loads(varparams[ix])
            params['lc'].append(local_dict['p']['lc'])
            params['t0'].append(local_dict['p']['t0'])
        params['lc'] = np.array(params['lc'])
        params['t0'] = np.array(params['t0'])

        mlt_obj = MLTflaringMixin()
        mlt_obj.photParams = PhotometricParameters()
        mlt_obj.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()
        mlt_obj._mlt_lc_file = cat._mlt_lc_file
        mlt_obj._actually_calculated_columns = ['delta_lsst_u', 'delta_lsst_g']
        valid_dex = [] # applyMLT does not actually use valid_dex; it looks for non-None params['lc']
        for ix in range(n_obj):
            if ix not in (1,2):
                params['lc'][ix] = None
        delta_mag_vector = mlt_obj.applyMLTflaring(valid_dex, params, mjd_arr,
                                                   parallax=parallax,
                                                   ebv=ebv,
                                                   quiescent_mags=quiescent_mags)
        n_time = len(mjd_arr)
        self.assertEqual(delta_mag_vector.shape, (6, n_obj, n_time))

        for i_time, mjd in enumerate(mjd_arr):
            for i_obj in range(n_obj):
                if i_obj not in (1,2):
                    for i_band in range(6):
                        self.assertEqual(delta_mag_vector[i_band][i_obj][i_time], 0.0)
                    continue

                for i_band in range(6):
                    if i_band==0:
                        self.assertEqual(delta_mag_vector[i_band][i_obj][i_time],
                                         delta_u[i_time*n_obj + i_obj])
                    elif i_band==1:
                        self.assertEqual(delta_mag_vector[i_band][i_obj][i_time],
                                         delta_g[i_time*n_obj + i_obj])
                    else:
                        self.assertEqual(delta_mag_vector[i_band][i_obj][i_time], 0.0)

    def test_mlt_clean_up(self):
        """
        Test that the MLT cache is correctly loaded after sims_clean_up is
        called.
        """
        db = MLT_test_DB(database=self.db_name, driver='sqlite')
        obs = ObservationMetaData(mjd=60000.0)
        cat = FlaringCatalog(db, obs_metadata=obs)
        cat._mlt_lc_file = self.mlt_lc_name
        cat_name_1 = os.path.join(self.scratch_dir,'mlt_clean_test_cat_1.txt')
        cat.write_catalog(cat_name_1)
        sims_clean_up()

        # re-generate the same catalog and verify that its
        # contents are unchanged
        db = MLT_test_DB(database=self.db_name, driver='sqlite')
        obs = ObservationMetaData(mjd=60000.0)
        cat = FlaringCatalog(db, obs_metadata=obs)
        cat._mlt_lc_file = self.mlt_lc_name
        cat_name_2 = os.path.join(self.scratch_dir,'mlt_clean_test_cat_2.txt')
        cat.write_catalog(cat_name_2)
        with open(cat_name_1, 'r') as in_file_1:
            lines_1 = in_file_1.readlines()
        with open(cat_name_2, 'r') as in_file_2:
            lines_2 = in_file_2.readlines()
        self.assertGreater(len(lines_1), 1)
        self.assertEqual(len(lines_1), len(lines_2))
        for line in lines_1:
            self.assertIn(line, lines_2)

        if os.path.exists(cat_name_1):
            os.unlink(cat_name_1)
        if os.path.exists(cat_name_2):
            os.unlink(cat_name_2)


class MLT_flare_mixed_with_none_model_test_case(unittest.TestCase):
    """
    This test class duplicates MLT_flare_model_test_case, except
    that some of the objects in the database have 'None' for their
    varParamStr.
    """

    longMessage = True

    @classmethod
    def setUpClass(cls):
        cls.scratch_dir = tempfile.mkdtemp(dir=ROOT, prefix='MLT_flare_mixed_with_none_model_test_case-')

        # Create a dummy mlt light curve file
        cls.mlt_lc_name = os.path.join(cls.scratch_dir,
                                       'test_mlt_mixed_with_none_lc_file.npz')

        lc_files = {}
        amp = 1.0e36
        lc_files['lc_1_time'] = np.arange(0.0, 3652.51, 0.1)
        lc_files['lc_1_u'] = amp*(1.0+np.power(np.sin(lc_files['lc_1_time']/100.0), 2))
        lc_files['lc_1_g'] = amp*(1.0+np.power(np.cos(lc_files['lc_1_time']/100.0), 2))

        amp = 2.0e35
        lc_files['lc_2_time'] = np.arange(0.0, 365.251, 0.01)
        lc_files['lc_2_u'] = amp*(1.0+np.power(np.sin(lc_files['lc_2_time']/50.0), 2))
        lc_files['lc_2_g'] = amp*(1.0+np.power(np.cos(lc_files['lc_2_time']/50.0), 2))

        with open(cls.mlt_lc_name, 'wb') as file_handle:
            np.savez(file_handle, **lc_files)

        # Create a database of stars using these light curves
        cls.db_name = os.path.join(cls.scratch_dir, 'test_mlt_mixed_with_none_db.db')

        conn = sqlite3.connect(cls.db_name)
        cursor = conn.cursor()
        cursor.execute('''CREATE TABLE mlt_test
                       (simobjid int, ra real, decl real, sedfilename text,
                        varParamStr text, parallax real, galacticAv real,
                        magNorm real)''')
        conn.commit()

        cursor.execute('''INSERT INTO mlt_test VALUES( 0, 25.0, 31.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz',
                       '{"m": "MLT",
                         "p": {"lc": "lc_1.txt", "t0": 456.2}}',
                       0.25, 2.432, 17.1)''')

        cursor.execute('''INSERT INTO mlt_test VALUES( 1, 25.2, 32.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz',
                       '{"m": "MLT",
                         "p": {"lc": "lc_1.txt", "t0": 41006.2}}',
                       0.15, 1.876, 17.2)''')

        cursor.execute('''INSERT INTO mlt_test VALUES( 2, 25.3, 10.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz',
                       '{"m": "MLT",
                         "p": {"lc": "lc_2.txt", "t0": 117.2}}',
                       0.3, 2.654, 17.3)''')

        cursor.execute('''INSERT INTO mlt_test VALUES( 3, 25.4, 11.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz', NULL,
                       0.22, 2.364, 17.4)''')
        conn.commit()
        conn.close()


    @classmethod
    def tearDownClass(cls):
        sims_clean_up()
        if os.path.exists(cls.mlt_lc_name):
            os.unlink(cls.mlt_lc_name)

        if os.path.exists(cls.db_name):
            os.unlink(cls.db_name)

        if os.path.exists(cls.scratch_dir):
            shutil.rmtree(cls.scratch_dir)

    def test_flare_magnitudes_mixed_with_none(self):
        """
        Test that we get the expected magnitudes out
        """
        db = MLT_test_DB(database=self.db_name, driver='sqlite')

        # load the quiescent SEDs of the objects in our catalog
        sed_list = SedList(['lte028-5.0+0.5a+0.0.BT-Settl.spec.gz']*4,
                           [17.1, 17.2, 17.3, 17.4],
                           galacticAvList = [2.432, 1.876, 2.654, 2.364])

        bp_dict = BandpassDict.loadTotalBandpassesFromFiles()

        # calculate the quiescent fluxes of the objects in our catalog
        baseline_fluxes = bp_dict.fluxListForSedList(sed_list)

        bb_wavelen = np.arange(100.0, 1600.0, 0.1)
        bb_flambda = blackbody_lambda(bb_wavelen*10.0, 9000.0)

        # this data is taken from the setUpClass() classmethod above
        t0_list = [456.2, 41006.2, 117.2, 10456.2]
        av_list = [2.432, 1.876, 2.654, 2.364]
        parallax_list = np.array([0.25, 0.15, 0.3, 0.22])
        distance_list = 1.0/(206265.0*radiansFromArcsec(0.001*parallax_list))
        distance_list *= 3.0857e18  # convert to cm

        dtype = np.dtype([('id', int), ('u', float), ('g', float)])

        photParams = PhotometricParameters()

        ss = Sed()

        quiet_cat_name = os.path.join(self.scratch_dir, 'mlt_mixed_with_none_quiet_cat.txt')
        flare_cat_name = os.path.join(self.scratch_dir, 'mlt_mixed_with_none_flaring_cat.txt')

        # loop over several MJDs and verify that, to within a
        # milli-mag, our flaring model gives us the magnitudes
        # expected, given the light curves specified in
        # setUpClass()
        for mjd in (59580.0, 60000.0, 70000.0, 80000.0):

            obs = ObservationMetaData(mjd=mjd)

            quiet_cat = QuiescentCatalog(db, obs_metadata=obs)
            quiet_cat.write_catalog(quiet_cat_name)

            flare_cat = FlaringCatalog(db, obs_metadata=obs)
            flare_cat._mlt_lc_file = self.mlt_lc_name
            flare_cat.write_catalog(flare_cat_name)

            quiescent_data = np.genfromtxt(quiet_cat_name, dtype=dtype, delimiter=',')
            flaring_data = np.genfromtxt(flare_cat_name, dtype=dtype, delimiter=',')

            for ix in range(len(flaring_data)):
                obj_id = flaring_data['id'][ix]
                self.assertEqual(obj_id, ix)

                # the models below are as specified in the
                # setUpClass() method
                if obj_id == 0 or obj_id == 1:
                    amp = 1.0e36
                    dt = 3652.5
                    t_min = flare_cat._survey_start - t0_list[obj_id]

                    tt = mjd - t_min
                    while tt > dt:
                        tt -= dt

                    u_flux = amp*(1.0+np.power(np.sin(tt/100.0), 2))
                    g_flux = amp*(1.0+np.power(np.cos(tt/100.0), 2))
                elif obj_id==2:
                    amp = 2.0e35
                    dt = 365.25
                    t_min = flare_cat._survey_start - t0_list[obj_id]

                    tt = mjd - t_min
                    while tt > dt:
                        tt -= dt
                    u_flux = amp*(1.0+np.power(np.sin(tt/50.0), 2))
                    g_flux = amp*(1.0+np.power(np.cos(tt/50.0), 2))
                else:
                    u_flux = 0.0
                    g_flux = 0.0

                # calculate the multiplicative effect of dust on a 9000K
                # black body
                bb_sed = Sed(wavelen=bb_wavelen, flambda=bb_flambda)
                u_bb_flux = bb_sed.calcFlux(bp_dict['u'])
                g_bb_flux = bb_sed.calcFlux(bp_dict['g'])
                a_x, b_x = bb_sed.setupCCMab()
                bb_sed.addCCMDust(a_x, b_x, A_v=av_list[obj_id])
                u_bb_dusty_flux = bb_sed.calcFlux(bp_dict['u'])
                g_bb_dusty_flux = bb_sed.calcFlux(bp_dict['g'])

                dust_u = u_bb_dusty_flux/u_bb_flux
                dust_g = g_bb_dusty_flux/g_bb_flux

                area = 4.0*np.pi*np.power(distance_list[obj_id], 2)
                tot_u_flux = baseline_fluxes[obj_id][0] + u_flux*dust_u*photParams.effarea/area
                tot_g_flux = baseline_fluxes[obj_id][1] + g_flux*dust_g*photParams.effarea/area

                msg = ('failed on object %d; mjd %.2f\n u_quiet %e u_flare %e\n g_quiet %e g_flare %e' %
                       (obj_id, mjd, quiescent_data['u'][obj_id], flaring_data['u'][obj_id],
                        quiescent_data['g'][obj_id], flaring_data['g'][obj_id]))

                self.assertEqual(quiescent_data['id'][obj_id], flaring_data['id'][obj_id], msg=msg)
                self.assertAlmostEqual(ss.magFromFlux(baseline_fluxes[obj_id][0]),
                                       quiescent_data['u'][obj_id], 3, msg=msg)
                self.assertAlmostEqual(ss.magFromFlux(baseline_fluxes[obj_id][1]),
                                       quiescent_data['g'][obj_id], 3, msg=msg)
                self.assertAlmostEqual(ss.magFromFlux(tot_u_flux), flaring_data['u'][obj_id],
                                       3, msg=msg)
                self.assertAlmostEqual(ss.magFromFlux(tot_g_flux), flaring_data['g'][obj_id],
                                       3, msg=msg)
                if obj_id != 3:
                    self.assertGreater(np.abs(flaring_data['g'][obj_id]-quiescent_data['g'][obj_id]),
                                       0.001, msg=msg)
                    self.assertGreater(np.abs(flaring_data['u'][obj_id]-quiescent_data['u'][obj_id]),
                                       0.001, msg=msg)
                else:
                    self.assertEqual(flaring_data['g'][obj_id]-quiescent_data['g'][obj_id], 0.0, msg=msg)
                    self.assertEqual(flaring_data['u'][obj_id]-quiescent_data['u'][obj_id], 0.0, msg=msg)

        if os.path.exists(quiet_cat_name):
            os.unlink(quiet_cat_name)
        if os.path.exists(flare_cat_name):
            os.unlink(flare_cat_name)


from lsst.sims.catUtils.mixins import Variability
from lsst.sims.catalogs.decorators import register_method


class DummyVariabilityMixin(Variability):

    scratch_dir = None

    @register_method('dummy')
    def applyDummy(self, valid_dexes, params, expmjd):
        if len(params) == 0:
            return np.array([[],[],[],[],[],[]])
        dtype = np.dtype([('t', float), ('du', float), ('dg', float),
                          ('dr', float), ('di', float), ('dz', float),
                          ('dy', float)])
        dmag = np.zeros((6, self.num_variable_obj(params)))
        for ix in valid_dexes[0]:
            lc_name = os.path.join(self.scratch_dir, params['lc'][ix])
            lc_data = np.genfromtxt(lc_name, dtype=dtype)
            dmag[0][ix] = np.interp(expmjd, lc_data['t'], lc_data['du'])
            dmag[1][ix] = np.interp(expmjd, lc_data['t'], lc_data['dg'])
            dmag[2][ix] = np.interp(expmjd, lc_data['t'], lc_data['dr'])
            dmag[3][ix] = np.interp(expmjd, lc_data['t'], lc_data['di'])
            dmag[4][ix] = np.interp(expmjd, lc_data['t'], lc_data['dz'])
            dmag[5][ix] = np.interp(expmjd, lc_data['t'], lc_data['dy'])
        return dmag

class FlaringCatalogDummy(PhotometryStars, VariabilityStars,
                          MLTflaringMixin, DummyVariabilityMixin,
                          InstanceCatalog):
    column_outputs = ['id', 'lsst_u', 'lsst_g']


class MLT_flare_mixed_with_dummy_model_test_case(unittest.TestCase):
    """
    This test class duplicates MLT_flare_model_test_case, except
    that one of the objects in the database has a different
    variability method applied to it to make sure that
    applyVariability properly handles cases where two variability
    models have identically named params ('lc' in this case).
    """

    longMessage = True

    @classmethod
    def setUpClass(cls):
        cls.scratch_dir = tempfile.mkdtemp(dir=ROOT, prefix='MLT_flare_mixed_with_dummy_model_test_case-')

        # Create a dummy mlt light curve file
        cls.mlt_lc_name = os.path.join(cls.scratch_dir,
                                       'test_mlt_mixed_with_dummy_lc_file.npz')

        lc_files = {}
        amp = 1.0e36
        lc_files['lc_1_time'] = np.arange(0.0, 3652.51, 0.1)
        lc_files['lc_1_u'] = amp*(1.0+np.power(np.sin(lc_files['lc_1_time']/100.0), 2))
        lc_files['lc_1_g'] = amp*(1.0+np.power(np.cos(lc_files['lc_1_time']/100.0), 2))

        amp = 2.0e35
        lc_files['lc_2_time'] = np.arange(0.0, 365.251, 0.01)
        lc_files['lc_2_u'] = amp*(1.0+np.power(np.sin(lc_files['lc_2_time']/50.0), 2))
        lc_files['lc_2_g'] = amp*(1.0+np.power(np.cos(lc_files['lc_2_time']/50.0), 2))

        with open(cls.mlt_lc_name, 'wb') as file_handle:
            np.savez(file_handle, **lc_files)

        # Create a database of stars using these light curves
        cls.db_name = os.path.join(cls.scratch_dir, 'test_mlt_mixed_with_dummy_db.db')

        conn = sqlite3.connect(cls.db_name)
        cursor = conn.cursor()
        cursor.execute('''CREATE TABLE mlt_test
                       (simobjid int, ra real, decl real, sedfilename text,
                        varParamStr text, parallax real, galacticAv real,
                        magNorm real)''')
        conn.commit()

        cursor.execute('''INSERT INTO mlt_test VALUES( 0, 25.0, 31.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz',
                       '{"m": "MLT",
                         "p": {"lc": "lc_1.txt", "t0": 456.2}}',
                       0.25, 2.432, 17.1)''')

        cursor.execute('''INSERT INTO mlt_test VALUES( 1, 25.2, 32.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz',
                       '{"m": "MLT",
                         "p": {"lc": "lc_1.txt", "t0": 41006.2}}',
                       0.15, 1.876, 17.2)''')

        cursor.execute('''INSERT INTO mlt_test VALUES( 2, 25.3, 10.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz',
                       '{"m": "MLT",
                         "p": {"lc": "lc_2.txt", "t0": 117.2}}',
                       0.3, 2.654, 17.3)''')

        cursor.execute('''INSERT INTO mlt_test VALUES( 3, 25.4, 11.0,
                       'lte028-5.0+0.5a+0.0.BT-Settl.spec.gz',
                       '{"m": "dummy",
                         "p": {"lc": "dummy_lc.txt"}}',
                       0.22, 2.364, 17.4)''')
        conn.commit()
        conn.close()

        cls.dummy_lc_name = os.path.join(cls.scratch_dir, 'dummy_lc.txt')

        with open(cls.dummy_lc_name, 'w') as output_file:
            for tt in np.arange(59580.0, 82000.0, 1000.0):
                output_file.write('%e %e %e %e %e %e %e\n' %
                                  (tt, 2*(tt-59580.0)/10000.0,
                                   3*(tt-59580.0)/10000.0,
                                   4*(tt-59580.0)/10000.0,
                                   5*(tt-59580.0)/10000.0,
                                   6*(tt-59580.0)/10000.0,
                                   7*(tt-59580.0)/10000.0))

    @classmethod
    def tearDownClass(cls):
        sims_clean_up()
        if os.path.exists(cls.mlt_lc_name):
            os.unlink(cls.mlt_lc_name)

        if os.path.exists(cls.db_name):
            os.unlink(cls.db_name)

        if os.path.exists(cls.dummy_lc_name):
            os.unlink(cls.dummy_lc_name)

        if os.path.exists(cls.scratch_dir):
            shutil.rmtree(cls.scratch_dir)

    def test_flare_magnitudes_mixed_with_dummy(self):
        """
        Test that we get the expected magnitudes out
        """
        db = MLT_test_DB(database=self.db_name, driver='sqlite')

        # load the quiescent SEDs of the objects in our catalog
        sed_list = SedList(['lte028-5.0+0.5a+0.0.BT-Settl.spec.gz']*4,
                           [17.1, 17.2, 17.3, 17.4],
                           galacticAvList = [2.432, 1.876, 2.654, 2.364])

        bp_dict = BandpassDict.loadTotalBandpassesFromFiles()

        # calculate the quiescent fluxes of the objects in our catalog
        baseline_fluxes = bp_dict.fluxListForSedList(sed_list)

        bb_wavelen = np.arange(100.0, 1600.0, 0.1)
        bb_flambda = blackbody_lambda(bb_wavelen*10.0, 9000.0)

        # this data is taken from the setUpClass() classmethod above
        t0_list = [456.2, 41006.2, 117.2, 10456.2]
        av_list = [2.432, 1.876, 2.654, 2.364]
        parallax_list = np.array([0.25, 0.15, 0.3, 0.22])
        distance_list = 1.0/(206265.0*radiansFromArcsec(0.001*parallax_list))
        distance_list *= 3.0857e18  # convert to cm

        dtype = np.dtype([('id', int), ('u', float), ('g', float)])

        photParams = PhotometricParameters()

        ss = Sed()

        quiet_cat_name = os.path.join(self.scratch_dir, 'mlt_mixed_with_dummy_quiet_cat.txt')
        flare_cat_name = os.path.join(self.scratch_dir, 'mlt_mixed_with_dummy_flaring_cat.txt')

        # loop over several MJDs and verify that, to within a
        # milli-mag, our flaring model gives us the magnitudes
        # expected, given the light curves specified in
        # setUpClass()
        for mjd in (59580.0, 60000.0, 70000.0, 80000.0):

            obs = ObservationMetaData(mjd=mjd)

            quiet_cat = QuiescentCatalog(db, obs_metadata=obs)
            quiet_cat.write_catalog(quiet_cat_name)

            flare_cat = FlaringCatalogDummy(db, obs_metadata=obs)
            flare_cat.scratch_dir = self.scratch_dir
            flare_cat._mlt_lc_file = self.mlt_lc_name
            flare_cat.write_catalog(flare_cat_name)

            quiescent_data = np.genfromtxt(quiet_cat_name, dtype=dtype, delimiter=',')
            flaring_data = np.genfromtxt(flare_cat_name, dtype=dtype, delimiter=',')

            self.assertGreater(len(quiescent_data), 2)
            self.assertEqual(len(quiescent_data), len(flaring_data))
            self.assertIn(3, flaring_data['id'])

            for ix in range(len(flaring_data)):
                obj_id = flaring_data['id'][ix]
                self.assertEqual(obj_id, ix)


                msg = ('failed on object %d; mjd %.2f\n u_quiet %e u_flare %e\n g_quiet %e g_flare %e' %
                       (obj_id, mjd, quiescent_data['u'][obj_id], flaring_data['u'][obj_id],
                        quiescent_data['g'][obj_id], flaring_data['g'][obj_id]))

                self.assertEqual(quiescent_data['id'][obj_id], flaring_data['id'][obj_id], msg=msg)
                self.assertAlmostEqual(ss.magFromFlux(baseline_fluxes[obj_id][0]),
                                       quiescent_data['u'][obj_id], 3, msg=msg)
                self.assertAlmostEqual(ss.magFromFlux(baseline_fluxes[obj_id][1]),
                                       quiescent_data['g'][obj_id], 3, msg=msg)
                if obj_id != 3:

                    # the models below are as specified in the
                    # setUpClass() method
                    if obj_id == 0 or obj_id == 1:
                        amp = 1.0e36
                        dt = 3652.5
                        t_min = flare_cat._survey_start - t0_list[obj_id]

                        tt = mjd - t_min
                        while tt > dt:
                            tt -= dt

                        u_flux = amp*(1.0+np.power(np.sin(tt/100.0), 2))
                        g_flux = amp*(1.0+np.power(np.cos(tt/100.0), 2))
                    elif obj_id==2:
                        amp = 2.0e35
                        dt = 365.25
                        t_min = flare_cat._survey_start - t0_list[obj_id]

                        tt = mjd - t_min
                        while tt > dt:
                            tt -= dt
                        u_flux = amp*(1.0+np.power(np.sin(tt/50.0), 2))
                        g_flux = amp*(1.0+np.power(np.cos(tt/50.0), 2))

                    # calculate the multiplicative effect of dust on a 9000K
                    # black body
                    bb_sed = Sed(wavelen=bb_wavelen, flambda=bb_flambda)
                    u_bb_flux = bb_sed.calcFlux(bp_dict['u'])
                    g_bb_flux = bb_sed.calcFlux(bp_dict['g'])
                    a_x, b_x = bb_sed.setupCCMab()
                    bb_sed.addCCMDust(a_x, b_x, A_v=av_list[obj_id])
                    u_bb_dusty_flux = bb_sed.calcFlux(bp_dict['u'])
                    g_bb_dusty_flux = bb_sed.calcFlux(bp_dict['g'])

                    dust_u = u_bb_dusty_flux/u_bb_flux
                    dust_g = g_bb_dusty_flux/g_bb_flux

                    area = 4.0*np.pi*np.power(distance_list[obj_id], 2)
                    tot_u_flux = baseline_fluxes[obj_id][0] + u_flux*dust_u*photParams.effarea/area
                    tot_g_flux = baseline_fluxes[obj_id][1] + g_flux*dust_g*photParams.effarea/area

                    self.assertAlmostEqual(ss.magFromFlux(tot_u_flux), flaring_data['u'][obj_id],
                                           3, msg=msg)
                    self.assertAlmostEqual(ss.magFromFlux(tot_g_flux), flaring_data['g'][obj_id],
                                           3, msg=msg)

                    self.assertGreater(np.abs(flaring_data['g'][obj_id]-quiescent_data['g'][obj_id]),
                                       0.001, msg=msg)
                    self.assertGreater(np.abs(flaring_data['u'][obj_id]-quiescent_data['u'][obj_id]),
                                       0.001, msg=msg)
                else:
                    self.assertAlmostEqual(flaring_data['g'][obj_id],
                                           quiescent_data['g'][obj_id]+3*(mjd-59580.0)/10000.0,
                                           3, msg=msg)
                    self.assertAlmostEqual(flaring_data['u'][obj_id],
                                           quiescent_data['u'][obj_id]+2*(mjd-59580.0)/10000.0,
                                           3, msg=msg)

        if os.path.exists(quiet_cat_name):
            os.unlink(quiet_cat_name)
        if os.path.exists(flare_cat_name):
            os.unlink(flare_cat_name)



class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
