import unittest
import os
import numpy as np
import tempfile
import sqlite3
import shutil
import numbers
import gc
import lsst.utils.tests

from lsst.utils import getPackageDir
from lsst.sims.utils.CodeUtilities import sims_clean_up
from lsst.sims.catalogs.decorators import register_method
from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.catalogs.db import DBObject
from lsst.sims.coordUtils import lsst_camera
from lsst.sims.coordUtils import chipNameFromPupilCoordsLSST
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.utils import AlertStellarVariabilityCatalog
from lsst.sims.catUtils.utils import AlertDataGenerator
from lsst.sims.catUtils.utils import StellarAlertDBObjMixin

from lsst.sims.utils import applyProperMotion
from lsst.sims.utils import ModifiedJulianDate
from lsst.sims.utils import findHtmid
from lsst.sims.utils import angularSeparation
from lsst.sims.photUtils import Sed
from lsst.sims.coordUtils import chipNameFromRaDecLSST
from lsst.sims.coordUtils import pixelCoordsFromRaDecLSST
from lsst.sims.photUtils import calcSNR_m5, BandpassDict
from lsst.sims.photUtils import PhotometricParameters
from lsst.sims.catUtils.mixins import CameraCoordsLSST
from lsst.sims.catUtils.mixins import AstrometryStars
from lsst.sims.catUtils.mixins import Variability
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catalogs.decorators import compound, cached

from lsst.sims.coordUtils import focalPlaneCoordsFromPupilCoordsLSST
from lsst.sims.coordUtils import pupilCoordsFromFocalPlaneCoordsLSST
from lsst.sims.coordUtils import chipNameFromPupilCoordsLSST

from lsst.sims.coordUtils import clean_up_lsst_camera

ROOT = os.path.abspath(os.path.dirname(__file__))


def setup_module(module):
    lsst.utils.tests.init()



class StarAlertTestDBObj(StellarAlertDBObjMixin, CatalogDBObject):
    objid = 'star_alert'
    tableid = 'stars'
    idColKey = 'simobjid'
    raColName = 'ra'
    decColName = 'dec'
    objectTypeId = 0
    columns = [('raJ2000', 'ra*0.01745329252'),
               ('decJ2000', 'dec*0.01745329252'),
               ('parallax', 'px*0.01745329252/3600.0'),
               ('properMotionRa', 'pmra*0.01745329252/3600.0'),
               ('properMotionDec', 'pmdec*0.01745329252/3600.0'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 500)]


class TestAlertsVarCatMixin(object):

    @register_method('alert_test')
    def applyAlertTest(self, valid_dexes, params, expmjd, variability_cache=None):
        if len(params) == 0:
            return np.array([[], [], [], [], [], []])

        if isinstance(expmjd, numbers.Number):
            dmags_out = np.zeros((6, self.num_variable_obj(params)))
        else:
            dmags_out = np.zeros((6, self.num_variable_obj(params), len(expmjd)))

        for i_star in range(self.num_variable_obj(params)):
            if params['amp'][i_star] is not None:
                dmags = params['amp'][i_star]*np.cos(params['per'][i_star]*expmjd)
                for i_filter in range(6):
                    dmags_out[i_filter][i_star] = dmags

        return dmags_out


class TestAlertsVarCat(TestAlertsVarCatMixin, AlertStellarVariabilityCatalog):
    pass


class TestAlertsTruthCat(TestAlertsVarCatMixin, CameraCoordsLSST, AstrometryStars,
                         Variability, InstanceCatalog):
    column_outputs = ['uniqueId', 'chipName', 'dmagAlert', 'magAlert']

    @compound('delta_umag', 'delta_gmag', 'delta_rmag',
              'delta_imag', 'delta_zmag', 'delta_ymag')
    def get_TruthVariability(self):
        return self.applyVariability(self.column_by_name('varParamStr'))

    @cached
    def get_dmagAlert(self):
        return self.column_by_name('delta_%smag' % self.obs_metadata.bandpass)

    @cached
    def get_magAlert(self):
        return self.column_by_name('%smag' % self.obs_metadata.bandpass) + \
               self.column_by_name('dmagAlert')


class TestAlertsTruthCatSNR(TestAlertsTruthCat):

    column_outputs = ['uniqueId', 'chipName', 'dmagAlert', 'magAlert', 'snrAlert', 'q_snr', 'tot_snr']

    @compound('snrAlert', 'q_snr', 'tot_snr')
    def get_snrVals(self):
        if not hasattr(self, 'lsstBandpassDict'):
            self.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()

        if not hasattr(self, 'photParams'):
            self.photParams = PhotometricParameters()

        q_m5 = {'u': 23.68, 'g': 24.89, 'r': 24.43,
                'i': 24.0, 'z': 24.45, 'y': 22.60}

        dummy_sed = Sed()

        q_mag = self.column_by_name('magAlert') - self.column_by_name('dmagAlert')
        q_flux = dummy_sed.fluxFromMag(q_mag)
        q_snr, gamma = calcSNR_m5(q_mag, self.lsstBandpassDict[self.obs_metadata.bandpass],
                                  q_m5[self.obs_metadata.bandpass], self.photParams)
        q_noise = q_flux/q_snr

        a_flux = dummy_sed.fluxFromMag(self.column_by_name('magAlert'))
        a_snr, gamma = calcSNR_m5(self.column_by_name('magAlert'),
                                  self.lsstBandpassDict[self.obs_metadata.bandpass],
                                  self.obs_metadata.m5[self.obs_metadata.bandpass], self.photParams)
        a_noise = a_flux/a_snr

        d_noise = np.sqrt(a_noise**2 + q_noise**2)
        d_flux = a_flux-q_flux
        d_snr = np.abs(d_flux/d_noise)

        return np.array([d_snr, q_snr, a_snr])


class AlertDataGeneratorTestCase(unittest.TestCase):

    longMessage = True

    @classmethod
    def setUpClass(cls):
        print('setting up %s' % sims_clean_up.targets)

        # These represent the dimmest magnitudes at which objects
        # are considered visible in each of the LSST filters
        # (taken from Table 2 of the overview paper)
        cls.obs_mag_cutoff = (23.68, 24.89, 24.43, 24.0, 24.45, 22.60)

        cls.opsim_db = os.path.join(getPackageDir('sims_data'),
                                    'OpSimData',
                                    'opsimblitz1_1133_sqlite.db')

        rng = np.random.RandomState(8123)

        obs_gen = ObservationMetaDataGenerator(database=cls.opsim_db)
        cls.obs_list = obs_gen.getObservationMetaData(night=(0, 2))
        cls.obs_list = rng.choice(cls.obs_list, 10, replace=False)
        fieldid_list = []
        for obs in cls.obs_list:
            fieldid_list.append(obs.OpsimMetaData['fieldID'])

        # make sure we have selected observations such that the
        # same field is revisited more than once
        assert len(np.unique(fieldid_list)) < len(fieldid_list)

        cls.input_dir = tempfile.mkdtemp(prefix='alertDataGen',
                                         dir=ROOT)

        cls.star_db_name = tempfile.mktemp(prefix='alertDataGen_star_db',
                                           dir=cls.input_dir,
                                           suffix='.db')

        conn = sqlite3.connect(cls.star_db_name)
        cursor = conn.cursor()
        cursor.execute('''CREATE TABLE stars
                          (simobjid int, htmid int, ra real, dec real,
                           umag real, gmag real, rmag real,
                           imag real, zmag real, ymag real,
                           px real, pmra real, pmdec real,
                           vrad real, varParamStr text)''')
        conn.commit()

        n_stars = 10

        cls.ra_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)
        cls.dec_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)
        u_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)
        g_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)
        r_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)
        i_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)
        z_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)
        y_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)
        cls.px_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)
        cls.pmra_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)
        cls.pmdec_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)
        cls.vrad_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)
        cls.amp_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)
        cls.period_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)

        id_offset = -n_stars
        for obs in cls.obs_list:
            id_offset += n_stars
            ra_0 = obs.pointingRA
            dec_0 = obs.pointingDec
            rr = rng.random_sample(n_stars)
            theta = rng.random_sample(n_stars)*2.0*np.pi
            ra = ra_0 + rr*np.cos(theta)
            dec = dec_0 + rr*np.sin(theta)
            var_period = rng.random_sample(n_stars)*0.25
            var_amp = rng.random_sample(n_stars)*1.0 + 0.01

            subset = rng.randint(0, high=len(var_amp)-1, size=3)
            var_amp[subset[:2]] = 0.0
            var_amp[subset[-1]] = -1.0

            umag = rng.random_sample(n_stars)*5.0 + 15.0
            gmag = rng.random_sample(n_stars)*5.0 + 15.0
            rmag = rng.random_sample(n_stars)*5.0 + 15.0
            imag = rng.random_sample(n_stars)*5.0 + 15.0
            zmag = rng.random_sample(n_stars)*5.0 + 15.0
            ymag = rng.random_sample(n_stars)*5.0 + 15.0
            px = rng.random_sample(n_stars)*0.1  # say it is arcsec
            pmra = rng.random_sample(n_stars)*50.0+100.0  # say it is arcsec/yr
            pmdec = rng.random_sample(n_stars)*50.0+100.0  # say it is arcsec/yr
            vrad = rng.random_sample(n_stars)*600.0 - 300.0

            subset = rng.randint(0, high=n_stars-1, size=3)
            umag[subset] = 40.0
            gmag[subset] = 40.0
            rmag[subset] = 40.0
            imag[subset] = 40.0
            zmag[subset] = 40.0
            ymag[subset] = 40.0

            cls.ra_truth[id_offset:id_offset+n_stars] = np.round(ra, decimals=6)
            cls.dec_truth[id_offset:id_offset+n_stars] = np.round(dec, decimals=6)
            u_truth[id_offset:id_offset+n_stars] = np.round(umag, decimals=4)
            g_truth[id_offset:id_offset+n_stars] = np.round(gmag, decimals=4)
            r_truth[id_offset:id_offset+n_stars] = np.round(rmag, decimals=4)
            i_truth[id_offset:id_offset+n_stars] = np.round(imag, decimals=4)
            z_truth[id_offset:id_offset+n_stars] = np.round(zmag, decimals=4)
            y_truth[id_offset:id_offset+n_stars] = np.round(ymag, decimals=4)
            cls.px_truth[id_offset:id_offset+n_stars] = np.round(px, decimals=4)
            cls.pmra_truth[id_offset:id_offset+n_stars] = np.round(pmra, decimals=4)
            cls.pmdec_truth[id_offset:id_offset+n_stars] = np.round(pmdec, decimals=4)
            cls.vrad_truth[id_offset:id_offset+n_stars] = np.round(vrad, decimals=4)
            cls.amp_truth[id_offset:id_offset+n_stars] = np.round(var_amp, decimals=4)
            cls.period_truth[id_offset:id_offset+n_stars] = np.round(var_period, decimals=4)

            max_str_len = -1

            for i_star in range(n_stars):
                if var_amp[i_star] >= -0.1:
                    varParamStr = ('{"m":"alert_test", "p":{"amp":%.4f, "per": %.4f}}'
                                   % (var_amp[i_star], var_period[i_star]))
                else:
                    varParamStr = 'None'

                if len(varParamStr) > max_str_len:
                    max_str_len = len(varParamStr)

                htmid = findHtmid(ra[i_star], dec[i_star], 21)

                query = ('''INSERT INTO stars VALUES(%d, %d, %.6f, %.6f,
                                                    %.4f, %.4f, %.4f, %.4f, %.4f, %.4f,
                                                    %.4f, %.4f, %.4f, %.4f, '%s')'''
                         % (i_star+id_offset+1, htmid, ra[i_star], dec[i_star],
                            umag[i_star], gmag[i_star], rmag[i_star],
                            imag[i_star], zmag[i_star], ymag[i_star],
                            px[i_star], pmra[i_star], pmdec[i_star],
                            vrad[i_star], varParamStr))

                cursor.execute(query)
        conn.commit()
        conn.close()

        assert max_str_len<500

        cls.mag0_truth_dict = {}
        cls.mag0_truth_dict[0] = u_truth
        cls.mag0_truth_dict[1] = g_truth
        cls.mag0_truth_dict[2] = r_truth
        cls.mag0_truth_dict[3] = i_truth
        cls.mag0_truth_dict[4] = z_truth
        cls.mag0_truth_dict[5] = y_truth


    @classmethod
    def tearDownClass(cls):
        sims_clean_up()
        if os.path.exists(cls.star_db_name):
            os.unlink(cls.star_db_name)
        if os.path.exists(cls.input_dir):
            shutil.rmtree(cls.input_dir)

        clean_up_lsst_camera()


    def test_alert_data_generation(self):
        """
        Test that the AlertDataGenerator generates the alerts
        it is supposed to by comparing to InstanceCatalogs of
        the objects being simulated and doing a brute force comparison
        """

        dmag_cutoff = 0.005
        mag_name_to_int = {'u': 0, 'g': 1, 'r': 2, 'i': 3, 'z' : 4, 'y': 5}

        star_db = StarAlertTestDBObj(database=self.star_db_name, driver='sqlite')

        # assemble the true light curves for each object; we need to figure out
        # if their np.max(dMag) ever goes over dmag_cutoff; then we will know if
        # we are supposed to simulate them
        true_lc_dict = {}
        true_lc_obshistid_dict = {}
        is_visible_dict = {}
        obs_dict = {}
        max_obshistid = -1
        n_total_observations = 0
        for obs in self.obs_list:
            obs_dict[obs.OpsimMetaData['obsHistID']] = obs
            obshistid = obs.OpsimMetaData['obsHistID']
            if obshistid > max_obshistid:
                max_obshistid = obshistid
            cat = TestAlertsTruthCat(star_db, obs_metadata=obs)

            for line in cat.iter_catalog():
                if line[1] is None:
                    continue

                n_total_observations += 1
                if line[0] not in true_lc_dict:
                    true_lc_dict[line[0]] = {}
                    true_lc_obshistid_dict[line[0]] = []

                true_lc_dict[line[0]][obshistid] = line[2]
                true_lc_obshistid_dict[line[0]].append(obshistid)

                if line[0] not in is_visible_dict:
                    is_visible_dict[line[0]] = False

                if line[3] <= self.obs_mag_cutoff[mag_name_to_int[obs.bandpass]]:
                    is_visible_dict[line[0]] = True

        obshistid_bits = int(np.ceil(np.log(max_obshistid)/np.log(2)))

        skipped_due_to_mag = 0

        objects_to_simulate = []
        obshistid_unqid_set = set()
        for obj_id in true_lc_dict:

            dmag_max = -1.0
            for obshistid in true_lc_dict[obj_id]:
                if np.abs(true_lc_dict[obj_id][obshistid]) > dmag_max:
                    dmag_max = np.abs(true_lc_dict[obj_id][obshistid])

            if dmag_max >= dmag_cutoff:
                if not is_visible_dict[obj_id]:
                    skipped_due_to_mag += 1
                    continue

                objects_to_simulate.append(obj_id)
                for obshistid in true_lc_obshistid_dict[obj_id]:
                    obshistid_unqid_set.add((obj_id << obshistid_bits) + obshistid)

        self.assertGreater(len(objects_to_simulate), 10)
        self.assertGreater(skipped_due_to_mag, 0)

        output_dir = tempfile.mkdtemp(dir=ROOT, prefix='alert_gen_output')
        log_file_name = tempfile.mktemp(dir=output_dir, suffix='log.txt')
        alert_gen = AlertDataGenerator(testing=True)

        alert_gen.subdivide_obs(self.obs_list, htmid_level=6)

        for htmid in alert_gen.htmid_list:
            alert_gen.alert_data_from_htmid(htmid, star_db,
                                            photometry_class=TestAlertsVarCat,
                                            output_prefix='alert_test',
                                            output_dir=output_dir,
                                            dmag_cutoff=dmag_cutoff,
                                            log_file_name=log_file_name)

        dummy_sed = Sed()

        bp_dict = BandpassDict.loadTotalBandpassesFromFiles()

        phot_params = PhotometricParameters()

        # First, verify that the contents of the sqlite files are all correct

        n_tot_simulated = 0

        alert_query = 'SELECT alert.uniqueId, alert.obshistId, meta.TAI, '
        alert_query += 'meta.band, quiescent.flux, alert.dflux, '
        alert_query += 'quiescent.snr, alert.snr, '
        alert_query += 'alert.ra, alert.dec, alert.chipNum, '
        alert_query += 'alert.xPix, alert.yPix, ast.pmRA, ast.pmDec, '
        alert_query += 'ast.parallax '
        alert_query += 'FROM alert_data AS alert '
        alert_query += 'INNER JOIN metadata AS meta ON meta.obshistId=alert.obshistId '
        alert_query += 'INNER JOIN quiescent_flux AS quiescent '
        alert_query += 'ON quiescent.uniqueId=alert.uniqueId '
        alert_query += 'AND quiescent.band=meta.band '
        alert_query += 'INNER JOIN baseline_astrometry AS ast '
        alert_query += 'ON ast.uniqueId=alert.uniqueId'

        alert_dtype = np.dtype([('uniqueId', int), ('obshistId', int),
                                ('TAI', float), ('band', int),
                                ('q_flux', float), ('dflux', float),
                                ('q_snr', float), ('tot_snr', float),
                                ('ra', float), ('dec', float),
                                ('chipNum', int), ('xPix', float), ('yPix', float),
                                ('pmRA', float), ('pmDec', float), ('parallax', float)])

        sqlite_file_list = os.listdir(output_dir)

        n_tot_simulated = 0
        obshistid_unqid_simulated_set = set()
        for file_name in sqlite_file_list:
            if not file_name.endswith('db'):
                continue
            full_name = os.path.join(output_dir, file_name)
            self.assertTrue(os.path.exists(full_name))
            alert_db = DBObject(full_name, driver='sqlite')
            alert_data = alert_db.execute_arbitrary(alert_query, dtype=alert_dtype)
            if len(alert_data) == 0:
                continue

            mjd_list = ModifiedJulianDate.get_list(TAI=alert_data['TAI'])
            for i_obj in range(len(alert_data)):
                n_tot_simulated += 1
                obshistid_unqid_simulated_set.add((alert_data['uniqueId'][i_obj] << obshistid_bits) +
                                                  alert_data['obshistId'][i_obj])

                unq = alert_data['uniqueId'][i_obj]
                obj_dex = (unq//1024)-1
                self.assertAlmostEqual(self.pmra_truth[obj_dex], 0.001*alert_data['pmRA'][i_obj], 4)
                self.assertAlmostEqual(self.pmdec_truth[obj_dex], 0.001*alert_data['pmDec'][i_obj], 4)
                self.assertAlmostEqual(self.px_truth[obj_dex], 0.001*alert_data['parallax'][i_obj], 4)

                ra_truth, dec_truth = applyProperMotion(self.ra_truth[obj_dex], self.dec_truth[obj_dex],
                                                        self.pmra_truth[obj_dex], self.pmdec_truth[obj_dex],
                                                        self.px_truth[obj_dex], self.vrad_truth[obj_dex],
                                                        mjd=mjd_list[i_obj])
                distance = angularSeparation(ra_truth, dec_truth,
                                             alert_data['ra'][i_obj], alert_data['dec'][i_obj])

                distance_arcsec = 3600.0*distance
                msg = '\ntruth: %e %e\nalert: %e %e\n' % (ra_truth, dec_truth,
                                                          alert_data['ra'][i_obj],
                                                          alert_data['dec'][i_obj])

                self.assertLess(distance_arcsec, 0.0005, msg=msg)

                obs = obs_dict[alert_data['obshistId'][i_obj]]

                chipname = chipNameFromRaDecLSST(self.ra_truth[obj_dex], self.dec_truth[obj_dex],
                                                 pm_ra=self.pmra_truth[obj_dex],
                                                 pm_dec=self.pmdec_truth[obj_dex],
                                                 parallax=self.px_truth[obj_dex],
                                                 v_rad=self.vrad_truth[obj_dex],
                                                 obs_metadata=obs,
                                                 band=obs.bandpass)

                chipnum = int(chipname.replace('R', '').replace('S', '').
                              replace(' ', '').replace(';', '').replace(',', '').
                              replace(':', ''))

                self.assertEqual(chipnum, alert_data['chipNum'][i_obj])

                xpix, ypix = pixelCoordsFromRaDecLSST(self.ra_truth[obj_dex], self.dec_truth[obj_dex],
                                                      pm_ra=self.pmra_truth[obj_dex],
                                                      pm_dec=self.pmdec_truth[obj_dex],
                                                      parallax=self.px_truth[obj_dex],
                                                      v_rad=self.vrad_truth[obj_dex],
                                                      obs_metadata=obs,
                                                      band=obs.bandpass)

                self.assertAlmostEqual(alert_data['xPix'][i_obj], xpix, 4)
                self.assertAlmostEqual(alert_data['yPix'][i_obj], ypix, 4)

                dmag_sim = -2.5*np.log10(1.0+alert_data['dflux'][i_obj]/alert_data['q_flux'][i_obj])
                self.assertAlmostEqual(true_lc_dict[alert_data['uniqueId'][i_obj]][alert_data['obshistId'][i_obj]],
                                       dmag_sim, 3)

                mag_name = ('u', 'g', 'r', 'i', 'z', 'y')[alert_data['band'][i_obj]]
                m5 = obs.m5[mag_name]

                q_mag = dummy_sed.magFromFlux(alert_data['q_flux'][i_obj])
                self.assertAlmostEqual(self.mag0_truth_dict[alert_data['band'][i_obj]][obj_dex],
                                       q_mag, 4)

                snr, gamma = calcSNR_m5(self.mag0_truth_dict[alert_data['band'][i_obj]][obj_dex],
                                        bp_dict[mag_name],
                                        self.obs_mag_cutoff[alert_data['band'][i_obj]],
                                        phot_params)

                self.assertAlmostEqual(snr/alert_data['q_snr'][i_obj], 1.0, 4)

                tot_mag = self.mag0_truth_dict[alert_data['band'][i_obj]][obj_dex] + \
                          true_lc_dict[alert_data['uniqueId'][i_obj]][alert_data['obshistId'][i_obj]]

                snr, gamma = calcSNR_m5(tot_mag, bp_dict[mag_name],
                                        m5, phot_params)
                self.assertAlmostEqual(snr/alert_data['tot_snr'][i_obj], 1.0, 4)

        for val in obshistid_unqid_set:
            self.assertIn(val, obshistid_unqid_simulated_set)
        self.assertEqual(len(obshistid_unqid_set), len(obshistid_unqid_simulated_set))

        astrometry_query = 'SELECT uniqueId, ra, dec, TAI '
        astrometry_query += 'FROM baseline_astrometry'
        astrometry_dtype = np.dtype([('uniqueId', int),
                                     ('ra', float),
                                     ('dec', float),
                                     ('TAI', float)])

        tai_list = []
        for obs in self.obs_list:
            tai_list.append(obs.mjd.TAI)
        tai_list = np.array(tai_list)

        n_tot_ast_simulated = 0
        for file_name in sqlite_file_list:
            if not file_name.endswith('db'):
                continue
            full_name = os.path.join(output_dir, file_name)
            self.assertTrue(os.path.exists(full_name))
            alert_db = DBObject(full_name, driver='sqlite')
            astrometry_data = alert_db.execute_arbitrary(astrometry_query, dtype=astrometry_dtype)

            if len(astrometry_data) == 0:
                continue

            mjd_list = ModifiedJulianDate.get_list(TAI=astrometry_data['TAI'])
            for i_obj in range(len(astrometry_data)):
                n_tot_ast_simulated += 1
                obj_dex = (astrometry_data['uniqueId'][i_obj]//1024) - 1
                ra_truth, dec_truth = applyProperMotion(self.ra_truth[obj_dex], self.dec_truth[obj_dex],
                                                        self.pmra_truth[obj_dex], self.pmdec_truth[obj_dex],
                                                        self.px_truth[obj_dex], self.vrad_truth[obj_dex],
                                                        mjd=mjd_list[i_obj])

                distance = angularSeparation(ra_truth, dec_truth,
                                             astrometry_data['ra'][i_obj],
                                             astrometry_data['dec'][i_obj])

                self.assertLess(3600.0*distance, 0.0005)

        del alert_gen
        gc.collect()
        self.assertGreater(n_tot_simulated, 10)
        self.assertGreater(len(obshistid_unqid_simulated_set), 10)
        self.assertLess(len(obshistid_unqid_simulated_set), n_total_observations)
        self.assertGreater(n_tot_ast_simulated, 0)

        out_file_list = os.listdir(output_dir)
        for file_name in out_file_list:
            os.unlink(os.path.join(output_dir, file_name))
        shutil.rmtree(output_dir)


    def test_alert_data_generation_snr(self):
        """
        Test that the AlertDataGenerator generates the alerts
        it is supposed to by comparing to InstanceCatalogs of
        the objects being simulated and doing a brute force comparison.

        Add a signal-to-noise ratio cut-off in the alert generation step.
        """

        dmag_cutoff = 0.005
        snr_cutoff = 15.0   # the simulated light curves in this test
                            # have ridiculous SNR

        mag_name_to_int = {'u': 0, 'g': 1, 'r': 2, 'i': 3, 'z' : 4, 'y': 5}

        star_db = StarAlertTestDBObj(database=self.star_db_name, driver='sqlite')

        # assemble the true light curves for each object; we need to figure out
        # if their np.max(dMag) ever goes over dmag_cutoff; then we will know if
        # we are supposed to simulate them
        true_lc_dict = {}
        true_lc_snr_dict ={}
        true_lc_q_snr_dict = {}
        true_lc_tot_snr_dict = {}
        true_lc_obshistid_dict = {}
        is_visible_dict = {}
        makes_snr_cut_dict = {}
        max_snr_dict = {}
        obs_dict = {}
        max_obshistid = -1
        n_total_observations = 0
        for obs in self.obs_list:
            obs_dict[obs.OpsimMetaData['obsHistID']] = obs
            obshistid = obs.OpsimMetaData['obsHistID']
            if obshistid > max_obshistid:
                max_obshistid = obshistid
            cat = TestAlertsTruthCatSNR(star_db, obs_metadata=obs)

            for line in cat.iter_catalog():
                if line[1] is None:
                    continue

                n_total_observations += 1
                if line[0] not in true_lc_dict:
                    true_lc_dict[line[0]] = {}
                    true_lc_snr_dict[line[0]] = {}
                    true_lc_q_snr_dict[line[0]] ={}
                    true_lc_tot_snr_dict[line[0]] = {}
                    true_lc_obshistid_dict[line[0]] = []

                true_lc_dict[line[0]][obshistid] = line[2]
                true_lc_snr_dict[line[0]][obshistid] = line[4]
                true_lc_q_snr_dict[line[0]][obshistid] = line[5]
                true_lc_tot_snr_dict[line[0]][obshistid] = line[6]
                true_lc_obshistid_dict[line[0]].append(obshistid)

                if line[0] not in is_visible_dict:
                    is_visible_dict[line[0]] = False

                if line[3] <= self.obs_mag_cutoff[mag_name_to_int[obs.bandpass]]:
                    is_visible_dict[line[0]] = True

                if line[0] not in makes_snr_cut_dict:
                    makes_snr_cut_dict[line[0]] = False

                if line[4] > snr_cutoff:
                    makes_snr_cut_dict[line[0]] = True

        obshistid_bits = int(np.ceil(np.log(max_obshistid)/np.log(2)))

        skipped_due_to_mag = 0
        skipped_due_to_snr = 0

        objects_to_simulate = []
        obshistid_unqid_set = set()
        for obj_id in true_lc_dict:

            dmag_max = -1.0
            for obshistid in true_lc_dict[obj_id]:
                if np.abs(true_lc_dict[obj_id][obshistid]) > dmag_max:
                    dmag_max = np.abs(true_lc_dict[obj_id][obshistid])

            if dmag_max >= dmag_cutoff:
                if not is_visible_dict[obj_id]:
                    skipped_due_to_mag += 1
                    continue

                if not makes_snr_cut_dict[obj_id]:
                    skipped_due_to_snr += 1
                    continue

                objects_to_simulate.append(obj_id)
                for obshistid in true_lc_obshistid_dict[obj_id]:
                    obshistid_unqid_set.add((obj_id << obshistid_bits) + obshistid)

        self.assertGreater(len(objects_to_simulate), 10)
        self.assertGreater(skipped_due_to_snr, 0)
        self.assertGreater(skipped_due_to_mag, 0)

        output_dir = tempfile.mkdtemp(dir=ROOT, prefix='alert_gen_output')
        log_file_name = tempfile.mktemp(dir=output_dir, suffix='log.txt')
        alert_gen = AlertDataGenerator(testing=True)

        alert_gen.subdivide_obs(self.obs_list, htmid_level=6)

        for htmid in alert_gen.htmid_list:
            alert_gen.alert_data_from_htmid(htmid, star_db,
                                            photometry_class=TestAlertsVarCat,
                                            output_prefix='alert_test',
                                            output_dir=output_dir,
                                            dmag_cutoff=dmag_cutoff,
                                            snr_cutoff=snr_cutoff,
                                            log_file_name=log_file_name)

        dummy_sed = Sed()

        bp_dict = BandpassDict.loadTotalBandpassesFromFiles()

        phot_params = PhotometricParameters()

        # First, verify that the contents of the sqlite files are all correct

        n_tot_simulated = 0

        alert_query = 'SELECT alert.uniqueId, alert.obshistId, meta.TAI, '
        alert_query += 'meta.band, quiescent.flux, alert.dflux, '
        alert_query += 'quiescent.snr, alert.snr, '
        alert_query += 'alert.ra, alert.dec, alert.chipNum, '
        alert_query += 'alert.xPix, alert.yPix, ast.pmRA, ast.pmDec, '
        alert_query += 'ast.parallax '
        alert_query += 'FROM alert_data AS alert '
        alert_query += 'INNER JOIN metadata AS meta ON meta.obshistId=alert.obshistId '
        alert_query += 'INNER JOIN quiescent_flux AS quiescent '
        alert_query += 'ON quiescent.uniqueId=alert.uniqueId '
        alert_query += 'AND quiescent.band=meta.band '
        alert_query += 'INNER JOIN baseline_astrometry AS ast '
        alert_query += 'ON ast.uniqueId=alert.uniqueId'

        alert_dtype = np.dtype([('uniqueId', int), ('obshistId', int),
                                ('TAI', float), ('band', int),
                                ('q_flux', float), ('dflux', float),
                                ('q_snr', float), ('tot_snr', float),
                                ('ra', float), ('dec', float),
                                ('chipNum', int), ('xPix', float), ('yPix', float),
                                ('pmRA', float), ('pmDec', float), ('parallax', float)])

        sqlite_file_list = os.listdir(output_dir)

        n_tot_simulated = 0
        obshistid_unqid_simulated_set = set()
        for file_name in sqlite_file_list:
            if not file_name.endswith('db'):
                continue
            full_name = os.path.join(output_dir, file_name)
            self.assertTrue(os.path.exists(full_name))
            alert_db = DBObject(full_name, driver='sqlite')
            alert_data = alert_db.execute_arbitrary(alert_query, dtype=alert_dtype)
            if len(alert_data) == 0:
                continue

            mjd_list = ModifiedJulianDate.get_list(TAI=alert_data['TAI'])
            for i_obj in range(len(alert_data)):
                n_tot_simulated += 1
                obshistid_unqid_simulated_set.add((alert_data['uniqueId'][i_obj] << obshistid_bits) +
                                                  alert_data['obshistId'][i_obj])

                unq = alert_data['uniqueId'][i_obj]
                obj_dex = (unq//1024)-1
                self.assertAlmostEqual(self.pmra_truth[obj_dex], 0.001*alert_data['pmRA'][i_obj], 4)
                self.assertAlmostEqual(self.pmdec_truth[obj_dex], 0.001*alert_data['pmDec'][i_obj], 4)
                self.assertAlmostEqual(self.px_truth[obj_dex], 0.001*alert_data['parallax'][i_obj], 4)

                ra_truth, dec_truth = applyProperMotion(self.ra_truth[obj_dex], self.dec_truth[obj_dex],
                                                        self.pmra_truth[obj_dex], self.pmdec_truth[obj_dex],
                                                        self.px_truth[obj_dex], self.vrad_truth[obj_dex],
                                                        mjd=mjd_list[i_obj])
                distance = angularSeparation(ra_truth, dec_truth,
                                             alert_data['ra'][i_obj], alert_data['dec'][i_obj])

                distance_arcsec = 3600.0*distance
                msg = '\ntruth: %e %e\nalert: %e %e\n' % (ra_truth, dec_truth,
                                                          alert_data['ra'][i_obj],
                                                          alert_data['dec'][i_obj])

                self.assertLess(distance_arcsec, 0.0005, msg=msg)

                obs = obs_dict[alert_data['obshistId'][i_obj]]

                chipname = chipNameFromRaDecLSST(self.ra_truth[obj_dex], self.dec_truth[obj_dex],
                                                 pm_ra=self.pmra_truth[obj_dex],
                                                 pm_dec=self.pmdec_truth[obj_dex],
                                                 parallax=self.px_truth[obj_dex],
                                                 v_rad=self.vrad_truth[obj_dex],
                                                 obs_metadata=obs)

                chipnum = int(chipname.replace('R', '').replace('S', '').
                              replace(' ', '').replace(';', '').replace(',', '').
                              replace(':', ''))

                self.assertEqual(chipnum, alert_data['chipNum'][i_obj])

                xpix, ypix = pixelCoordsFromRaDecLSST(self.ra_truth[obj_dex], self.dec_truth[obj_dex],
                                                      pm_ra=self.pmra_truth[obj_dex],
                                                      pm_dec=self.pmdec_truth[obj_dex],
                                                      parallax=self.px_truth[obj_dex],
                                                      v_rad=self.vrad_truth[obj_dex],
                                                      obs_metadata=obs)

                self.assertAlmostEqual(alert_data['xPix'][i_obj], xpix, 4)
                self.assertAlmostEqual(alert_data['yPix'][i_obj], ypix, 4)

                dmag_sim = -2.5*np.log10(1.0+alert_data['dflux'][i_obj]/alert_data['q_flux'][i_obj])
                self.assertAlmostEqual(true_lc_dict[alert_data['uniqueId'][i_obj]][alert_data['obshistId'][i_obj]],
                                       dmag_sim, 3)

                q_noise = alert_data['q_flux'][i_obj]/alert_data['q_snr'][i_obj]
                tot_noise = (alert_data['dflux'][i_obj]+alert_data['q_flux'][i_obj])/alert_data['tot_snr'][i_obj]
                d_noise = np.sqrt(q_noise**2 + tot_noise**2)
                d_snr = np.abs(alert_data['dflux'][i_obj]/d_noise)

                i0 = alert_data['uniqueId'][i_obj]
                i1 = alert_data['obshistId'][i_obj]
                msg = '\nq_snr %e %e\n' % (alert_data['q_snr'][i_obj], true_lc_q_snr_dict[i0][i1])
                msg += 'tot_snr %e %e\n' % (alert_data['tot_snr'][i_obj], true_lc_tot_snr_dict[i0][i1])

                self.assertAlmostEqual(true_lc_snr_dict[alert_data['uniqueId'][i_obj]][alert_data['obshistId'][i_obj]],
                                       d_snr, 3, msg=msg)

                mag_name = ('u', 'g', 'r', 'i', 'z', 'y')[alert_data['band'][i_obj]]
                m5 = obs.m5[mag_name]

                q_mag = dummy_sed.magFromFlux(alert_data['q_flux'][i_obj])
                self.assertAlmostEqual(self.mag0_truth_dict[alert_data['band'][i_obj]][obj_dex],
                                       q_mag, 4)

                snr, gamma = calcSNR_m5(self.mag0_truth_dict[alert_data['band'][i_obj]][obj_dex],
                                        bp_dict[mag_name],
                                        self.obs_mag_cutoff[alert_data['band'][i_obj]],
                                        phot_params)

                self.assertAlmostEqual(snr/alert_data['q_snr'][i_obj], 1.0, 4)

                tot_mag = self.mag0_truth_dict[alert_data['band'][i_obj]][obj_dex] + \
                          true_lc_dict[alert_data['uniqueId'][i_obj]][alert_data['obshistId'][i_obj]]

                snr, gamma = calcSNR_m5(tot_mag, bp_dict[mag_name],
                                        m5, phot_params)
                self.assertAlmostEqual(snr/alert_data['tot_snr'][i_obj], 1.0, 4)

        for val in obshistid_unqid_set:
            self.assertIn(val, obshistid_unqid_simulated_set)
        self.assertEqual(len(obshistid_unqid_set), len(obshistid_unqid_simulated_set))

        astrometry_query = 'SELECT uniqueId, ra, dec, TAI '
        astrometry_query += 'FROM baseline_astrometry'
        astrometry_dtype = np.dtype([('uniqueId', int),
                                     ('ra', float),
                                     ('dec', float),
                                     ('TAI', float)])

        tai_list = []
        for obs in self.obs_list:
            tai_list.append(obs.mjd.TAI)
        tai_list = np.array(tai_list)

        n_tot_ast_simulated = 0
        for file_name in sqlite_file_list:
            if not file_name.endswith('db'):
                continue
            full_name = os.path.join(output_dir, file_name)
            self.assertTrue(os.path.exists(full_name))
            alert_db = DBObject(full_name, driver='sqlite')
            astrometry_data = alert_db.execute_arbitrary(astrometry_query, dtype=astrometry_dtype)

            if len(astrometry_data) == 0:
                continue

            mjd_list = ModifiedJulianDate.get_list(TAI=astrometry_data['TAI'])
            for i_obj in range(len(astrometry_data)):
                n_tot_ast_simulated += 1
                obj_dex = (astrometry_data['uniqueId'][i_obj]//1024) - 1
                ra_truth, dec_truth = applyProperMotion(self.ra_truth[obj_dex], self.dec_truth[obj_dex],
                                                        self.pmra_truth[obj_dex], self.pmdec_truth[obj_dex],
                                                        self.px_truth[obj_dex], self.vrad_truth[obj_dex],
                                                        mjd=mjd_list[i_obj])

                distance = angularSeparation(ra_truth, dec_truth,
                                             astrometry_data['ra'][i_obj],
                                             astrometry_data['dec'][i_obj])

                self.assertLess(3600.0*distance, 0.0005)

        del alert_gen
        gc.collect()
        self.assertGreater(n_tot_simulated, 10)
        self.assertGreater(len(obshistid_unqid_simulated_set), 10)
        self.assertLess(len(obshistid_unqid_simulated_set), n_total_observations)
        self.assertGreater(n_tot_ast_simulated, 0)

        out_file_list = os.listdir(output_dir)
        for file_name in out_file_list:
            os.unlink(os.path.join(output_dir, file_name))
        shutil.rmtree(output_dir)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
