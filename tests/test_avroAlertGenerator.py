import unittest
import os
import tempfile
import shutil
import gc
import numpy as np
import sqlite3
import numbers
import lsst.utils.tests

from lsst.utils import getPackageDir
from lsst.sims.utils.CodeUtilities import sims_clean_up
from lsst.sims.catalogs.decorators import register_method
from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.utils import AlertStellarVariabilityCatalog
from lsst.sims.catUtils.utils import StellarAlertDBObjMixin
from lsst.sims.utils import findHtmid
from lsst.sims.utils import applyProperMotion, ModifiedJulianDate
from lsst.sims.photUtils import Sed, calcSNR_m5, BandpassDict
from lsst.sims.photUtils import PhotometricParameters
from lsst.sims.coordUtils import lsst_camera
from lsst.sims.coordUtils import chipNameFromPupilCoordsLSST
from lsst.sims.catUtils.mixins import CameraCoordsLSST
from lsst.sims.catUtils.mixins import AstrometryStars
from lsst.sims.catUtils.mixins import Variability
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catalogs.decorators import compound, cached

from lsst.sims.catUtils.utils import AlertDataGenerator
from lsst.sims.catUtils.utils import AvroAlertGenerator

from lsst.sims.coordUtils import clean_up_lsst_camera

_avro_is_installed = True
try:
    from avro.io import DatumReader
    from avro.datafile import DataFileReader
except ImportError:
    _avro_is_installed = False
    pass


ROOT = os.path.abspath(os.path.dirname(__file__))


def setup_module(module):
    lsst.utils.tests.init()


class StarAlertTestDBObj_avro(StellarAlertDBObjMixin, CatalogDBObject):
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


class TestAlertsVarCatMixin_avro(object):

    @register_method('avro_test')
    def applyAlertTest(self, valid_dexes, params, expmjd, variability_cache=None):
        if len(params) == 0:
            return np.array([[], [], [], [], [], []])

        if isinstance(expmjd, numbers.Number):
            dMags_out = np.zeros((6, self.num_variable_obj(params)))
        else:
            dMags_out = np.zeros((6, self.num_variable_obj(params), len(expmjd)))

        for i_star in range(self.num_variable_obj(params)):
            if params['amp'][i_star] is not None:
                dmags = params['amp'][i_star]*np.cos(params['per'][i_star]*expmjd)
                for i_filter in range(6):
                    dMags_out[i_filter][i_star] = dmags

        return dMags_out


class TestAlertsVarCat_avro(TestAlertsVarCatMixin_avro, AlertStellarVariabilityCatalog):
    pass


class TestAlertsTruthCat_avro(TestAlertsVarCatMixin_avro, CameraCoordsLSST, AstrometryStars,
                              Variability, InstanceCatalog):
    column_outputs = ['uniqueId', 'chipName', 'dmagAlert', 'magAlert',
                      'raICRS', 'decICRS', 'xPix', 'yPix']

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


@unittest.skipIf(not _avro_is_installed, 'avro is not installed on this system')
class AvroAlertTestCase(unittest.TestCase):

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

        cls.input_dir = tempfile.mkdtemp(prefix='avroAlertGen',
                                         dir=ROOT)

        cls.star_db_name = tempfile.mktemp(prefix='avroAlertGen_star_db',
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
                    varParamStr = ('{"m":"avro_test", "p":{"amp":%.4f, "per": %.4f}}'
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

        cls.mag0_truth_dict = {}
        cls.mag0_truth_dict[0] = u_truth
        cls.mag0_truth_dict[1] = g_truth
        cls.mag0_truth_dict[2] = r_truth
        cls.mag0_truth_dict[3] = i_truth
        cls.mag0_truth_dict[4] = z_truth
        cls.mag0_truth_dict[5] = y_truth
        assert max_str_len < 500  # make sure varParamStr fits in the space alotted to it
                                  # in StarAlertTestDBObj_avro

    @classmethod
    def tearDownClass(cls):
        sims_clean_up()
        if os.path.exists(cls.star_db_name):
            os.unlink(cls.star_db_name)
        if os.path.exists(cls.input_dir):
            shutil.rmtree(cls.input_dir)

        clean_up_lsst_camera()

    def setUp(self):
        self.alert_data_output_dir = tempfile.mkdtemp(dir=ROOT, prefix='avro_gen_output')
        self.avro_out_dir = tempfile.mkdtemp(dir=ROOT, prefix='avroTestOut')

    def tearDown(self):
        for file_name in os.listdir(self.alert_data_output_dir):
            os.unlink(os.path.join(self.alert_data_output_dir, file_name))
        shutil.rmtree(self.alert_data_output_dir)

        for file_name in os.listdir(self.avro_out_dir):
            os.unlink(os.path.join(self.avro_out_dir, file_name))
        shutil.rmtree(self.avro_out_dir)

    def test_avro_alert_generation(self):
        dmag_cutoff = 0.005
        mag_name_to_int = {'u': 0, 'g': 1, 'r': 2, 'i': 3, 'z': 4, 'y': 5}

        star_db = StarAlertTestDBObj_avro(database=self.star_db_name, driver='sqlite')

        # assemble a dict of all of the alerts that need to be generated

        obshistid_list = []
        for obs in self.obs_list:
            obshistid_list.append(obs.OpsimMetaData['obsHistID'])
        obshistid_max = max(obshistid_list)
        obshistid_bits = int(np.ceil(np.log(obshistid_max)/np.log(2.0)))

        true_alert_dict = {}
        obs_dict = {}
        for obs in self.obs_list:
            obs_dict[obs.OpsimMetaData['obsHistID']] = obs
            obshistid = obs.OpsimMetaData['obsHistID']
            cat = TestAlertsTruthCat_avro(star_db, obs_metadata=obs)
            cat.camera = lsst_camera()

            for line in cat.iter_catalog():
                if line[1] is None:
                    continue

                dmag = line[2]
                mag = line[3]
                if (np.abs(dmag) > dmag_cutoff and
                    mag <= self.obs_mag_cutoff[mag_name_to_int[obs.bandpass]]):

                    alertId = (line[0] << obshistid_bits) + obshistid
                    self.assertNotIn(alertId, true_alert_dict)
                    true_alert_dict[alertId] = {}
                    true_alert_dict[alertId]['chipName'] = line[1]
                    true_alert_dict[alertId]['dmag'] = dmag
                    true_alert_dict[alertId]['mag'] = mag
                    true_alert_dict[alertId]['ra'] = np.degrees(line[4])
                    true_alert_dict[alertId]['decl'] = np.degrees(line[5])
                    true_alert_dict[alertId]['xPix'] = line[6]
                    true_alert_dict[alertId]['yPix'] = line[7]

        self.assertGreater(len(true_alert_dict), 10)

        log_file_name = tempfile.mktemp(dir=self.alert_data_output_dir, suffix='log.txt')
        alert_gen = AlertDataGenerator(testing=True)

        alert_gen.subdivide_obs(self.obs_list, htmid_level=6)

        for htmid in alert_gen.htmid_list:
            alert_gen.alert_data_from_htmid(htmid, star_db,
                                            photometry_class=TestAlertsVarCat_avro,
                                            output_prefix='alert_test',
                                            output_dir=self.alert_data_output_dir,
                                            dmag_cutoff=dmag_cutoff,
                                            log_file_name=log_file_name)

        obshistid_to_htmid = {}
        for htmid in alert_gen.htmid_list:
            for obs in alert_gen.obs_from_htmid(htmid):
                obshistid = obs.OpsimMetaData['obsHistID']
                if obshistid not in obshistid_to_htmid:
                    obshistid_to_htmid[obshistid] = []
                obshistid_to_htmid[obshistid].append(htmid)

        avro_gen = AvroAlertGenerator()
        avro_gen.load_schema(os.path.join(getPackageDir('sims_catUtils'), 'tests', 'testData', 'avroSchema'))
        sql_prefix_list = ['alert_test']
        out_prefix = 'test_avro'
        log_file_name = tempfile.mktemp(dir=self.avro_out_dir,
                                        prefix='test_avro',
                                        suffix='log.txt')
        for obshistid in obshistid_list:
            avro_gen.write_alerts(obshistid, self.alert_data_output_dir,
                                  sql_prefix_list,
                                  obshistid_to_htmid[obshistid],
                                  self.avro_out_dir, out_prefix,
                                  dmag_cutoff, lock=None,
                                  log_file_name=log_file_name)

        list_of_avro_files = os.listdir(self.avro_out_dir)
        self.assertGreater(len(list_of_avro_files), 2)
        alert_ct = 0
        dummy_sed = Sed()
        bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
        photParams = PhotometricParameters()
        diasourceId_set = set()
        for avro_file_name in list_of_avro_files:
            if avro_file_name.endswith('log.txt'):
                continue
            full_name = os.path.join(self.avro_out_dir, avro_file_name)
            with DataFileReader(open(full_name, 'rb'), DatumReader()) as data_reader:
                for alert in data_reader:
                    alert_ct += 1
                    obshistid = alert['alertId'] >> 20
                    obs = obs_dict[obshistid]
                    uniqueId = alert['diaObject']['diaObjectId']
                    true_alert_id = (uniqueId << obshistid_bits) + obshistid
                    self.assertIn(true_alert_id, true_alert_dict)
                    self.assertEqual(alert['l1dbId'], uniqueId)

                    true_alert = true_alert_dict[true_alert_id]

                    diaSource = alert['diaSource']
                    self.assertAlmostEqual(diaSource['ra'], true_alert['ra'], 10)
                    self.assertAlmostEqual(diaSource['decl'], true_alert['decl'], 10)
                    self.assertAlmostEqual(diaSource['x'], true_alert['xPix'], 3)
                    self.assertAlmostEqual(diaSource['y'], true_alert['yPix'], 3)
                    self.assertAlmostEqual(diaSource['midPointTai'], obs.mjd.TAI, 4)

                    true_tot_flux = dummy_sed.fluxFromMag(true_alert['mag'])
                    true_q_mag = true_alert['mag'] - true_alert['dmag']
                    true_q_flux = dummy_sed.fluxFromMag(true_q_mag)
                    true_dflux = true_tot_flux - true_q_flux
                    self.assertAlmostEqual(diaSource['psFlux']/true_dflux, 1.0, 6)
                    self.assertAlmostEqual(diaSource['totFlux']/true_tot_flux, 1.0, 6)
                    self.assertAlmostEqual(diaSource['diffFlux']/true_dflux, 1.0, 6)

                    true_tot_snr, gamma = calcSNR_m5(true_alert['mag'], bp_dict[obs.bandpass],
                                                     obs.m5[obs.bandpass], photParams)

                    true_q_snr, gamma = calcSNR_m5(true_q_mag, bp_dict[obs.bandpass],
                                                   self.obs_mag_cutoff[mag_name_to_int[obs.bandpass]],
                                                   photParams)

                    true_tot_err = true_tot_flux/true_tot_snr
                    true_q_err = true_q_flux/true_q_snr
                    true_diff_err = np.sqrt(true_tot_err**2 + true_q_err**2)

                    self.assertAlmostEqual(diaSource['snr']/np.abs(true_dflux/true_diff_err),
                                           1.0, 6)

                    self.assertAlmostEqual(diaSource['totFluxErr']/true_tot_err, 1.0, 6)
                    self.assertAlmostEqual(diaSource['diffFluxErr']/true_diff_err, 1.0, 6)

                    chipnum = int(true_alert['chipName'].replace('R', '').replace('S', '').
                                  replace(',', '').replace(':', '').replace(' ', ''))

                    true_ccdid = (chipnum*10**7)+obshistid
                    self.assertEqual(true_ccdid, diaSource['ccdVisitId'])
                    self.assertEqual(uniqueId, diaSource['diaObjectId'])

                    self.assertNotIn(diaSource['diaSourceId'], diasourceId_set)
                    diasourceId_set.add(diaSource['diaSourceId'])

                    diaObject = alert['diaObject']
                    obj_dex = (uniqueId//1024) - 1
                    self.assertAlmostEqual(0.001*diaObject['pmRa']/self.pmra_truth[obj_dex], 1.0, 5)
                    self.assertAlmostEqual(0.001*diaObject['pmDecl']/self.pmdec_truth[obj_dex], 1.0, 5)
                    self.assertAlmostEqual(0.001*diaObject['parallax']/self.px_truth[obj_dex], 1.0, 5)

                    (true_ra_base,
                     true_dec_base) = applyProperMotion(self.ra_truth[obj_dex],
                                                        self.dec_truth[obj_dex],
                                                        self.pmra_truth[obj_dex],
                                                        self.pmdec_truth[obj_dex],
                                                        self.px_truth[obj_dex],
                                                        self.vrad_truth[obj_dex],
                                                        mjd=ModifiedJulianDate(TAI=diaObject['radecTai']))

                    self.assertAlmostEqual(true_ra_base, diaObject['ra'], 7)
                    self.assertAlmostEqual(true_dec_base, diaObject['decl'], 7)

        self.assertEqual(alert_ct, len(true_alert_dict))

    def test_avro_alert_generation_diff_dmag(self):
        """
        Make sure everything works properly when the AlertDataGenerator
        and the AvroAlertGenerator have different dmag thresholds
        """
        dmag_cutoff_sqlite = 0.005
        dmag_cutoff_avro = 0.2
        mag_name_to_int = {'u': 0, 'g': 1, 'r': 2, 'i': 3, 'z': 4, 'y': 5}

        star_db = StarAlertTestDBObj_avro(database=self.star_db_name, driver='sqlite')

        # assemble a dict of all of the alerts that need to be generated

        obshistid_list = []
        for obs in self.obs_list:
            obshistid_list.append(obs.OpsimMetaData['obsHistID'])
        obshistid_max = max(obshistid_list)
        obshistid_bits = int(np.ceil(np.log(obshistid_max)/np.log(2.0)))

        true_alert_dict = {}
        obs_dict = {}
        ignored_sqlite = 0  # count number of alerts written to sqlite, but not avro
        for obs in self.obs_list:
            obs_dict[obs.OpsimMetaData['obsHistID']] = obs
            obshistid = obs.OpsimMetaData['obsHistID']
            cat = TestAlertsTruthCat_avro(star_db, obs_metadata=obs)
            cat.camera = lsst_camera()

            for line in cat.iter_catalog():
                if line[1] is None:
                    continue

                dmag = line[2]
                mag = line[3]
                if (np.abs(dmag) > dmag_cutoff_avro and
                    mag <= self.obs_mag_cutoff[mag_name_to_int[obs.bandpass]]):

                    alertId = (line[0] << obshistid_bits) + obshistid
                    self.assertNotIn(alertId, true_alert_dict)
                    true_alert_dict[alertId] = {}
                    true_alert_dict[alertId]['chipName'] = line[1]
                    true_alert_dict[alertId]['dmag'] = dmag
                    true_alert_dict[alertId]['mag'] = mag
                    true_alert_dict[alertId]['ra'] = np.degrees(line[4])
                    true_alert_dict[alertId]['decl'] = np.degrees(line[5])
                    true_alert_dict[alertId]['xPix'] = line[6]
                    true_alert_dict[alertId]['yPix'] = line[7]
                elif np.abs(dmag) > dmag_cutoff_sqlite:
                    ignored_sqlite += 1

        self.assertGreater(len(true_alert_dict), 10)

        self.assertGreater(ignored_sqlite, 50)  # just make sure that some sqlite
                                                # alerts were ignored by the more
                                                # stringent avro cut

        log_file_name = tempfile.mktemp(dir=self.alert_data_output_dir, suffix='log.txt')
        alert_gen = AlertDataGenerator(testing=True)

        alert_gen.subdivide_obs(self.obs_list, htmid_level=6)

        for htmid in alert_gen.htmid_list:
            alert_gen.alert_data_from_htmid(htmid, star_db,
                                            photometry_class=TestAlertsVarCat_avro,
                                            output_prefix='alert_test',
                                            output_dir=self.alert_data_output_dir,
                                            dmag_cutoff=dmag_cutoff_sqlite,
                                            log_file_name=log_file_name)

        obshistid_to_htmid = {}
        for htmid in alert_gen.htmid_list:
            for obs in alert_gen.obs_from_htmid(htmid):
                obshistid = obs.OpsimMetaData['obsHistID']
                if obshistid not in obshistid_to_htmid:
                    obshistid_to_htmid[obshistid] = []
                obshistid_to_htmid[obshistid].append(htmid)

        avro_gen = AvroAlertGenerator()
        avro_gen.load_schema(os.path.join(getPackageDir('sims_catUtils'), 'tests', 'testData', 'avroSchema'))
        sql_prefix_list = ['alert_test']
        out_prefix = 'test_avro'
        log_file_name = tempfile.mktemp(dir=self.avro_out_dir,
                                        prefix='test_avro',
                                        suffix='log.txt')
        for obshistid in obshistid_list:
            avro_gen.write_alerts(obshistid, self.alert_data_output_dir,
                                  sql_prefix_list,
                                  obshistid_to_htmid[obshistid],
                                  self.avro_out_dir, out_prefix,
                                  dmag_cutoff_avro, lock=None,
                                  log_file_name=log_file_name)

        list_of_avro_files = os.listdir(self.avro_out_dir)
        self.assertGreater(len(list_of_avro_files), 2)
        alert_ct = 0
        dummy_sed = Sed()
        bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
        photParams = PhotometricParameters()
        diasourceId_set = set()
        for avro_file_name in list_of_avro_files:
            if avro_file_name.endswith('log.txt'):
                continue
            full_name = os.path.join(self.avro_out_dir, avro_file_name)
            with DataFileReader(open(full_name, 'rb'), DatumReader()) as data_reader:
                for alert in data_reader:
                    alert_ct += 1
                    obshistid = alert['alertId'] >> 20
                    obs = obs_dict[obshistid]
                    uniqueId = alert['diaObject']['diaObjectId']
                    true_alert_id = (uniqueId << obshistid_bits) + obshistid
                    self.assertIn(true_alert_id, true_alert_dict)
                    self.assertEqual(alert['l1dbId'], uniqueId)

                    true_alert = true_alert_dict[true_alert_id]

                    diaSource = alert['diaSource']
                    self.assertAlmostEqual(diaSource['ra'], true_alert['ra'], 10)
                    self.assertAlmostEqual(diaSource['decl'], true_alert['decl'], 10)
                    self.assertAlmostEqual(diaSource['x'], true_alert['xPix'], 3)
                    self.assertAlmostEqual(diaSource['y'], true_alert['yPix'], 3)
                    self.assertAlmostEqual(diaSource['midPointTai'], obs.mjd.TAI, 4)

                    true_tot_flux = dummy_sed.fluxFromMag(true_alert['mag'])
                    true_q_mag = true_alert['mag'] - true_alert['dmag']
                    true_q_flux = dummy_sed.fluxFromMag(true_q_mag)
                    true_dflux = true_tot_flux - true_q_flux
                    self.assertAlmostEqual(diaSource['psFlux']/true_dflux, 1.0, 6)
                    self.assertAlmostEqual(diaSource['totFlux']/true_tot_flux, 1.0, 6)
                    self.assertAlmostEqual(diaSource['diffFlux']/true_dflux, 1.0, 6)

                    true_tot_snr, gamma = calcSNR_m5(true_alert['mag'], bp_dict[obs.bandpass],
                                                     obs.m5[obs.bandpass], photParams)

                    true_q_snr, gamma = calcSNR_m5(true_q_mag, bp_dict[obs.bandpass],
                                                   self.obs_mag_cutoff[mag_name_to_int[obs.bandpass]],
                                                   photParams)

                    true_tot_err = true_tot_flux/true_tot_snr
                    true_q_err = true_q_flux/true_q_snr
                    true_diff_err = np.sqrt(true_tot_err**2 + true_q_err**2)

                    self.assertAlmostEqual(diaSource['snr']/np.abs(true_dflux/true_diff_err),
                                           1.0, 6)

                    self.assertAlmostEqual(diaSource['totFluxErr']/true_tot_err, 1.0, 6)
                    self.assertAlmostEqual(diaSource['diffFluxErr']/true_diff_err, 1.0, 6)

                    chipnum = int(true_alert['chipName'].replace('R', '').replace('S', '').
                                  replace(',', '').replace(':', '').replace(' ', ''))

                    true_ccdid = (chipnum*10**7)+obshistid
                    self.assertEqual(true_ccdid, diaSource['ccdVisitId'])
                    self.assertEqual(uniqueId, diaSource['diaObjectId'])

                    self.assertNotIn(diaSource['diaSourceId'], diasourceId_set)
                    diasourceId_set.add(diaSource['diaSourceId'])

                    diaObject = alert['diaObject']
                    obj_dex = (uniqueId//1024) - 1
                    self.assertAlmostEqual(0.001*diaObject['pmRa']/self.pmra_truth[obj_dex], 1.0, 5)
                    self.assertAlmostEqual(0.001*diaObject['pmDecl']/self.pmdec_truth[obj_dex], 1.0, 5)
                    self.assertAlmostEqual(0.001*diaObject['parallax']/self.px_truth[obj_dex], 1.0, 5)

                    (true_ra_base,
                     true_dec_base) = applyProperMotion(self.ra_truth[obj_dex],
                                                        self.dec_truth[obj_dex],
                                                        self.pmra_truth[obj_dex],
                                                        self.pmdec_truth[obj_dex],
                                                        self.px_truth[obj_dex],
                                                        self.vrad_truth[obj_dex],
                                                        mjd=ModifiedJulianDate(TAI=diaObject['radecTai']))

                    self.assertAlmostEqual(true_ra_base, diaObject['ra'], 7)
                    self.assertAlmostEqual(true_dec_base, diaObject['decl'], 7)

        self.assertEqual(alert_ct, len(true_alert_dict))

    def test_avro_alert_generation_snr(self):
        """
        Test that the avroAlertGenerator selects the correct objects when you
        set a signal-to-noise cutoff
        """
        dmag_cutoff = 0.005
        snr_cutoff = 15.0
        mag_name_to_int = {'u': 0, 'g': 1, 'r': 2, 'i': 3, 'z': 4, 'y': 5}
        bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
        phot_params = PhotometricParameters()

        quiescent_m5 = {'u': 23.68, 'g': 24.89, 'r': 24.43,
                        'i': 24.0, 'z': 24.45, 'y': 22.60}

        dummy_sed = Sed()

        star_db = StarAlertTestDBObj_avro(database=self.star_db_name, driver='sqlite')

        # assemble a dict of all of the alerts that need to be generated

        obshistid_list = []
        for obs in self.obs_list:
            obshistid_list.append(obs.OpsimMetaData['obsHistID'])
        obshistid_max = max(obshistid_list)
        obshistid_bits = int(np.ceil(np.log(obshistid_max)/np.log(2.0)))

        skipped_by_snr = 0
        true_alert_dict = {}
        obs_dict = {}
        for obs in self.obs_list:
            obs_dict[obs.OpsimMetaData['obsHistID']] = obs
            obshistid = obs.OpsimMetaData['obsHistID']
            cat = TestAlertsTruthCat_avro(star_db, obs_metadata=obs)
            cat.camera = lsst_camera()

            for line in cat.iter_catalog():
                if line[1] is None:
                    continue

                dmag = line[2]
                mag = line[3]
                if (np.abs(dmag) > dmag_cutoff and
                    mag <= self.obs_mag_cutoff[mag_name_to_int[obs.bandpass]]):


                    q_mag = mag - dmag
                    q_flux = dummy_sed.fluxFromMag(q_mag)
                    q_snr, gamma = calcSNR_m5(q_mag,
                                              bp_dict[obs.bandpass],
                                              quiescent_m5[obs.bandpass],
                                              phot_params)
                    q_noise = q_flux/q_snr

                    a_flux = dummy_sed.fluxFromMag(mag)
                    a_snr, gamma = calcSNR_m5(mag,
                                              bp_dict[obs.bandpass],
                                              obs.m5[obs.bandpass],
                                              phot_params)
                    a_noise = a_flux/a_snr

                    d_noise = np.sqrt(a_noise**2+q_noise**2)
                    d_snr = np.abs((a_flux-q_flux)/d_noise)
                    if d_snr < snr_cutoff:
                        skipped_by_snr += 1
                        continue

                    alertId = (line[0] << obshistid_bits) + obshistid
                    self.assertNotIn(alertId, true_alert_dict)
                    true_alert_dict[alertId] = {}
                    true_alert_dict[alertId]['chipName'] = line[1]
                    true_alert_dict[alertId]['dmag'] = dmag
                    true_alert_dict[alertId]['mag'] = mag
                    true_alert_dict[alertId]['ra'] = np.degrees(line[4])
                    true_alert_dict[alertId]['decl'] = np.degrees(line[5])
                    true_alert_dict[alertId]['xPix'] = line[6]
                    true_alert_dict[alertId]['yPix'] = line[7]


        self.assertGreater(skipped_by_snr, 10)
        self.assertGreater(len(true_alert_dict), 10)

        log_file_name = tempfile.mktemp(dir=self.alert_data_output_dir, suffix='log.txt')
        alert_gen = AlertDataGenerator(testing=True)

        alert_gen.subdivide_obs(self.obs_list, htmid_level=6)

        for htmid in alert_gen.htmid_list:
            alert_gen.alert_data_from_htmid(htmid, star_db,
                                            photometry_class=TestAlertsVarCat_avro,
                                            output_prefix='alert_test',
                                            output_dir=self.alert_data_output_dir,
                                            dmag_cutoff=dmag_cutoff,
                                            log_file_name=log_file_name)

        obshistid_to_htmid = {}
        for htmid in alert_gen.htmid_list:
            for obs in alert_gen.obs_from_htmid(htmid):
                obshistid = obs.OpsimMetaData['obsHistID']
                if obshistid not in obshistid_to_htmid:
                    obshistid_to_htmid[obshistid] = []
                obshistid_to_htmid[obshistid].append(htmid)

        avro_gen = AvroAlertGenerator()
        avro_gen.load_schema(os.path.join(getPackageDir('sims_catUtils'), 'tests', 'testData', 'avroSchema'))
        sql_prefix_list = ['alert_test']
        out_prefix = 'test_avro'
        log_file_name = tempfile.mktemp(dir=self.avro_out_dir,
                                        prefix='test_avro',
                                        suffix='log.txt')
        for obshistid in obshistid_list:
            avro_gen.write_alerts(obshistid, self.alert_data_output_dir,
                                  sql_prefix_list,
                                  obshistid_to_htmid[obshistid],
                                  self.avro_out_dir, out_prefix,
                                  dmag_cutoff,
                                  snr_cutoff=snr_cutoff,
                                  lock=None,
                                  log_file_name=log_file_name)

        list_of_avro_files = os.listdir(self.avro_out_dir)
        self.assertGreater(len(list_of_avro_files), 2)
        alert_ct = 0
        dummy_sed = Sed()
        bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
        photParams = PhotometricParameters()
        diasourceId_set = set()
        for avro_file_name in list_of_avro_files:
            if avro_file_name.endswith('log.txt'):
                continue
            full_name = os.path.join(self.avro_out_dir, avro_file_name)
            with DataFileReader(open(full_name, 'rb'), DatumReader()) as data_reader:
                for alert in data_reader:
                    alert_ct += 1
                    obshistid = alert['alertId'] >> 20
                    obs = obs_dict[obshistid]
                    uniqueId = alert['diaObject']['diaObjectId']
                    true_alert_id = (uniqueId << obshistid_bits) + obshistid
                    self.assertIn(true_alert_id, true_alert_dict)
                    self.assertEqual(alert['l1dbId'], uniqueId)

                    true_alert = true_alert_dict[true_alert_id]

                    diaSource = alert['diaSource']
                    self.assertAlmostEqual(diaSource['ra'], true_alert['ra'], 10)
                    self.assertAlmostEqual(diaSource['decl'], true_alert['decl'], 10)
                    self.assertAlmostEqual(diaSource['x'], true_alert['xPix'], 3)
                    self.assertAlmostEqual(diaSource['y'], true_alert['yPix'], 3)
                    self.assertAlmostEqual(diaSource['midPointTai'], obs.mjd.TAI, 4)

                    true_tot_flux = dummy_sed.fluxFromMag(true_alert['mag'])
                    true_q_mag = true_alert['mag'] - true_alert['dmag']
                    true_q_flux = dummy_sed.fluxFromMag(true_q_mag)
                    true_dflux = true_tot_flux - true_q_flux
                    self.assertAlmostEqual(diaSource['psFlux']/true_dflux, 1.0, 6)
                    self.assertAlmostEqual(diaSource['totFlux']/true_tot_flux, 1.0, 6)
                    self.assertAlmostEqual(diaSource['diffFlux']/true_dflux, 1.0, 6)

                    true_tot_snr, gamma = calcSNR_m5(true_alert['mag'], bp_dict[obs.bandpass],
                                                     obs.m5[obs.bandpass], photParams)

                    true_q_snr, gamma = calcSNR_m5(true_q_mag, bp_dict[obs.bandpass],
                                                   self.obs_mag_cutoff[mag_name_to_int[obs.bandpass]],
                                                   photParams)

                    true_tot_err = true_tot_flux/true_tot_snr
                    true_q_err = true_q_flux/true_q_snr
                    true_diff_err = np.sqrt(true_tot_err**2 + true_q_err**2)

                    self.assertAlmostEqual(diaSource['snr']/np.abs(true_dflux/true_diff_err),
                                           1.0, 6)

                    self.assertAlmostEqual(diaSource['totFluxErr']/true_tot_err, 1.0, 6)
                    self.assertAlmostEqual(diaSource['diffFluxErr']/true_diff_err, 1.0, 6)

                    chipnum = int(true_alert['chipName'].replace('R', '').replace('S', '').
                                  replace(',', '').replace(':', '').replace(' ', ''))

                    true_ccdid = (chipnum*10**7)+obshistid
                    self.assertEqual(true_ccdid, diaSource['ccdVisitId'])
                    self.assertEqual(uniqueId, diaSource['diaObjectId'])

                    self.assertNotIn(diaSource['diaSourceId'], diasourceId_set)
                    diasourceId_set.add(diaSource['diaSourceId'])

                    diaObject = alert['diaObject']
                    obj_dex = (uniqueId//1024) - 1
                    self.assertAlmostEqual(0.001*diaObject['pmRa']/self.pmra_truth[obj_dex], 1.0, 5)
                    self.assertAlmostEqual(0.001*diaObject['pmDecl']/self.pmdec_truth[obj_dex], 1.0, 5)
                    self.assertAlmostEqual(0.001*diaObject['parallax']/self.px_truth[obj_dex], 1.0, 5)

                    (true_ra_base,
                     true_dec_base) = applyProperMotion(self.ra_truth[obj_dex],
                                                        self.dec_truth[obj_dex],
                                                        self.pmra_truth[obj_dex],
                                                        self.pmdec_truth[obj_dex],
                                                        self.px_truth[obj_dex],
                                                        self.vrad_truth[obj_dex],
                                                        mjd=ModifiedJulianDate(TAI=diaObject['radecTai']))

                    self.assertAlmostEqual(true_ra_base, diaObject['ra'], 7)
                    self.assertAlmostEqual(true_dec_base, diaObject['decl'], 7)

        self.assertEqual(alert_ct, len(true_alert_dict))

class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
