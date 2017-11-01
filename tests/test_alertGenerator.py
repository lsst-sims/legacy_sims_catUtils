import unittest
import os
import numpy as np
import tempfile
import sqlite3
import shutil
import numbers
import gc
import h5py
import lsst.utils.tests

from lsst.utils import getPackageDir
from lsst.sims.utils.CodeUtilities import sims_clean_up
from lsst.sims.catalogs.decorators import register_method
from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.coordUtils import lsst_camera
from lsst.sims.coordUtils import chipNameFromPupilCoordsLSST
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.utils import AlertStellarVariabilityCatalog
from lsst.sims.catUtils.utils import AlertDataGenerator

from lsst.sims.utils import applyProperMotion, ObservationMetaData
from lsst.sims.utils import radiansFromArcsec, arcsecFromRadians
from lsst.sims.utils import ModifiedJulianDate
from lsst.sims.utils import _angularSeparation, angularSeparation
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

ROOT = os.path.abspath(os.path.dirname(__file__))


def setup_module(module):
    lsst.utils.tests.init()


class AlertDataGeneratorTestCase(unittest.TestCase):

    longMessage = True

    @classmethod
    def setUpClass(cls):
        print('setting up %s' % sims_clean_up.targets)
        cls.opsim_db = os.path.join(getPackageDir('sims_data'),
                                    'OpSimData',
                                    'opsimblitz1_1133_sqlite.db')

        rng = np.random.RandomState(8123)

        obs_gen = ObservationMetaDataGenerator(database=cls.opsim_db)
        cls.obs_list = obs_gen.getObservationMetaData(night=(0,2))
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
                          (simobjid int, ra real, dec real,
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

            subset = rng.randint(0,high=len(var_amp)-1, size=3)
            var_amp[subset[:2]] = 0.0
            var_amp[subset[-1]] = -1.0

            umag = rng.random_sample(n_stars)*5.0 + 15.0
            gmag = rng.random_sample(n_stars)*5.0 + 15.0
            rmag = rng.random_sample(n_stars)*5.0 + 15.0
            imag = rng.random_sample(n_stars)*5.0 + 15.0
            zmag = rng.random_sample(n_stars)*5.0 + 15.0
            ymag = rng.random_sample(n_stars)*5.0 + 15.0
            px = rng.random_sample(n_stars)*0.1 # say it is arcsec
            pmra = rng.random_sample(n_stars)*50.0+100.0 # say it is arcsec/yr
            pmdec = rng.random_sample(n_stars)*50.0+100.0 # say it is arcsec/yr
            vrad = rng.random_sample(n_stars)*600.0 - 300.0

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

            cls.max_str_len = -1

            for i_star in range(n_stars):
                if var_amp[i_star] >=-0.1:
                    varParamStr = ('{"m":"alert_test", "p":{"amp":%.4f, "per": %.4f}}'
                                   % (var_amp[i_star], var_period[i_star]))
                else:
                    varParamStr = 'None'

                if len(varParamStr) > cls.max_str_len:
                    cls.max_str_len = len(varParamStr)

                query = ('''INSERT INTO stars VALUES(%d, %.6f, %.6f,
                                                    %.4f, %.4f, %.4f, %.4f, %.4f, %.4f,
                                                    %.4f, %.4f, %.4f, %.4f, '%s')'''
                         % (i_star+id_offset+1, ra[i_star], dec[i_star],
                            umag[i_star], gmag[i_star], rmag[i_star],
                            imag[i_star], zmag[i_star], ymag[i_star],
                            px[i_star], pmra[i_star], pmdec[i_star],
                            vrad[i_star], varParamStr))

                cursor.execute(query)
        conn.commit()
        conn.close()

        cls.output_dir = tempfile.mkdtemp(dir=ROOT, prefix='alert_gen_output')
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
        for file_name in os.listdir(cls.output_dir):
            os.unlink(os.path.join(cls.output_dir, file_name))
        shutil.rmtree(cls.output_dir)
        attr_list = list(chipNameFromPupilCoordsLSST.__dict__.keys())
        for attr in attr_list:
            chipNameFromPupilCoordsLSST.__delattr__(attr)
        gc.collect()

    def test_alert_generation(self):

        _max_var_param_str = self.max_str_len

        class StarAlertTestDBObj(CatalogDBObject):
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
                       ('variabilityParameters', 'varParamStr', str, _max_var_param_str)]


        class TestAlertsVarCatMixin(object):

            @register_method('alert_test')
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

        class TestAlertsVarCat(TestAlertsVarCatMixin, AlertStellarVariabilityCatalog):
            pass


        star_db = StarAlertTestDBObj(database=self.star_db_name, driver='sqlite')

        dmag_cutoff = 0.005
        output_root = os.path.join(self.output_dir, 'alert_test')
        alert_gen = AlertDataGenerator(n_proc_max=1,
                                       photometry_class=TestAlertsVarCat,
                                       output_prefix = output_root,
                                       dmag_cutoff=dmag_cutoff,
                                       testing=True)

        alert_gen.subdivide_obs(self.obs_list)

        for htmid in alert_gen.htmid_list:
            alert_gen.alert_data_from_htmid(htmid, star_db)

        dummy_sed = Sed()

        template_m5_dict = {'u': 23.9, 'g': 25.0, 'r': 24.7, 'i': 24.0,
                            'z': 23.3, 'y': 22.1}  # from Table 2 of the overview paper

        bp_dict = BandpassDict.loadTotalBandpassesFromFiles()

        photParams = PhotometricParameters()

        # First, verify that the contents of the hdf5 files are all correct

        # While doing that, keep track of the uniqueId of every event that
        # is simulated for each obsHistID.  Afteward, we will verify that we
        # got all of the objects that satisfy dmag_cutoff
        all_simulated_events_dict = {}

        hdf_file_list = os.listdir(self.output_dir)
        for file_name in hdf_file_list:
            if not file_name.endswith('hdf5'):
                continue
            full_name = os.path.join(self.output_dir, file_name)
            test_file = h5py.File(full_name, 'r')
            mjd_list = ModifiedJulianDate.get_list(TAI=test_file['TAI'].value)
            band_list = test_file['bandpass'].value
            obshistID_list = test_file['obshistID'].value
            for obshistID, mjd, bandpass in zip(obshistID_list, mjd_list, band_list):
                all_simulated_events_dict[obshistID] = []

                # get the current ObservationMetaData
                for obs in self.obs_list:
                    if obs.OpsimMetaData['obsHistID'] == obshistID:
                        current_obs = obs
                        break

                ct_list = test_file['%d_map' % obshistID].value
                self.assertGreater(len(ct_list), 0)
                for batch_ct in ct_list:
                    id_list = test_file['%d_%d_uniqueId' % (obshistID, batch_ct)].value
                    for id_val in id_list:
                        all_simulated_events_dict[obshistID].append(id_val)
                    ra_list = test_file['%d_%d_raICRS' % (obshistID, batch_ct)].value
                    dec_list = test_file['%d_%d_decICRS' % (obshistID, batch_ct)].value
                    flux_list = test_file['%d_%d_flux' % (obshistID, batch_ct)].value
                    dflux_list = test_file['%d_%d_dflux' % (obshistID, batch_ct)].value
                    chipnum_list = test_file['%d_%d_chipNum' % (obshistID, batch_ct)].value
                    xpix_list = test_file['%d_%d_xPix' % (obshistID, batch_ct)].value
                    ypix_list = test_file['%d_%d_yPix' % (obshistID, batch_ct)].value
                    snr_list = test_file['%d_%d_SNR' % (obshistID, batch_ct)].value
                    self.assertGreater(len(id_list), 0)
                    self.assertEqual(len(id_list), len(ra_list))
                    self.assertEqual(len(id_list), len(dec_list))
                    for i_obj in range(len(id_list)):
                        obj_dex = (id_list[i_obj]//1024)-1

                        # verify that ICRS positions are correct to within 0.001 arcsec
                        ra0 = self.ra_truth[obj_dex]
                        dec0 = self.dec_truth[obj_dex]
                        px = self.px_truth[obj_dex]
                        pmra = self.pmra_truth[obj_dex]
                        pmdec = self.pmdec_truth[obj_dex]
                        vrad = self.vrad_truth[obj_dex]
                        raICRS, decICRS = applyProperMotion(ra0, dec0,
                                                            pmra,pmdec,px,
                                                            vrad, mjd=mjd)

                        dd = _angularSeparation(ra_list[i_obj], dec_list[i_obj],
                                                np.radians(raICRS), np.radians(decICRS))

                        dd_moved = angularSeparation(ra0, dec0, raICRS, decICRS)*3600.0

                        msg = '\nPosition (hdf5): %e %e\n' % (ra_list[i_obj], dec_list[i_obj])
                        msg += 'Position (truth): %e %e\n' % (np.radians(raICRS), np.radians(decICRS))
                        msg += 'diff %e arcsec; moved %e arsec\n' % (arcsecFromRadians(dd), dd_moved)
                        msg += 'pmra %e pmdec %e px %e vrad %e\n' % (pmra, pmdec, px, vrad)

                        self.assertLess(arcsecFromRadians(dd), 0.001, msg=msg)

                        # verify that flux calculations are correct
                        amp = self.amp_truth[obj_dex]
                        period = self.period_truth[obj_dex]
                        mag0 = self.mag0_truth_dict[bandpass][obj_dex]
                        dmag = amp*np.cos(period*mjd.TAI)
                        flux = dummy_sed.fluxFromMag(mag0+dmag)
                        dflux = flux - dummy_sed.fluxFromMag(mag0)
                        msg = ('\nuniqueID %d TAI %.4f\ndFlux (hdf5): %e\ndFlux (truth): %e\ndmag (truth): %e\n' %
                               (id_list[i_obj], mjd.TAI, dflux_list[i_obj], dflux, dmag))
                        msg+="amp %e period %e\n" % (amp, period)
                        msg+="obsHistID %d\n" % obshistID
                        msg +="flux %.6f %.6f\n" % (flux, flux_list[i_obj])
                        msg += "mag %e\n" % mag0

                        self.assertGreater(np.abs(dmag), dmag_cutoff)
                        self.assertAlmostEqual(dflux/dflux_list[i_obj], 1.0, 1, msg=msg)

                        msg = ('\nuniqueID %d TAI %.4f\nFlux (hdf5): %e\nFlux (truth): %e\ndmag (truth): %e\n' %
                               (id_list[i_obj], mjd.TAI, flux_list[i_obj], flux, dmag))
                        msg+="amp %e period %e\n" % (amp, period)

                        self.assertAlmostEqual(flux/flux_list[i_obj], 1.0, 1, msg=msg)

                        # verify that chipNum and pixel positions are correct
                        chipname = chipNameFromRaDecLSST(ra0, dec0, pm_ra=pmra, pm_dec=pmdec,
                                                         parallax=px, v_rad=vrad, obs_metadata=current_obs)

                        chipnum = int(chipname.replace('R','').replace('S','').replace(':','').
                                      replace(',','').replace(' ',''))
                        self.assertEqual(chipnum, chipnum_list[i_obj])

                        xpix, ypix = pixelCoordsFromRaDecLSST(ra0, dec0, pm_ra=pmra, pm_dec=pmdec,
                                                              parallax=px, v_rad=vrad,
                                                              obs_metadata=current_obs)

                        msg = '\nPixel position (hdf5): %.6f %.6f\n' % (xpix_list[i_obj], ypix_list[i_obj])
                        msg += 'Pixel position (true): %.6f %.6f\n' % (xpix, ypix)
                        self.assertAlmostEqual(xpix, xpix_list[i_obj], 4, msg=msg)
                        self.assertAlmostEqual(ypix, ypix_list[i_obj], 4, msg=msg)

                        # verify that signal to noise calculation is correct
                        bp = bp_dict[current_obs.bandpass]
                        q_m5 = template_m5_dict[current_obs.bandpass]
                        m5 = current_obs.m5[current_obs.bandpass]

                        q_snr, gamma = calcSNR_m5(mag0, bp, q_m5, photParams)

                        obs_snr, gamma = calcSNR_m5(dummy_sed.magFromFlux(flux), bp, m5, photParams)

                        q_sigma = dummy_sed.fluxFromMag(mag0)/q_snr
                        obs_sigma = flux/obs_snr

                        sigma = np.sqrt(q_sigma**2 + obs_sigma**2)
                        snr = dflux/sigma

                        self.assertAlmostEqual(snr_list[i_obj]/snr, 1.0, 4)

        # now verify that we simulated all of the events we were supposed to
        class TestAlertsTruthCat(TestAlertsVarCatMixin, CameraCoordsLSST, AstrometryStars,
                                 Variability, InstanceCatalog):
            column_outputs = ['uniqueId', 'chipName', 'dmagAlert']

            @compound('delta_umag', 'delta_gmag', 'delta_rmag',
                      'delta_imag', 'delta_zmag', 'delta_ymag')
            def get_TruthVariability(self):
                return self.applyVariability(self.column_by_name('varParamStr'))

            @cached
            def get_dmagAlert(self):
                return self.column_by_name('delta_%smag' % self.obs_metadata.bandpass)

        for obs in self.obs_list:
            # make a much larger ObservationMetaData to make sure we get all of
            # the objects
            obs.boundLength = 10.0

            obshistid = obs.OpsimMetaData['obsHistID']
            cat = TestAlertsTruthCat(star_db, obs_metadata=obs)
            id_list = []
            for line in cat.iter_catalog():
                id_val = line[0]
                chip_name = line[1]
                dmag = line[2]
                if chip_name is not None and np.abs(dmag)>=dmag_cutoff:
                    id_list.append(id_val)
                    msg = '\nchipName: %s\n' % chip_name
                    msg += 'dmag: %e\n' % dmag
                    msg += 'obshistid: %d\n' % obshistid
                    self.assertIn(id_val, all_simulated_events_dict[obshistid], msg=msg)

            for id_val in all_simulated_events_dict[obshistid]:
                self.assertIn(id_val, id_list)

        del alert_gen
        del cat
        gc.collect()


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
