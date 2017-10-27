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
from lsst.sims.coordUtils import lsst_camera
from lsst.sims.coordUtils import chipNameFromPupilCoordsLSST
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.utils import AlertStellarVariabilityCatalog
from lsst.sims.catUtils.utils import AlertDataGenerator


ROOT = os.path.abspath(os.path.dirname(__file__))


def setup_module(module):
    lsst.utils.tests.init()


class AlertDataGeneratorTestCase(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        print('setting up %s' % sims_clean_up.targets)
        cls.opsim_db = os.path.join(getPackageDir('sims_data'),
                                    'OpSimData',
                                    'opsimblitz1_1133_sqlite.db')

        rng = np.random.RandomState(8123)

        obs_gen = ObservationMetaDataGenerator(database=cls.opsim_db)
        cls.obs_list = obs_gen.getObservationMetaData(night=(0,2))
        cls.obs_list = rng.choice(cls.obs_list, 10)
        fieldid_list = []
        for obs in cls.obs_list:
            fieldid_list.append(obs.OpsimMetaData['fieldID'])

        # make sure we have selected observations such that the
        # same field is revisited more than once
        assert len(np.unique(fieldid_list)) < len(fieldid_list)

        cls.temp_dir = tempfile.mkdtemp(prefix='alertDataGen',
                                        dir=ROOT)

        cls.star_db_name = tempfile.mktemp(prefix='alertDataGen_star_db',
                                           dir=cls.temp_dir,
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
        cls.u_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)
        cls.g_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)
        cls.r_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)
        cls.i_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)
        cls.z_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)
        cls.y_truth = np.zeros(n_stars*len(cls.obs_list), dtype=float)
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
            pmra = rng.random_sample(n_stars)*0.5 # say it is arcsec/yr
            pmdec = rng.random_sample(n_stars)*0.5 # say it is arcsec/yr
            vrad = rng.random_sample(n_stars)*600.0 - 300.0

            cls.ra_truth[id_offset:id_offset+n_stars] = ra
            cls.dec_truth[id_offset:id_offset+n_stars] = dec
            cls.u_truth[id_offset:id_offset+n_stars] = umag
            cls.g_truth[id_offset:id_offset+n_stars] = gmag
            cls.r_truth[id_offset:id_offset+n_stars] = rmag
            cls.i_truth[id_offset:id_offset+n_stars] = imag
            cls.z_truth[id_offset:id_offset+n_stars] = zmag
            cls.y_truth[id_offset:id_offset+n_stars] = ymag
            cls.px_truth[id_offset:id_offset+n_stars] = px
            cls.pmra_truth[id_offset:id_offset+n_stars] = pmra
            cls.pmdec_truth[id_offset:id_offset+n_stars] = pmdec
            cls.vrad_truth[id_offset:id_offset+n_stars] = vrad
            cls.amp_truth[id_offset:id_offset+n_stars] = var_amp
            cls.period_truth[id_offset:id_offset+n_stars] = var_period

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


    @classmethod
    def tearDownClass(cls):
        sims_clean_up()
        if os.path.exists(cls.star_db_name):
            os.unlink(cls.star_db_name)
        if os.path.exists(cls.temp_dir):
            shutil.rmtree(cls.temp_dir)
        for file_name in os.listdir(cls.output_dir):
            os.unlink(os.path.join(cls.output_dir, file_name))
        shutil.rmtree(cls.output_dir)
        del lsst_camera._lsst_camera
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


        class TestAlertVarCat(AlertStellarVariabilityCatalog):

            @register_method('alert_test')
            def applyAlertTest(self, valid_dexes, params, expmjd, variability_cache=None):
                if len(params) == 0:
                    return np.array([[], [], [], [], [], []])

                if isinstance(expmjd, numbers.Number):
                    raise RuntimeError("did not mean for this method to operate "
                                       "on a single MJD")

                dMags_out = np.zeros((6, self.num_variable_obj(params), len(expmjd)))

                for i_star in range(self.num_variable_obj(params)):
                    if params['amp'][i_star] is not None:
                        dmags = params['amp'][i_star]*np.cos(params['per'][i_star]*expmjd)
                        for i_filter in range(6):
                            dMags_out[i_filter][i_star] = dmags

                return dMags_out


        star_db = StarAlertTestDBObj(database=self.star_db_name, driver='sqlite')

        output_root = os.path.join(self.output_dir, 'alert_test')
        alert_gen = AlertDataGenerator(n_proc_max=1,
                                       photometry_class=TestAlertVarCat,
                                       output_prefix = output_root)

        alert_gen.subdivide_obs(self.obs_list)

        for htmid in alert_gen.htmid_list:
            alert_gen.alert_data_from_htmid(htmid, star_db)

        del alert_gen
        gc.collect()


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
