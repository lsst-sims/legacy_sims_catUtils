from __future__ import with_statement
import unittest
import os
import numpy as np
import sqlite3
from astropy.analytic_functions import blackbody_lambda
import lsst.utils.tests
from lsst.utils import getPackageDir
from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.catUtils.mixins import VariabilityStars
from lsst.sims.catUtils.mixins import MLTflaringMixin
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catUtils.mixins import PhotometryStars
from lsst.sims.utils import ObservationMetaData
from lsst.sims.photUtils import SedList, BandpassDict, Sed
from lsst.sims.utils import radiansFromArcsec
from lsst.sims.photUtils import PhotometricParameters

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

class MLT_flare_test_case(unittest.TestCase):

    longMessage = True

    @classmethod
    def setUpClass(cls):
        cls.scratch_dir = os.path.join(getPackageDir('sims_catUtils'),
                                       'tests', 'scratchSpace')
        # Create a dummy mlt light curve file
        cls.mlt_lc_name = os.path.join(cls.scratch_dir,
                                       'test_mlt_lc_file.npz')

        lc_files = {}
        amp = 1.0e32
        lc_files['lc_1_time'] = np.arange(0.0, 3652.51, 0.1)
        lc_files['lc_1_u'] = amp*(1.0+np.power(np.sin(lc_files['lc_1_time']/100.0), 2))
        lc_files['lc_1_g'] = amp*(1.0+np.power(np.cos(lc_files['lc_1_time']/100.0), 2))

        amp = 2.0e31
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
        if os.path.exists(cls.mlt_lc_name):
            os.unlink(cls.mlt_lc_name)

        if os.path.exists(cls.db_name):
            os.unlink(cls.db_name)

    def test_flare_magnitudes(self):
        """
        Test that we get the expected magnitudes out
        """
        db = MLT_test_DB(database=self.db_name, driver='sqlite')

        class QuiescentCatalog(PhotometryStars, InstanceCatalog):
            column_outputs = ['id', 'lsst_u', 'lsst_g']

        class FlaringCatalog(PhotometryStars, VariabilityStars,
                             MLTflaringMixin, InstanceCatalog):
            column_outputs = ['id', 'lsst_u', 'lsst_g']


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
        distance_list *= 3.0857e16  # convert to cm

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
                    amp = 1.0e32
                    dt = 3652.5
                    t_min = flare_cat._survey_start - t0_list[obj_id]

                    tt = mjd - t_min
                    while tt > dt:
                        tt -= dt

                    u_flux = amp*(1.0+np.power(np.sin(tt/100.0), 2))
                    g_flux = amp*(1.0+np.power(np.cos(tt/100.0), 2))
                else:
                    amp = 2.0e31
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


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
