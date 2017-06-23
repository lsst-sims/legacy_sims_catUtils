import unittest

import os
import numpy as np
import json

import lsst.utils.tests
from lsst.utils import getPackageDir

from lsst.sims.catalogs.db import fileDBObject
from lsst.sims.catUtils.utils import FastStellarLightCurveGenerator
from lsst.sims.catUtils.utils import StellarLightCurveGenerator

from lsst.sims.utils.CodeUtilities import sims_clean_up

class FastStellar_lc_gen_case(unittest.TestCase):

    longMessage = True

    @classmethod
    def setUpClass(cls):
        """
        Create a fake catalog of RR Lyrae stars and MLT dwarves with flaring
        light curves. Store it in cls.stellar_db
        """
        cls.scratchDir = os.path.join(getPackageDir("sims_catUtils"))
        cls.scratchDir = os.path.join(cls.scratchDir, "tests", "scratchSpace")

        cls.raRange = (78.0, 85.0)
        cls.decRange = (-69.0, -65.0)

        rng = np.random.RandomState(88)
        cls.n_stars = 20
        sed_dir = os.path.join(getPackageDir("sims_sed_library"))
        sed_dir = os.path.join(sed_dir, "starSED", "kurucz")
        list_of_seds = os.listdir(sed_dir)

        lc_dir = os.path.join(getPackageDir("sims_sed_library"), "rrly_lc")
        lc_dir = os.path.join(lc_dir, "RRab")
        list_of_rrly_lc = ['rrly_lc/RRab/%s' % ww for ww in os.listdir(lc_dir) if "per.txt" in ww]

        cls.mlt_lc_file_name = os.path.join(cls.scratchDir, "fast_lc_mlt_file.npz")
        if os.path.exists(cls.mlt_lc_file_name):
            os.unlink(cls.mlt_lc_file_name)

        mlt_lc_files = {}
        mlt_lc_files['lc_1_time'] = np.arange(0.0, 3652.51, 0.1)
        mlt_lc_files['lc_1_g'] = 2.2e32*np.power(np.cos(mlt_lc_files['lc_1_time']/100.0-5.0),2)
        mlt_lc_files['lc_1_r'] = 1.3e32*(1.0+np.sin(mlt_lc_files['lc_1_time']/100.0-3.0))

        mlt_lc_files['lc_2_time'] = np.arange(0.0, 3652.51, 0.1)
        mlt_lc_files['lc_2_g'] = 5.1e33*(1.0+np.cos(mlt_lc_files['lc_2_time']/300.0-10.0))
        mlt_lc_files['lc_2_r'] = 4.3e32*(1.0+np.sin(mlt_lc_files['lc_2_time']/50.0-71.0))

        with open(cls.mlt_lc_file_name, 'wb') as file_handle:
            np.savez(file_handle, **mlt_lc_files)

        cls.dtype = np.dtype([('id', np.int),
                             ('raDeg', np.float),
                             ('decDeg', np.float),
                             ('raJ2000', np.float),
                             ('decJ2000', np.float),
                             ('magNorm', np.float),
                             ('galacticAv', np.float),
                             ('sedFilename', str, 300),
                             ('varParamStr', str, 300),
                             ('parallax', np.float),
                             ('ebv', np.float)])

        # write the catalog as a text file to be ingested with fileDBObject
        cls.txt_name = os.path.join(cls.scratchDir, "fast_stellar_lc_catalog.txt")
        with open(cls.txt_name, "w") as output_file:
            output_file.write('# a silly header\n')
            sed_dex = rng.random_integers(0, len(list_of_seds)-1, size=cls.n_stars/2)
            lc_dex = rng.random_integers(0, len(list_of_rrly_lc)-1, size=cls.n_stars/2)
            mjd0 = rng.random_sample(cls.n_stars/2)*10000.0+40000.0
            raList = rng.random_sample(cls.n_stars/2)*(cls.raRange[1]-cls.raRange[0])+cls.raRange[0]
            decList = cls.decRange[0] + rng.random_sample(cls.n_stars/2)*(cls.decRange[1]-cls.decRange[1])
            magNormList = rng.random_sample(cls.n_stars/2)*3.0+14.0
            AvList = rng.random_sample(cls.n_stars/2)*0.2+0.1
            pxList = rng.random_sample(cls.n_stars/2)*0.1
            for ix in range(cls.n_stars/2):
                varparams = {'varMethodName': 'applyRRly',
                             'pars': {'tStartMjd': mjd0[ix],
                                      'filename': list_of_rrly_lc[lc_dex[ix]]}}
                varparamstr = json.dumps(varparams)
                output_file.write("%d;%lf;%lf;%lf;%lf;%lf;%lf;%s;%s;%lf;%lf\n"
                                  % (ix, raList[ix], decList[ix],
                                     np.radians(raList[ix]),
                                     np.radians(decList[ix]),
                                     magNormList[ix], AvList[ix],
                                     list_of_seds[sed_dex[ix]],
                                     varparamstr,pxList[ix],
                                     AvList[ix]/3.1))

            sed_dex = rng.random_integers(0, len(list_of_seds)-1, size=cls.n_stars/2)
            lc_dex = rng.random_integers(1, 2, size=cls.n_stars/2)
            mjd0 = rng.random_sample(cls.n_stars/2)*10000.0+40000.0
            raList = rng.random_sample(cls.n_stars/2)*(cls.raRange[1]-cls.raRange[0])+cls.raRange[0]
            decList = cls.decRange[0] + rng.random_sample(cls.n_stars/2)*(cls.decRange[1]-cls.decRange[1])
            magNormList = rng.random_sample(cls.n_stars/2)*3.0+14.0
            AvList = rng.random_sample(cls.n_stars/2)*0.2+0.1
            pxList = rng.random_sample(cls.n_stars/2)*0.1
            for ix in range(cls.n_stars/2):
                varparams = {'m':'MLT', 'p':{'lc':'lc_%d' % lc_dex[ix], 't0':rng.random_sample()*1000.0}}
                varparamstr = json.dumps(varparams)
                output_file.write("%d;%lf;%lf;%lf;%lf;%lf;%lf;%s;%s;%lf;%lf\n"
                                  % (ix+cls.n_stars/2, raList[ix], decList[ix],
                                     np.radians(raList[ix]),
                                     np.radians(decList[ix]),
                                     magNormList[ix], AvList[ix],
                                     list_of_seds[sed_dex[ix]],
                                     varparamstr,pxList[ix],
                                     AvList[ix]/3.1))

        cls.stellar_db = fileDBObject(cls.txt_name, delimiter=';',
                                      runtable='test', dtype=cls.dtype,
                                      idColKey='id')

        cls.stellar_db.raColName = 'raDeg'
        cls.stellar_db.decColName = 'decDeg'
        cls.stellar_db.objectTypeId = 32

        cls.opsimDb = os.path.join(getPackageDir("sims_data"), "OpSimData")
        cls.opsimDb = os.path.join(cls.opsimDb, "opsimblitz1_1133_sqlite.db")

    @classmethod
    def tearDownClass(cls):
        sims_clean_up()
        if os.path.exists(cls.txt_name):
            os.unlink(cls.txt_name)
        if os.path.exists(cls.mlt_lc_file_name):
            os.unlink(cls.mlt_lc_file_name)

    def test_fast_stellar_lc_gen(self):

        bandpass = ('r', 'g')

        lc_slow = StellarLightCurveGenerator(self.stellar_db, self.opsimDb)
        lc_slow._lightCurveCatalogClass._mlt_lc_file = self.mlt_lc_file_name
        ptngs = lc_slow.get_pointings((68.0, 95.0), (-69.0, -55.0), bandpass=bandpass)
        slow_lc, slow_truth = lc_slow.light_curves_from_pointings(ptngs, chunk_size=10)
        self.assertEqual(len(slow_truth), self.n_stars)
        self.assertEqual(len(slow_lc), self.n_stars)

        lc_fast = FastStellarLightCurveGenerator(self.stellar_db, self.opsimDb)
        lc_fast._lightCurveCatalogClass._mlt_lc_file = self.mlt_lc_file_name
        ptngs = lc_fast.get_pointings((68.0, 95.0), (-69.0, -55.0), bandpass=bandpass)
        fast_lc, fast_truth = lc_fast.light_curves_from_pointings(ptngs, chunk_size=10)

        self.assertEqual(len(fast_lc), len(slow_lc))
        self.assertEqual(len(fast_truth), len(slow_truth))

        for obj_id in fast_lc:
            self.assertEqual(len(fast_lc[obj_id]), len(slow_lc[obj_id]))
            for bp in fast_lc[obj_id]:
                self.assertEqual(len(fast_lc[obj_id][bp]), len(slow_lc[obj_id][bp]))
                for data_key in fast_lc[obj_id][bp]:
                    self.assertEqual(len(slow_lc[obj_id][bp][data_key]), len(fast_lc[obj_id][bp][data_key]))
                    self.assertEqual(slow_lc[obj_id][bp][data_key].shape, slow_lc[obj_id][bp][data_key].shape)
                    self.assertLess(np.abs(fast_lc[obj_id][bp][data_key]-slow_lc[obj_id][bp][data_key]).max(),
                                    1.0e-10)

        for obj_id in fast_truth:
            self.assertEqual(fast_truth[obj_id], slow_truth[obj_id])


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
