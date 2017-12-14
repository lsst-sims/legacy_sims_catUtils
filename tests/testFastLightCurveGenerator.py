import unittest

import os
import numpy as np
import json
import tempfile
import shutil

import lsst.utils.tests
from lsst.utils import getPackageDir

from lsst.sims.catalogs.db import fileDBObject
from lsst.sims.catUtils.utils import FastStellarLightCurveGenerator
from lsst.sims.catUtils.utils import StellarLightCurveGenerator
from lsst.sims.catUtils.utils import AgnLightCurveGenerator
from lsst.sims.catUtils.utils import FastAgnLightCurveGenerator

from lsst.sims.utils.CodeUtilities import sims_clean_up
from lsst.sims.utils import ModifiedJulianDate

ROOT = os.path.abspath(os.path.dirname(__file__))


class FastStellar_stellar_lc_gen_case(unittest.TestCase):

    longMessage = True

    @classmethod
    def setUpClass(cls):
        """
        Create a fake catalog of RR Lyrae stars and MLT dwarves with flaring
        light curves. Store it in cls.stellar_db
        """
        cls.scratchDir = tempfile.mkdtemp(dir=ROOT, prefix='FastStellar_stellar_lc_gen_case-')

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
            sed_dex = rng.randint(0, len(list_of_seds), size=cls.n_stars//2)
            lc_dex = rng.randint(0, len(list_of_rrly_lc), size=cls.n_stars//2)
            mjd0 = rng.random_sample(cls.n_stars//2)*10000.0+40000.0
            raList = rng.random_sample(cls.n_stars//2)*(cls.raRange[1]-cls.raRange[0])+cls.raRange[0]
            decList = cls.decRange[0] + rng.random_sample(cls.n_stars//2)*(cls.decRange[1]-cls.decRange[1])
            magNormList = rng.random_sample(cls.n_stars//2)*3.0+14.0
            AvList = rng.random_sample(cls.n_stars//2)*0.2+0.1
            pxList = rng.random_sample(cls.n_stars//2)*0.1
            for ix in range(cls.n_stars//2):
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

            sed_dex = rng.randint(0, len(list_of_seds), size=cls.n_stars//2)
            lc_dex = rng.randint(1, 3, size=cls.n_stars//2)
            mjd0 = rng.random_sample(cls.n_stars//2)*10000.0+40000.0
            raList = rng.random_sample(cls.n_stars//2)*(cls.raRange[1]-cls.raRange[0])+cls.raRange[0]
            decList = cls.decRange[0] + rng.random_sample(cls.n_stars//2)*(cls.decRange[1]-cls.decRange[1])
            magNormList = rng.random_sample(cls.n_stars//2)*3.0+14.0
            AvList = rng.random_sample(cls.n_stars//2)*0.2+0.1
            pxList = rng.random_sample(cls.n_stars//2)*0.1
            for ix in range(cls.n_stars//2):
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
        if os.path.exists(cls.scratchDir):
            shutil.rmtree(cls.scratchDir)

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


class Fast_agn_lc_gen_test_case(unittest.TestCase):

    longMessge = True

    @classmethod
    def setUpClass(cls):
        rng = np.random.RandomState(119)

        n_galaxies = 20

        sed_dir = os.path.join(getPackageDir("sims_sed_library"), "galaxySED")
        list_of_seds = os.listdir(sed_dir)
        disk_sed_dexes = rng.randint(0, len(list_of_seds), size=n_galaxies)
        bulge_sed_dexes = rng.randint(0, len(list_of_seds), size=n_galaxies)

        avBulge = rng.random_sample(n_galaxies)*0.3+0.1
        avDisk = rng.random_sample(n_galaxies)*0.3+0.1

        mjdList = rng.random_sample(n_galaxies)*10.0+49330.0
        redshiftList = rng.random_sample(n_galaxies)*1.5+0.01

        tauList = rng.random_sample(n_galaxies)*1.0+1.0
        sfuList = rng.random_sample(n_galaxies)*2.0+1.0
        sfgList = rng.random_sample(n_galaxies)*2.0+1.0
        sfrList = rng.random_sample(n_galaxies)*2.0+1.0
        sfiList = rng.random_sample(n_galaxies)*2.0+1.0
        sfzList = rng.random_sample(n_galaxies)*2.0+1.0
        sfyList = rng.random_sample(n_galaxies)*2.0+1.0

        raList = rng.random_sample(n_galaxies)*7.0+78.0
        decList = rng.random_sample(n_galaxies)*4.0-69.0

        normDisk = rng.random_sample(n_galaxies)*5.0+20.0
        normBulge = rng.random_sample(n_galaxies)*5.0+20.0
        normAgn = rng.random_sample(n_galaxies)*5.0+20.0

        with lsst.utils.tests.getTempFilePath('.txt') as txt_cat_name:
            with open(txt_cat_name, "w") as output_file:
                for ix in range(n_galaxies):
                    varParam = {'varMethodName': 'applyAgn',
                                'pars': {'agn_tau': tauList[ix], 'agn_sfu': sfuList[ix],
                                         'agn_sfg': sfgList[ix], 'agn_sfr': sfrList[ix],
                                         'agn_sfi': sfiList[ix], 'agn_sfz': sfzList[ix],
                                         'agn_sfy': sfyList[ix], 't0_mjd': mjdList[ix],
                                         'seed': rng.randint(0, 200000)}}

                    paramStr = json.dumps(varParam)

                    output_file.write("%d;%f;%f;" % (ix, raList[ix], decList[ix])
                                      + "%f;%f;" % (np.radians(raList[ix]), np.radians(decList[ix]))
                                      + "%f;" % (redshiftList[ix])
                                      + "%s;%f;%f;" % (list_of_seds[disk_sed_dexes[ix]],
                                                       avDisk[ix], normDisk[ix])
                                      + "%s;%f;%f;" % (list_of_seds[bulge_sed_dexes[ix]],
                                                       avBulge[ix], normBulge[ix])
                                      + "agn.spec;%s;%f\n" % (paramStr, normAgn[ix]))

            dtype = np.dtype([
                             ('galid', np.int),
                             ('raDeg', np.float), ('decDeg', np.float),
                             ('raJ2000', np.float), ('decJ2000', np.float),
                             ('redshift', np.float),
                             ('sedFilenameDisk', str, 300), ('internalAvDisk', np.float),
                             ('magNormDisk', np.float),
                             ('sedFilenameBulge', str, 300), ('internalAvBulge', np.float),
                             ('magNormBulge', np.float),
                             ('sedFilenameAgn', str, 300), ('varParamStr', str, 600),
                             ('magNormAgn', np.float)
                             ])

            cls.agn_db = fileDBObject(txt_cat_name, delimiter=';',
                                      runtable='test', dtype=dtype,
                                      idColKey='galid')

        cls.agn_db.raColName = 'raDeg'
        cls.agn_db.decColName = 'decDeg'
        cls.agn_db.objectTypeId = 112

        # what follows is a hack to deal with the fact thar
        # our varParamStr values are longer than 256 characters
        # which is the default maximum length that a
        # CatalogDBObject expects a string to be
        #
        cls.agn_db.dbTypeMap['STRING'] = (str, 600)
        cls.agn_db.columns = None
        cls.agn_db._make_default_columns()
        cls.agn_db._make_column_map()
        cls.agn_db._make_type_map()

        cls.opsimDb = os.path.join(getPackageDir("sims_data"), "OpSimData")
        cls.opsimDb = os.path.join(cls.opsimDb, "opsimblitz1_1133_sqlite.db")

    @classmethod
    def tearDownClass(cls):
        sims_clean_up()

    def test_fast_agn_light_curves(self):
        raRange = (78.0, 85.0)
        decRange = (-69.0, -65.0)
        bandpass = ('g', 'r')

        slow_lc_gen = AgnLightCurveGenerator(self.agn_db, self.opsimDb)
        pointings = slow_lc_gen.get_pointings(raRange, decRange, bandpass=bandpass)
        for row in pointings:
            for obs in row:
                mjd = ModifiedJulianDate(TAI=obs.mjd.TAI-49000.0+59580.0)
                obs.mjd = mjd

        slow_lc, slow_truth = slow_lc_gen.light_curves_from_pointings(pointings)

        self.assertGreater(len(slow_lc), 2)  # make sure we got some light curves

        fast_lc_gen = FastAgnLightCurveGenerator(self.agn_db, self.opsimDb)
        pointings = fast_lc_gen.get_pointings(raRange, decRange, bandpass=bandpass)
        for row in pointings:
            for obs in row:
                mjd = ModifiedJulianDate(TAI=obs.mjd.TAI-49000.0+59580.0)
                obs.mjd = mjd

        fast_lc, fast_truth = fast_lc_gen.light_curves_from_pointings(pointings)

        self.assertEqual(len(slow_lc), len(fast_lc))
        self.assertEqual(len(slow_truth), len(fast_truth))

        for obj_id in slow_lc:
            self.assertEqual(len(fast_lc[obj_id]), len(slow_lc[obj_id]))
            for bp in fast_lc[obj_id]:
                self.assertEqual(len(fast_lc[obj_id][bp]), len(slow_lc[obj_id][bp]))
                for data_key in fast_lc[obj_id][bp]:
                    self.assertEqual(fast_lc[obj_id][bp][data_key].shape, slow_lc[obj_id][bp][data_key].shape)
                    self.assertLess(np.abs(fast_lc[obj_id][bp][data_key]-slow_lc[obj_id][bp][data_key]).max(),
                                    2.0e-10, msg='failed on %d, %s, %s' % (obj_id, bp, data_key))


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
