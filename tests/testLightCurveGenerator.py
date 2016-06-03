from __future__ import with_statement
import unittest
import os
import json
import numpy as np
import lsst.utils.tests as utilsTests
from lsst.utils import getPackageDir

from lsst.sims.catalogs.generation.db import fileDBObject

from lsst.sims.catalogs.measures.instance import InstanceCatalog, compound
from lsst.sims.catUtils.mixins import PhotometryStars, VariabilityStars

from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.catUtils.utils import StellarLightCurveGenerator

class stellarControlCatalog(InstanceCatalog,
                            PhotometryStars, VariabilityStars):

    column_outputs = ["uniqueId", "raJ2000", "decJ2000", "mag", "sigma_mag"]

    @compound("mag", "sigma_mag")
    def get_phot(self):
        mm = self.column_by_name("lsst_%s" % self.obs_metadata.bandpass)
        sig = self.column_by_name("sigma_lsst_%s" % self.obs_metadata.bandpass)
        return np.array([mm, sig])


class StellarLightCurveTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        """
        Create a fake catalog of RR Lyrae stars.  Store it in cls.stellar_db
        """
        cls.scratchDir = os.path.join(getPackageDir("sims_catUtils"))
        cls.scratchDir = os.path.join(cls.scratchDir, "tests", "scratchSpace")

        rng = np.random.RandomState(88)
        n_stars = 10000
        sed_dir = os.path.join(getPackageDir("sims_sed_library"))
        sed_dir = os.path.join(sed_dir, "starSED", "kurucz")
        list_of_seds = os.listdir(sed_dir)

        lc_dir = os.path.join(getPackageDir("sims_sed_library"), "rrly_lc")
        lc_dir = os.path.join(lc_dir, "RRab")
        list_of_lc = ['rrly_lc/RRab/%s' % ww for ww in os.listdir(lc_dir) if "per.txt" in ww]

        cls.dtype = np.dtype([('id', np.int),
                             ('raDeg', np.float),
                             ('decDeg', np.float),
                             ('raJ2000', np.float),
                             ('decJ2000', np.float),
                             ('magNorm', np.float),
                             ('galacticAv', np.float),
                             ('sedFilename', str, 300),
                             ('varParamStr', str, 300)])

        # write the catalog as a text file to be ingested with fileDBObject
        cls.txt_name = os.path.join(cls.scratchDir, "stellar_lc_catalog.txt")
        with open(cls.txt_name, "w") as output_file:
            sed_dex = rng.random_integers(0, len(list_of_seds)-1, size=n_stars)
            lc_dex = rng.random_integers(0, len(list_of_lc)-1, size=n_stars)
            mjd0 = rng.random_sample(n_stars)*10000.0+40000.0
            raList = rng.random_sample(n_stars)*360.0
            decList = -90.0 + rng.random_sample(n_stars)*120.0
            magNormList = rng.random_sample(n_stars)*3.0+14.0
            AvList = rng.random_sample(n_stars)*0.2+0.1
            for ix in range(n_stars):
                varparams = {'varMethodName':'applyRRly',
                             'pars':{'tStartMjd':mjd0[ix],
                                     'filename':list_of_lc[lc_dex[ix]]}}
                varparamstr = json.dumps(varparams)
                output_file.write("%d;%lf;%lf;%lf;%lf;%lf;%lf;%s;%s\n"
                                  % (ix, raList[ix], decList[ix],
                                     np.radians(raList[ix]),
                                     np.radians(decList[ix]),
                                     magNormList[ix], AvList[ix],
                                     list_of_seds[sed_dex[ix]],
                                     varparamstr))

        cls.stellar_db = fileDBObject(cls.txt_name, delimiter=';',
                                      runtable='test', dtype=cls.dtype,
                                      idColKey='id')

        cls.stellar_db.raColName='raDeg'
        cls.stellar_db.decColName='decDeg'
        cls.stellar_db.objectTypeId=32

        cls.opsimDb = os.path.join(getPackageDir("sims_data"), "OpSimData")
        cls.opsimDb = os.path.join(cls.opsimDb, "opsimblitz1_1133_sqlite.db")


    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.txt_name):
            os.unlink(cls.txt_name)


    def test_stellar_light_curves(self):
        """
        Test the StellarLightCurveGenerator by generating some RR Lyrae light
        curves and comparing them to the results obtained by generating a
        series of InstanceCatalogs containing the same objects at the same
        MJDs
        """

        raRange = (78.0, 85.0)
        decRange = (-69.0, -65.0)
        bandpass = 'r'


        lc_gen = StellarLightCurveGenerator(self.stellar_db, self.opsimDb)
        test_light_curves = lc_gen.generate_light_curves(raRange, decRange, bandpass)

        self.assertGreater(len(test_light_curves), 2) # make sure we got some light curves

        for unique_id in test_light_curves:
            # verify that the sources returned all do vary by making sure that the
            # np.diff run on the magnitudes reutrns something non-zero
            self.assertGreater(np.abs(np.diff(test_light_curves[unique_id][1])).max(), 0.0)
            self.assertGreater(len(test_light_curves[unique_id][0]), 0)

        # Now test that specifying a small chunk_size does not change the output
        # light curves
        chunk_light_curves = lc_gen.generate_light_curves(raRange, decRange, bandpass, chunk_size=1)

        for unique_id in test_light_curves:
            self.assertEqual(len(test_light_curves[unique_id][0]), len(chunk_light_curves[unique_id][0]))
            np.testing.assert_array_equal(test_light_curves[unique_id][0], chunk_light_curves[unique_id][0])
            np.testing.assert_array_equal(test_light_curves[unique_id][1], chunk_light_curves[unique_id][1])
            np.testing.assert_array_equal(test_light_curves[unique_id][2], chunk_light_curves[unique_id][2])

        # Now find all of the ObservationMetaData that were included in our
        # light curves, generate InstanceCatalogs from them separately,
        # and verify that the contents of the InstanceCatalogs agree with
        # the contents of the light curves.

        gen = ObservationMetaDataGenerator(database=self.opsimDb,
                                           driver='sqlite')

        obs_list = gen.getObservationMetaData(fieldRA=raRange,
                                              fieldDec=decRange,
                                              telescopeFilter=bandpass,
                                              boundLength=1.75)

        for obs in obs_list:
            cat = stellarControlCatalog(self.stellar_db,
                                        obs_metadata=obs)

            for star_obj in cat.iter_catalog():
                lc = test_light_curves[star_obj[0]]
                dex = np.argmin(np.abs(lc[0]-obs.mjd.TAI))
                self.assertLess(np.abs(lc[0][dex]-obs.mjd.TAI), 1.0e-7)
                self.assertLess(np.abs(lc[1][dex]-star_obj[3]), 1.0e-7)
                self.assertLess(np.abs(lc[2][dex]-star_obj[4]), 1.0e-7)


    def test_date_range(self):
        """
        Run test_stellar_light_curves, this time specifying a range in MJD.
        """

        raRange = (0.0, 110.0)
        decRange = (-90.0, -50.0)
        bandpass = 'g'
        mjdRange = (49356.0, 49357.0)

        lc_gen = StellarLightCurveGenerator(self.stellar_db, self.opsimDb)
        test_light_curves = lc_gen.generate_light_curves(raRange, decRange,
                                                         bandpass,
                                                         expMJD=mjdRange)

        self.assertGreater(len(test_light_curves), 2)

        for unique_id in test_light_curves:
            # verify that the sources returned all do vary by making sure that the
            # np.diff run on the magnitudes reutrns something non-zero
            self.assertGreater(np.abs(np.diff(test_light_curves[unique_id][1])).max(), 0.0)
            self.assertGreater(len(test_light_curves[unique_id][0]), 0)
            self.assertGreater(test_light_curves[unique_id][0].min(), mjdRange[0]-1.0e-12)
            self.assertLess(test_light_curves[unique_id][0].max(), mjdRange[1]+1.0e-12)

        # Now test that specifying a small chunk_size does not change the output
        # light curves
        chunk_light_curves = lc_gen.generate_light_curves(raRange, decRange, bandpass,
                                                          expMJD=mjdRange, chunk_size=1)

        for unique_id in test_light_curves:
            self.assertEqual(len(test_light_curves[unique_id][0]), len(chunk_light_curves[unique_id][0]))
            np.testing.assert_array_equal(test_light_curves[unique_id][0], chunk_light_curves[unique_id][0])
            np.testing.assert_array_equal(test_light_curves[unique_id][1], chunk_light_curves[unique_id][1])
            np.testing.assert_array_equal(test_light_curves[unique_id][2], chunk_light_curves[unique_id][2])

        # Now find all of the ObservationMetaData that were included in our
        # light curves, generate InstanceCatalogs from them separately,
        # and verify that the contents of the InstanceCatalogs agree with
        # the contents of the light curves.

        gen = ObservationMetaDataGenerator(database=self.opsimDb,
                                           driver='sqlite')

        obs_list = gen.getObservationMetaData(fieldRA=raRange,
                                              fieldDec=decRange,
                                              telescopeFilter=bandpass,
                                              expMJD=mjdRange,
                                              boundLength=1.75)


        for obs in obs_list:
            cat = stellarControlCatalog(self.stellar_db,
                                        obs_metadata=obs)

            for star_obj in cat.iter_catalog():
                lc = test_light_curves[star_obj[0]]
                dex = np.argmin(np.abs(lc[0]-obs.mjd.TAI))
                self.assertLess(np.abs(lc[0][dex]-obs.mjd.TAI), 1.0e-7)
                self.assertLess(np.abs(lc[1][dex]-star_obj[3]), 1.0e-7)
                self.assertLess(np.abs(lc[2][dex]-star_obj[4]), 1.0e-7)


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(StellarLightCurveTest)

    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)
if __name__ == "__main__":
    run(True)
