from __future__ import with_statement
import unittest
import os
import json
import numpy as np
import lsst.utils.tests as utilsTests
from lsst.utils import getPackageDir

from lsst.sims.catalogs.generation.db import fileDBObject

from lsst.sims.catalogs.measures.instance import InstanceCatalog, compound
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator

from lsst.sims.catUtils.mixins import PhotometryStars, VariabilityStars
from lsst.sims.catUtils.utils import StellarLightCurveGenerator

from lsst.sims.catUtils.mixins import PhotometryGalaxies, VariabilityGalaxies
from lsst.sims.catUtils.utils import AgnLightCurveGenerator


class stellarControlCatalog(InstanceCatalog,
                            PhotometryStars, VariabilityStars):

    column_outputs = ["uniqueId", "raJ2000", "decJ2000", "mag", "sigma_mag"]

    @compound("mag", "sigma_mag")
    def get_phot(self):
        mm = self.column_by_name("lsst_%s" % self.obs_metadata.bandpass)
        sig = self.column_by_name("sigma_lsst_%s" % self.obs_metadata.bandpass)
        return np.array([mm, sig])


class agnControlCatalog(InstanceCatalog,
                        PhotometryGalaxies, VariabilityGalaxies):

    column_outputs = ["uniqueId", "raJ2000", "decJ2000", "mag", "sigma_mag"]

    @compound("mag", "sigma_mag")
    def get_phot(self):
        return np.array([
                         self.column_by_name("%sAgn" % self.obs_metadata.bandpass),
                         self.column_by_name("sigma_%sAgn" % self.obs_metadata.bandpass)
                        ])


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
                varparams = {'varMethodName': 'applyRRly',
                             'pars': {'tStartMjd': mjd0[ix],
                                      'filename': list_of_lc[lc_dex[ix]]}}
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

        cls.stellar_db.raColName = 'raDeg'
        cls.stellar_db.decColName = 'decDeg'
        cls.stellar_db.objectTypeId = 32

        cls.opsimDb = os.path.join(getPackageDir("sims_data"), "OpSimData")
        cls.opsimDb = os.path.join(cls.opsimDb, "opsimblitz1_1133_sqlite.db")

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.txt_name):
            os.unlink(cls.txt_name)

    def test_get_pointings(self):
        """
        Test that the get_pointings method does, in fact, return ObservationMetaData
        that are grouped appropriately.
        """

        raRange = (78.0, 89.0)
        decRange = (-74.0, -60.0)
        bandpass = 'g'

        lc_gen = StellarLightCurveGenerator(self.stellar_db, self.opsimDb)

        pointings = lc_gen.get_pointings(raRange, decRange, bandpass=bandpass)

        self.assertGreater(len(pointings), 1)

        for group in pointings:
            for ix, obs in enumerate(group):
                self.assertAlmostEqual(obs.pointingRA, group[0].pointingRA, 12)
                self.assertAlmostEqual(obs.pointingDec, group[0].pointingDec, 12)
                self.assertEqual(obs.bandpass, bandpass)
                if ix > 0:
                    self.assertGreater(obs.mjd.TAI, group[ix-1].mjd.TAI)

    def test_get_pointings_multiband(self):
        """
        Test that the get_pointings method does, in fact, return ObservationMetaData
        that are grouped appropriately.  Test on more than one filter.
        """

        raRange = (78.0, 89.0)
        decRange = (-74.0, -60.0)
        bandpass = ('g', 'z')

        lc_gen = StellarLightCurveGenerator(self.stellar_db, self.opsimDb)

        pointings = lc_gen.get_pointings(raRange, decRange, bandpass=bandpass)

        self.assertGreater(len(pointings), 1)

        ct_g = 0
        ct_z = 0

        for group in pointings:
            for ix, obs in enumerate(group):
                self.assertAlmostEqual(obs.pointingRA, group[0].pointingRA, 12)
                self.assertAlmostEqual(obs.pointingDec, group[0].pointingDec, 12)

                if obs.bandpass == 'g':
                    ct_g += 1
                elif obs.bandpass == 'z':
                    ct_z += 1
                else:
                    raise RuntimeError("Asked for filters (g,z) but got %s" % obs.bandpass)

                if ix > 0:
                    self.assertGreater(obs.mjd.TAI, group[ix-1].mjd.TAI)

        self.assertGreater(ct_g, 0)
        self.assertGreater(ct_z, 0)

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
        pointings = lc_gen.get_pointings(raRange, decRange, bandpass=bandpass)
        test_light_curves, truth_info = lc_gen.light_curves_from_pointings(pointings)

        self.assertGreater(len(test_light_curves), 2)  # make sure we got some light curves

        for unique_id in test_light_curves:
            # verify that the sources returned all do vary by making sure that the
            # np.diff run on the magnitudes reutrns something non-zero
            self.assertGreater(np.abs(np.diff(test_light_curves[unique_id][bandpass]['mag'])).max(), 0.0)
            self.assertGreater(len(test_light_curves[unique_id][bandpass]['mjd']), 0)

        # Now test that specifying a small chunk_size does not change the output
        # light curves
        chunk_light_curves, truth_info = lc_gen.light_curves_from_pointings(pointings, chunk_size=1)

        for unique_id in test_light_curves:
            self.assertEqual(len(test_light_curves[unique_id][bandpass]['mjd']),
                             len(chunk_light_curves[unique_id][bandpass]['mjd']))
            np.testing.assert_array_equal(test_light_curves[unique_id][bandpass]['mjd'],
                                          chunk_light_curves[unique_id][bandpass]['mjd'])
            np.testing.assert_array_equal(test_light_curves[unique_id][bandpass]['mag'],
                                          chunk_light_curves[unique_id][bandpass]['mag'])
            np.testing.assert_array_equal(test_light_curves[unique_id][bandpass]['error'],
                                          chunk_light_curves[unique_id][bandpass]['error'])

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
                lc = test_light_curves[star_obj[0]][bandpass]
                dex = np.argmin(np.abs(lc['mjd']-obs.mjd.TAI))
                self.assertLess(np.abs(lc['mjd'][dex]-obs.mjd.TAI), 1.0e-7)
                self.assertLess(np.abs(lc['mag'][dex]-star_obj[3]), 1.0e-7)
                self.assertLess(np.abs(lc['error'][dex]-star_obj[4]), 1.0e-7)

    def test_date_range(self):
        """
        Run test_stellar_light_curves, this time specifying a range in MJD.
        """

        raRange = (0.0, 110.0)
        decRange = (-90.0, -50.0)
        bandpass = 'g'
        mjdRange = (49356.0, 49357.0)

        lc_gen = StellarLightCurveGenerator(self.stellar_db, self.opsimDb)
        pointings = lc_gen.get_pointings(raRange, decRange, bandpass=bandpass, expMJD=mjdRange)
        test_light_curves, truth_info = lc_gen.light_curves_from_pointings(pointings)

        self.assertGreater(len(test_light_curves), 2)

        for unique_id in test_light_curves:
            # verify that the sources returned all do vary by making sure that the
            # np.diff run on the magnitudes reutrns something non-zero
            self.assertGreater(np.abs(np.diff(test_light_curves[unique_id][bandpass]['mag'])).max(), 0.0)
            self.assertGreater(len(test_light_curves[unique_id][bandpass]['mjd']), 0)
            self.assertGreater(test_light_curves[unique_id][bandpass]['mjd'].min(), mjdRange[0]-1.0e-12)
            self.assertLess(test_light_curves[unique_id][bandpass]['mjd'].max(), mjdRange[1]+1.0e-12)

        # Now test that specifying a small chunk_size does not change the output
        # light curves
        chunk_light_curves, truth_info = lc_gen.light_curves_from_pointings(pointings, chunk_size=1)

        for unique_id in test_light_curves:
            self.assertEqual(len(test_light_curves[unique_id][bandpass]['mjd']),
                             len(chunk_light_curves[unique_id][bandpass]['mjd']))
            np.testing.assert_array_equal(test_light_curves[unique_id][bandpass]['mjd'],
                                          chunk_light_curves[unique_id][bandpass]['mjd'])
            np.testing.assert_array_equal(test_light_curves[unique_id][bandpass]['mag'],
                                          chunk_light_curves[unique_id][bandpass]['mag'])
            np.testing.assert_array_equal(test_light_curves[unique_id][bandpass]['error'],
                                          chunk_light_curves[unique_id][bandpass]['error'])

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
                lc = test_light_curves[star_obj[0]][bandpass]
                dex = np.argmin(np.abs(lc['mjd']-obs.mjd.TAI))
                self.assertLess(np.abs(lc['mjd'][dex]-obs.mjd.TAI), 1.0e-7)
                self.assertLess(np.abs(lc['mag'][dex]-star_obj[3]), 1.0e-7)
                self.assertLess(np.abs(lc['error'][dex]-star_obj[4]), 1.0e-7)

    def test_multiband_light_curves(self):
        """
        Check that multi-band light curves are returned correctly.
        """

        raRange = (78.0, 82.0)
        decRange = (-69.0, -65.0)
        bandpass = ('r', 'g')

        gen = StellarLightCurveGenerator(self.stellar_db, self.opsimDb)
        pointings = gen.get_pointings(raRange, decRange, bandpass=bandpass)
        lc_dict, truth_info = gen.light_curves_from_pointings(pointings)

        obs_gen = ObservationMetaDataGenerator(database=self.opsimDb, driver='sqlite')
        control_pointings_r = obs_gen.getObservationMetaData(fieldRA=raRange, fieldDec=decRange,
                                                             telescopeFilter='r', boundLength=1.75)

        control_pointings_g = obs_gen.getObservationMetaData(fieldRA=raRange, fieldDec=decRange,
                                                             telescopeFilter='g', boundLength=1.75)

        self.assertGreater(len(control_pointings_g), 0)
        self.assertGreater(len(control_pointings_r), 0)

        for obs in control_pointings_r:
            cat = stellarControlCatalog(self.stellar_db,
                                        obs_metadata=obs)

            for star_obj in cat.iter_catalog():
                lc = lc_dict[star_obj[0]]['r']
                dex = np.argmin(np.abs(lc['mjd']-obs.mjd.TAI))
                self.assertLess(np.abs(lc['mjd'][dex]-obs.mjd.TAI), 1.0e-7)
                self.assertLess(np.abs(lc['mag'][dex]-star_obj[3]), 1.0e-7)
                self.assertLess(np.abs(lc['error'][dex]-star_obj[4]), 1.0e-7)

        for obs in control_pointings_g:
            cat = stellarControlCatalog(self.stellar_db,
                                        obs_metadata=obs)

            for star_obj in cat.iter_catalog():
                lc = lc_dict[star_obj[0]]['g']
                dex = np.argmin(np.abs(lc['mjd']-obs.mjd.TAI))
                self.assertLess(np.abs(lc['mjd'][dex]-obs.mjd.TAI), 1.0e-7)
                self.assertLess(np.abs(lc['mag'][dex]-star_obj[3]), 1.0e-7)
                self.assertLess(np.abs(lc['error'][dex]-star_obj[4]), 1.0e-7)


class AgnLightCurveTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        rng = np.random.RandomState(119)

        cls.txt_cat_name = os.path.join(getPackageDir("sims_catUtils"), "tests")
        cls.txt_cat_name = os.path.join(cls.txt_cat_name, "scratchSpace", "agn_lc_cat.txt")

        n_galaxies = 20

        sed_dir = os.path.join(getPackageDir("sims_sed_library"), "galaxySED")
        list_of_seds = os.listdir(sed_dir)
        disk_sed_dexes = rng.random_integers(0, len(list_of_seds)-1, size=n_galaxies)
        bulge_sed_dexes = rng.random_integers(0, len(list_of_seds)-1, size=n_galaxies)

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

        with open(cls.txt_cat_name, "w") as output_file:
            for ix in range(n_galaxies):
                varParam = {'varMethodName': 'applyAgn',
                            'pars': {'agn_tau': tauList[ix], 'agn_sfu': sfuList[ix],
                                     'agn_sfg': sfgList[ix], 'agn_sfr': sfrList[ix],
                                     'agn_sfi': sfiList[ix], 'agn_sfz': sfzList[ix],
                                     'agn_sfy': sfyList[ix], 't0_mjd': mjdList[ix],
                                     'seed': rng.randint(0, 200000)
                                    }}

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

        cls.agn_db = fileDBObject(cls.txt_cat_name, delimiter=';',
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
        if os.path.exists(cls.txt_cat_name):
            os.unlink(cls.txt_cat_name)

    def test_agn_light_curves(self):
        """
        Test the AgnLightCurveGenerator by generating some AGN light
        curves and comparing them to the results obtained by generating a
        series of InstanceCatalogs containing the same objects at the same
        MJDs
        """

        raRange = (78.0, 85.0)
        decRange = (-69.0, -65.0)
        bandpass = 'g'

        lc_gen = AgnLightCurveGenerator(self.agn_db, self.opsimDb)
        pointings = lc_gen.get_pointings(raRange, decRange, bandpass=bandpass)
        test_light_curves, truth_info = lc_gen.light_curves_from_pointings(pointings)

        self.assertGreater(len(test_light_curves), 2)  # make sure we got some light curves

        for unique_id in test_light_curves:
            # verify that the sources returned all do vary by making sure that the
            # np.diff run on the magnitudes reutrns something non-zero
            self.assertGreater(np.abs(np.diff(test_light_curves[unique_id][bandpass]['mag'])).max(), 0.0)
            self.assertGreater(len(test_light_curves[unique_id][bandpass]['mjd']), 0)

        # Now test that specifying a small chunk_size does not change the output
        # light curves
        chunk_light_curves, truth_info = lc_gen.light_curves_from_pointings(pointings, chunk_size=1)

        for unique_id in test_light_curves:
            self.assertEqual(len(test_light_curves[unique_id][bandpass]['mjd']),
                             len(chunk_light_curves[unique_id][bandpass]['mjd']))
            np.testing.assert_array_equal(test_light_curves[unique_id][bandpass]['mjd'],
                                          chunk_light_curves[unique_id][bandpass]['mjd'])
            np.testing.assert_array_equal(test_light_curves[unique_id][bandpass]['mag'],
                                          chunk_light_curves[unique_id][bandpass]['mag'])
            np.testing.assert_array_equal(test_light_curves[unique_id][bandpass]['error'],
                                          chunk_light_curves[unique_id][bandpass]['error'])

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
            cat = agnControlCatalog(self.agn_db,
                                    obs_metadata=obs)

            for agn_obj in cat.iter_catalog():
                lc = test_light_curves[agn_obj[0]][bandpass]
                dex = np.argmin(np.abs(lc['mjd']-obs.mjd.TAI))
                self.assertLess(np.abs(lc['mjd'][dex]-obs.mjd.TAI), 1.0e-7)
                self.assertLess(np.abs(lc['mag'][dex]-agn_obj[3]), 1.0e-7)
                self.assertLess(np.abs(lc['error'][dex]-agn_obj[4]), 1.0e-7)

    def test_multiband_light_curves(self):
        """
        Check that multi-band light curves are returned correctly.
        """

        raRange = (78.0, 82.0)
        decRange = (-69.0, -65.0)
        bandpass = ('r', 'g')

        gen = AgnLightCurveGenerator(self.agn_db, self.opsimDb)
        pointings = gen.get_pointings(raRange, decRange, bandpass=bandpass)
        lc_dict, truth_info = gen.light_curves_from_pointings(pointings)

        obs_gen = ObservationMetaDataGenerator(database=self.opsimDb, driver='sqlite')
        control_pointings_r = obs_gen.getObservationMetaData(fieldRA=raRange, fieldDec=decRange,
                                                             telescopeFilter='r', boundLength=1.75)

        control_pointings_g = obs_gen.getObservationMetaData(fieldRA=raRange, fieldDec=decRange,
                                                             telescopeFilter='g', boundLength=1.75)

        self.assertGreater(len(control_pointings_g), 0)
        self.assertGreater(len(control_pointings_r), 0)

        for obs in control_pointings_r:
            cat = agnControlCatalog(self.agn_db,
                                    obs_metadata=obs)

            for star_obj in cat.iter_catalog():
                lc = lc_dict[star_obj[0]]['r']
                dex = np.argmin(np.abs(lc['mjd']-obs.mjd.TAI))
                self.assertLess(np.abs(lc['mjd'][dex]-obs.mjd.TAI), 1.0e-7)
                self.assertLess(np.abs(lc['mag'][dex]-star_obj[3]), 1.0e-7)
                self.assertLess(np.abs(lc['error'][dex]-star_obj[4]), 1.0e-7)

        for obs in control_pointings_g:
            cat = agnControlCatalog(self.agn_db,
                                    obs_metadata=obs)

            for star_obj in cat.iter_catalog():
                lc = lc_dict[star_obj[0]]['g']
                dex = np.argmin(np.abs(lc['mjd']-obs.mjd.TAI))
                self.assertLess(np.abs(lc['mjd'][dex]-obs.mjd.TAI), 1.0e-7)
                self.assertLess(np.abs(lc['mag'][dex]-star_obj[3]), 1.0e-7)
                self.assertLess(np.abs(lc['error'][dex]-star_obj[4]), 1.0e-7)


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(StellarLightCurveTest)
    suites += unittest.makeSuite(AgnLightCurveTest)

    return unittest.TestSuite(suites)


def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)


if __name__ == "__main__":
    run(True)
