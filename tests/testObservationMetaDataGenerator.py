from __future__ import with_statement
import os
import unittest
import sqlite3
import numpy as np
import lsst.utils.tests
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.utils import CircleBounds, BoxBounds, altAzPaFromRaDec
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogSersic2D
from lsst.sims.catUtils.utils import testGalaxyBulgeDBObj
from lsst.sims.catalogs.utils import makePhoSimTestDB
from lsst.utils import getPackageDir


def setup_module(module):
    lsst.utils.tests.init()


def get_val_from_obs(tag, obs):
    """
    tag is the name of a data column

    obs is an ObservationMetaData

    returns the value of 'tag' for this obs
    """
    if tag == 'fieldRA':
        return obs.pointingRA
    elif tag == 'fieldDec':
        return obs.pointingDec
    elif tag == 'rotSkyPos':
        return obs.rotSkyPos
    elif tag == 'expMJD':
        return obs.mjd.TAI
    elif tag == 'airmass':
        alt, az, pa = altAzPaFromRaDec(obs.pointingRA, obs.pointingDec, obs)
        return 1.0/(np.cos(np.pi/2.0-np.radians(alt)))
    elif tag == 'm5':
        return obs.m5[obs.bandpass]
    elif tag == 'skyBrightness':
        return obs.skyBrightness
    elif tag == 'seeing':
        return obs.seeing[obs.bandpass]
    elif tag == 'telescopeFilter':
        return obs.bandpass

    transforms = {'dist2Moon': np.degrees,
                  'moonAlt': np.degrees,
                  'sunAlt': np.degrees,
                  'moonRA': np.degrees,
                  'moonDec': np.degrees}

    if tag in transforms:
        return transforms[tag](obs.OpsimMetaData[tag])

    return obs.OpsimMetaData[tag]


def get_val_from_rec(tag, rec):
    """
    tag is the name of a data column

    rec is a record returned by a query to an OpSim database

    returns the value of 'tag' for this rec
    """
    if tag == 'telescopeFilter':
        return rec['filter']
    elif tag == 'seeing':
        return rec['finSeeing']  # because the test opsim database uses the old schema
    elif tag == 'm5':
        return rec['fiveSigmaDepth']
    elif tag == 'skyBrightness':
        return rec['filtSkyBrightness']
    elif tag in ('fieldRA', 'fieldDec', 'moonRA', 'moonDec',
                 'rotSkyPos', 'sunAlt', 'moonAlt', 'dist2Moon', 'altitude',
                 'azimuth'):

        return np.degrees(rec[tag])

    return rec[tag]


class ObservationMetaDataGeneratorTest(unittest.TestCase):

    longMessage = True

    def setUp(self):

        dbPath = os.path.join(getPackageDir('sims_data'),
                              'OpSimData/opsimblitz1_1133_sqlite.db')

        self.gen = ObservationMetaDataGenerator(database=dbPath,
                                                driver='sqlite')

    def tearDown(self):
        del self.gen

    def testExceptions(self):
        """
        Make sure that RuntimeErrors get raised when they should
        """
        gen = self.gen
        self.assertRaises(RuntimeError, gen.getObservationMetaData)
        self.assertRaises(RuntimeError, gen.getObservationMetaData, fieldRA=(1.0, 2.0, 3.0))

    def testOrderBy(self):
        """
        Test that the ObservationMetaData really do get ordered on expMJD
        """
        results = self.gen.getObservationMetaData(altitude=(50.0, 80.0), limit=100)
        self.assertEqual(len(results), 100)
        for ix, obs in enumerate(results):
            self.assertGreaterEqual(obs.OpsimMetaData['altitude'], np.radians(50.0))
            self.assertLessEqual(obs.OpsimMetaData['altitude'], np.radians(80.0))
            if ix>1:
                self.assertGreater(obs.mjd.TAI, results[ix-1].mjd.TAI)

    def testQueryOnRanges(self):
        """
        Test that ObservationMetaData objects returned by queries of the form
        min < value < max
        are, in fact, within that range.

        Test when querying on both a single and two columns.
        """
        gen = self.gen

        # An list containing the bounds of our queries.
        # The order of the tuples must correspond to the order of
        # self.columnMapping in ObservationMetaDataGenerator.
        # This was generated with a separate script which printed
        # the median and maximum values of all of the quantities
        # in our test opsim database
        bounds = [('obsHistID', (5973, 7000)),
                  ('fieldRA', (np.degrees(1.370916), np.degrees(1.40))),
                  ('rawSeeing', (0.728562, 0.9)),
                  ('seeing', (0.7, 0.9)),
                  ('dist2Moon', (np.degrees(1.570307), np.degrees(1.9))),
                  ('expMJD', (49367.129396, 49370.0)),
                  ('m5', (22.815249, 23.0)),
                  ('skyBrightness', (19.017605, 19.5))]

        # test querying on a single column
        for line in bounds:
            tag = line[0]

            args = {tag: line[1]}

            results = gen.getObservationMetaData(**args)
            msg = "failed querying on %s" % tag
            self.assertGreater(len(results), 0, msg=msg)

            for obs in results:
                val = get_val_from_obs(tag, obs)
                self.assertGreaterEqual(val, line[1][0], msg=msg)
                self.assertLessEqual(val, line[1][1], msg=msg)

        # test querying on two columns at once
        for ix in range(len(bounds)):
            tag1 = bounds[ix][0]

            for jx in range(ix+1, len(bounds)):
                tag2 = bounds[jx][0]

                args = {}
                args[tag1] = bounds[ix][1]
                args[tag2] = bounds[jx][1]
                results = gen.getObservationMetaData(**args)
                msg = "failed querying %s and %s" % (tag1, tag2)
                self.assertGreater(len(results), 0, msg=msg)
                for obs in results:
                    v1 = get_val_from_obs(tag1, obs)
                    v2 = get_val_from_obs(tag2, obs)
                    self.assertGreaterEqual(v1, bounds[ix][1][0], msg=msg)
                    self.assertLessEqual(v1, bounds[ix][1][1], msg=msg)
                    self.assertGreaterEqual(v2, bounds[jx][1][0], msg=msg)
                    self.assertLessEqual(v2, bounds[jx][1][1], msg=msg)

    def testOpSimQueryOnRanges(self):
        """
        Test that getOpimRecords() returns correct results
        """
        bounds = [('obsHistID', (5973, 7000)),
                  ('fieldRA', (np.degrees(1.370916), np.degrees(1.40))),
                  ('rawSeeing', (0.728562, 0.9)),
                  ('seeing', (0.7, 0.9)),
                  ('dist2Moon', (np.degrees(1.570307), np.degrees(1.9))),
                  ('expMJD', (49367.129396, 49370.0)),
                  ('m5', (22.815249, 23.0)),
                  ('skyBrightness', (19.017605, 19.5))]

        for line in bounds:
            tag = line[0]
            args = {tag: line[1]}
            results = self.gen.getOpSimRecords(**args)
            msg = 'failed querying %s ' % tag
            self.assertGreater(len(results), 0)
            for rec in results:
                val = get_val_from_rec(tag, rec)
                self.assertGreaterEqual(val, line[1][0], msg=msg)
                self.assertLessEqual(val, line[1][1], msg=msg)

        for ix in range(len(bounds)):
            tag1 = bounds[ix][0]
            for jx in range(ix+1, len(bounds)):
                tag2 = bounds[jx][0]
                args = {tag1: bounds[ix][1], tag2: bounds[jx][1]}
                results = self.gen.getOpSimRecords(**args)
                msg = 'failed while querying %s and %s' % (tag1, tag2)
                self.assertGreater(len(results), 0)
                for rec in results:
                    v1 = get_val_from_rec(tag1, rec)
                    v2 = get_val_from_rec(tag2, rec)
                    self.assertGreaterEqual(v1, bounds[ix][1][0], msg=msg)
                    self.assertLessEqual(v1, bounds[ix][1][1], msg=msg)
                    self.assertGreaterEqual(v2, bounds[jx][1][0], msg=msg)
                    self.assertLessEqual(v2, bounds[jx][1][1], msg=msg)

    def testQueryExactValues(self):
        """
        Test that ObservationMetaData returned by a query demanding an exact value do,
        in fact, adhere to that requirement.
        """
        gen = self.gen

        bounds = [('obsHistID', 5973),
                  ('expDate', 1220779),
                  ('fieldRA', np.degrees(1.370916)),
                  ('fieldDec', np.degrees(-0.456238)),
                  ('moonRA', np.degrees(2.914132)),
                  ('moonDec', np.degrees(0.06305)),
                  ('rotSkyPos', np.degrees(3.116656)),
                  ('telescopeFilter', 'i'),
                  ('rawSeeing', 0.728562),
                  ('seeing', 0.88911899999999999),
                  ('sunAlt', np.degrees(-0.522905)),
                  ('moonAlt', np.degrees(0.099096)),
                  ('dist2Moon', np.degrees(1.570307)),
                  ('moonPhase', 52.2325),
                  ('expMJD', 49367.129396),
                  ('visitExpTime', 30.0),
                  ('m5', 22.815249),
                  ('skyBrightness', 19.017605)]

        for ii in range(len(bounds)):
            tag = bounds[ii][0]
            args = {}
            args[tag] = bounds[ii][1]
            results = gen.getObservationMetaData(**args)
            msg = 'failed querying %s' % tag
            self.assertGreater(len(results), 0, msg=msg)
            for obs in results:
                self.assertEqual(get_val_from_obs(tag, obs), bounds[ii][1], msg=msg)

    def testOpSimQueryExact(self):
        """
        Test that querying OpSim records for exact values works
        """

        bounds = [('obsHistID', 5973),
                  ('expDate', 1220779),
                  ('fieldRA', np.degrees(1.370916)),
                  ('fieldDec', np.degrees(-0.456238)),
                  ('moonRA', np.degrees(2.914132)),
                  ('moonDec', np.degrees(0.06305)),
                  ('rotSkyPos', np.degrees(3.116656)),
                  ('telescopeFilter', 'i'),
                  ('rawSeeing', 0.728562),
                  ('seeing', 0.88911899999999999),
                  ('sunAlt', np.degrees(-0.522905)),
                  ('moonAlt', np.degrees(0.099096)),
                  ('dist2Moon', np.degrees(1.570307)),
                  ('moonPhase', 52.2325),
                  ('expMJD', 49367.129396),
                  ('visitExpTime', 30.0),
                  ('m5', 22.815249),
                  ('skyBrightness', 19.017605)]

        for line in bounds:
            tag = line[0]
            args = {tag: line[1]}
            results = self.gen.getOpSimRecords(**args)
            msg = 'failed while querying %s' % tag
            self.assertGreater(len(results), 0, msg=msg)
            for rec in results:
                self.assertEqual(get_val_from_rec(tag, rec), line[1], msg=msg)

    def testPassInOtherQuery(self):
        """
        Test that you can pass OpSim pointings generated from another source
        into an ObservationMetaDataGenerator and still get ObservationMetaData
        out
        """

        pointing_list = self.gen.getOpSimRecords(fieldRA=np.degrees(1.370916))
        self.assertGreater(len(pointing_list), 1)
        local_gen = ObservationMetaDataGenerator()
        obs_list = local_gen.ObservationMetaDataFromPointingArray(pointing_list)
        self.assertEqual(len(obs_list), len(pointing_list))

        for pp in pointing_list:
            obs = local_gen.ObservationMetaDataFromPointing(pp)
            self.assertIsInstance(obs, ObservationMetaData)

    def testQueryLimit(self):
        """
        Test that, when we specify a limit on the number of ObservationMetaData we want returned,
        that limit is respected
        """
        gen = self.gen
        results = gen.getObservationMetaData(fieldRA=(np.degrees(1.370916), np.degrees(1.5348635)),
                                             limit=20)
        self.assertEqual(len(results), 20)

    def testQueryOnFilter(self):
        """
        Test that queries on the filter work.
        """
        gen = self.gen
        results = gen.getObservationMetaData(fieldRA=np.degrees(1.370916), telescopeFilter='i')
        ct = 0
        for obs_metadata in results:
            self.assertAlmostEqual(obs_metadata._pointingRA, 1.370916)
            self.assertEqual(obs_metadata.bandpass, 'i')
            ct += 1

        # Make sure that more than zero ObservationMetaData were returned
        self.assertGreater(ct, 0)

    def testObsMetaDataBounds(self):
        """
        Make sure that the bound specifications (i.e. a circle or a box on the
        sky) are correctly passed through to the resulting ObservationMetaData
        """

        gen = self.gen

        # Test a cirlce with a specified radius
        results = gen.getObservationMetaData(fieldRA=np.degrees(1.370916),
                                             telescopeFilter='i',
                                             boundLength=0.9)
        ct = 0
        for obs_metadata in results:
            self.assertTrue(isinstance(obs_metadata.bounds, CircleBounds),
                            msg='obs_metadata.bounds is not an intance of '
                            'CircleBounds')

            # include some wiggle room, in case ObservationMetaData needs to
            # adjust the boundLength to accommodate the transformation between
            # ICRS and observed coordinates
            self.assertGreaterEqual(obs_metadata.bounds.radiusdeg, 0.9)
            self.assertLessEqual(obs_metadata.bounds.radiusdeg, 0.95)

            self.assertAlmostEqual(obs_metadata.bounds.RA,
                                   np.radians(obs_metadata.pointingRA), 5)
            self.assertAlmostEqual(obs_metadata.bounds.DEC,
                                   np.radians(obs_metadata.pointingDec), 5)
            ct += 1

        # Make sure that some ObservationMetaData were tested
        self.assertGreater(ct, 0)

        boundLengthList = [1.2, (1.2, 0.6)]
        for boundLength in boundLengthList:
            results = gen.getObservationMetaData(fieldRA=np.degrees(1.370916),
                                                 telescopeFilter='i',
                                                 boundType='box',
                                                 boundLength=boundLength)

            if hasattr(boundLength, '__len__'):
                dra = boundLength[0]
                ddec = boundLength[1]
            else:
                dra = boundLength
                ddec = boundLength

            ct = 0
            for obs_metadata in results:
                RAdeg = obs_metadata.pointingRA
                DECdeg = obs_metadata.pointingDec
                self.assertTrue(isinstance(obs_metadata.bounds, BoxBounds),
                                msg='obs_metadata.bounds is not an instance of '
                                'BoxBounds')

                self.assertAlmostEqual(obs_metadata.bounds.RAminDeg, RAdeg-dra, 10)

                self.assertAlmostEqual(obs_metadata.bounds.RAmaxDeg, RAdeg+dra, 10)

                self.assertAlmostEqual(obs_metadata.bounds.DECminDeg, DECdeg-ddec, 10)

                self.assertAlmostEqual(obs_metadata.bounds.DECmaxDeg, DECdeg+ddec, 10)

                self.assertAlmostEqual(obs_metadata.bounds.RA, np.radians(obs_metadata.pointingRA), 5)
                self.assertAlmostEqual(obs_metadata.bounds.DEC, np.radians(obs_metadata.pointingDec), 5)

                ct += 1

            # Make sure that some ObservationMetaData were tested
            self.assertGreater(ct, 0)

    def testQueryOnNight(self):
        """
        Check that the ObservationMetaDataGenerator can query on the 'night'
        column in the OpSim summary table
        """

        # the test database goes from night=0 to night=28
        # corresponding to 49353.032079 <= mjd <= 49381.38533
        night0 = 49353.032079

        results = self.gen.getObservationMetaData(night=(11, 13))

        self.assertGreater(len(results), 1800)
        # there should be about 700 observations a night;
        # make sure we get at least 600

        for obs in results:
            self.assertGreaterEqual(obs.mjd.TAI, night0+11.0)
            self.assertLessEqual(obs.mjd.TAI, night0+13.5)
            # the 0.5 is there because the last observation on night 13 could be
            # 13 days and 8 hours after the first observation on night 0
            self.assertGreaterEqual(obs.OpsimMetaData['night'], 11)
            self.assertLessEqual(obs.OpsimMetaData['night'], 13)

        # query for an exact night
        results = self.gen.getObservationMetaData(night=15)

        self.assertGreater(len(results), 600)
        # there should be about 700 observations a night;
        # make sure we get at least 600

        for obs in results:
            self.assertEqual(obs.OpsimMetaData['night'], 15)
            self.assertGreaterEqual(obs.mjd.TAI, night0+14.9)
            self.assertLessEqual(obs.mjd.TAI, night0+15.9)

    def testQueryOnConstraint(self):
        """
        Test that we can create a list of ObservationMetaData with an
        arbitary SQL 'where' clause
        """

        constraint = 'where night<3 and altitude>1.3'
        results = self.gen.getObservationMetaDataFromConstraint(constraint)
        self.assertGreater(len(results), 0)
        for obs in results:
            self.assertLess(obs.OpsimMetaData['night'], 3)
            self.assertGreater(obs.OpsimMetaData['altitude'], 1.3)

        constraint = 'where altitude<1.2 limit 10'
        results = self.gen.getObservationMetaDataFromConstraint(constraint)
        self.assertEqual(len(results), 10)
        for obs in results:
            self.assertLess(obs.OpsimMetaData['altitude'], 1.2)

    def testCreationOfPhoSimCatalog(self):
        """
        Make sure that we can create PhoSim input catalogs using the returned
        ObservationMetaData. This test will just make sure that all of the
        expected header entries are there.
        """

        scratch_dir = os.path.join(getPackageDir('sims_catUtils'), 'tests',
                                   'scratchSpace')
        dbName = os.path.join(scratch_dir, 'obsMetaDataGeneratorTest.db')
        catName = os.path.join(scratch_dir, 'testPhoSimFromObsMetaDataGenerator.txt')
        if os.path.exists(dbName):
            os.unlink(dbName)
        makePhoSimTestDB(filename=dbName)
        bulgeDB = testGalaxyBulgeDBObj(driver='sqlite', database=dbName)
        gen = self.gen
        results = gen.getObservationMetaData(fieldRA=np.degrees(1.370916),
                                             telescopeFilter='i')
        testCat = PhoSimCatalogSersic2D(bulgeDB, obs_metadata=results[0])
        testCat.phoSimHeaderMap = {}
        testCat.write_catalog(catName)

        if os.path.exists(catName):
            os.unlink(catName)

        if os.path.exists(dbName):
            os.unlink(dbName)


class ObsMetaDataGenMockOpsimTest(unittest.TestCase):
    """
    This class will test the performance of the ObservationMetaDataGenerator
    on a 'mock OpSim' database (i.e. a database of pointings whose Summary
    table contains only a subset of the official OpSim schema)
    """

    @classmethod
    def setUpClass(cls):
        scratch_dir = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace')
        cls.opsim_db_name = os.path.join(scratch_dir, 'mock_opsim_sqlite.db')

        if os.path.exists(cls.opsim_db_name):
            os.unlink(cls.opsim_db_name)

        conn = sqlite3.connect(cls.opsim_db_name)
        c = conn.cursor()
        c.execute('''CREATE TABLE Summary (obsHistID int, expMJD real, '''
                  '''fieldRA real, fieldDec real, filter text)''')
        conn.commit()
        rng = np.random.RandomState(77)
        n_pointings = 100
        ra_data = rng.random_sample(n_pointings)*2.0*np.pi
        dec_data = (rng.random_sample(n_pointings)-0.5)*np.pi
        mjd_data = rng.random_sample(n_pointings)*1000.0 + 59580.0
        filter_dexes = rng.randint(0, 6, n_pointings)
        bands = ('u', 'g', 'r', 'i', 'z', 'y')
        filter_data = []
        for ii in filter_dexes:
            filter_data.append(bands[ii])

        for ii in range(n_pointings):
            cmd = '''INSERT INTO Summary VALUES(%i, %f, %f, %f, '%s')''' % \
                  (ii, mjd_data[ii], ra_data[ii], dec_data[ii], filter_data[ii])
            c.execute(cmd)
        conn.commit()
        conn.close()

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.opsim_db_name):
            os.unlink(cls.opsim_db_name)

    def setUp(self):
        self.obs_meta_gen = ObservationMetaDataGenerator(database=self.opsim_db_name)

    def tearDown(self):
        del self.obs_meta_gen

    def testSpatialQuery(self):
        """
        Test that when we run a spatial query on the mock opsim database, we get expected results.
        """

        raBounds = (45.0, 100.0)
        results = self.obs_meta_gen.getObservationMetaData(fieldRA=raBounds)
        self.assertGreater(len(results), 0)
        for obs in results:
            self.assertGreater(obs.pointingRA, raBounds[0])
            self.assertLessEqual(obs.pointingDec, raBounds[1])

    def testSelectException(self):
        """
        Test that an exception is raised if you try to SELECT pointings on a column that does not exist
        """
        with self.assertRaises(RuntimeError) as context:
            self.obs_meta_gen.getObservationMetaData(rotSkyPos=(27.0, 112.0))
        self.assertIn("You have asked ObservationMetaDataGenerator to SELECT",
                      context.exception.args[0])

    def testIncompletDB(self):
        """
        Test that if the mock OpSim database does not have all required columns, an exception
        is raised.
        """
        scratch_dir = os.path.join(getPackageDir('sims_catUtils'), 'tests', 'scratchSpace')
        opsim_db_name = os.path.join(scratch_dir, 'incomplete_mock_opsim_sqlite.db')

        if os.path.exists(opsim_db_name):
            os.unlink(opsim_db_name)

        conn = sqlite3.connect(opsim_db_name)
        c = conn.cursor()
        c.execute('''CREATE TABLE Summary (obsHistID int, expMJD real, '''
                  '''fieldRA real, filter text)''')
        conn.commit()

        rng = np.random.RandomState(77)
        n_pointings = 100
        ra_data = rng.random_sample(n_pointings)*2.0*np.pi
        mjd_data = rng.random_sample(n_pointings)*1000.0 + 59580.0
        filter_dexes = rng.randint(0, 6, n_pointings)
        bands = ('u', 'g', 'r', 'i', 'z', 'y')
        filter_data = []
        for ii in filter_dexes:
            filter_data.append(bands[ii])

        for ii in range(n_pointings):
            cmd = '''INSERT INTO Summary VALUES(%i, %f, %f, '%s')''' % \
                  (ii, mjd_data[ii], ra_data[ii], filter_data[ii])
            c.execute(cmd)
        conn.commit()
        conn.close()

        incomplete_obs_gen = ObservationMetaDataGenerator(database=opsim_db_name)

        with self.assertRaises(RuntimeError) as context:
            incomplete_obs_gen.getObservationMetaData(telescopeFilter='r')
        self.assertIn("ObservationMetaDataGenerator requires that the database",
                      context.exception.args[0])

        if os.path.exists(opsim_db_name):
            os.unlink(opsim_db_name)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
