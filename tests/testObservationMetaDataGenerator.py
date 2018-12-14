from __future__ import with_statement
from builtins import range
import os
import unittest
import sqlite3
import tempfile
import numpy as np
import lsst.utils.tests
import lsst.sims.utils.htmModule as htm
from lsst.sims.utils.CodeUtilities import sims_clean_up
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.utils import CircleBounds, BoxBounds, altAzPaFromRaDec
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogSersic2D
from lsst.sims.catUtils.exampleCatalogDefinitions import DefaultPhoSimHeaderMap
from lsst.sims.catUtils.utils import testGalaxyBulgeDBObj
from lsst.sims.catUtils.utils import makePhoSimTestDB
from lsst.utils import getPackageDir

ROOT = os.path.abspath(os.path.dirname(__file__))


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

    @classmethod
    def tearDownClass(cls):
        sims_clean_up()

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
            self.assertGreaterEqual(obs._OpsimMetaData['night'], 11)
            self.assertLessEqual(obs._OpsimMetaData['night'], 13)

        # query for an exact night
        results = self.gen.getObservationMetaData(night=15)

        self.assertGreater(len(results), 600)
        # there should be about 700 observations a night;
        # make sure we get at least 600

        for obs in results:
            self.assertEqual(obs._OpsimMetaData['night'], 15)
            self.assertGreaterEqual(obs.mjd.TAI, night0+14.9)
            self.assertLessEqual(obs.mjd.TAI, night0+15.9)

    def testCreationOfPhoSimCatalog(self):
        """
        Make sure that we can create PhoSim input catalogs using the returned
        ObservationMetaData. This test will just make sure that all of the
        expected header entries are there.
        """

        dbName = tempfile.mktemp(dir=ROOT, prefix='obsMetaDataGeneratorTest-', suffix='.db')
        makePhoSimTestDB(filename=dbName)
        bulgeDB = testGalaxyBulgeDBObj(driver='sqlite', database=dbName)
        gen = self.gen
        results = gen.getObservationMetaData(fieldRA=np.degrees(1.370916),
                                             telescopeFilter='i')
        testCat = PhoSimCatalogSersic2D(bulgeDB, obs_metadata=results[0])
        testCat.phoSimHeaderMap = {}
        with lsst.utils.tests.getTempFilePath('.txt') as catName:
            testCat.write_catalog(catName)

        if os.path.exists(dbName):
            os.unlink(dbName)

    def testCreationOfPhoSimCatalog_2(self):
        """
        Make sure that we can create PhoSim input catalogs using the returned
        ObservationMetaData.

        Use the actual DefaultPhoSimHeader map; make sure that opsim_version
        does not make it into the header.
        """

        dbName = tempfile.mktemp(dir=ROOT, prefix='obsMetaDataGeneratorTest-', suffix='.db')
        makePhoSimTestDB(filename=dbName)
        bulgeDB = testGalaxyBulgeDBObj(driver='sqlite', database=dbName)
        gen = self.gen
        results = gen.getObservationMetaData(fieldRA=np.degrees(1.370916),
                                             telescopeFilter='i')
        testCat = PhoSimCatalogSersic2D(bulgeDB, obs_metadata=results[0])
        testCat.phoSimHeaderMap = DefaultPhoSimHeaderMap
        with lsst.utils.tests.getTempFilePath('.txt') as catName:
            testCat.write_catalog(catName)
            ct_lines = 0
            with open(catName, 'r') as in_file:
                for line in in_file:
                    ct_lines += 1
                    self.assertNotIn('opsim_version', line)
                self.assertGreater(ct_lines, 10)  # check that some lines did get written

        if os.path.exists(dbName):
            os.unlink(dbName)

    def testCreationOfPhoSimCatalog_3(self):
        """
        Make sure that we can create PhoSim input catalogs using the returned
        ObservationMetaData.

        Test that an error is actually raised if we try to build a PhoSim catalog
        with a v3 header map using a v4 ObservationMetaData
        """

        dbName = tempfile.mktemp(dir=ROOT, prefix='obsMetaDataGeneratorTest-', suffix='.db')
        makePhoSimTestDB(filename=dbName)
        bulgeDB = testGalaxyBulgeDBObj(driver='sqlite', database=dbName)
        opsim_db = os.path.join(getPackageDir('sims_data'), 'OpSimData',
                                'astro-lsst-01_2014.db')
        assert os.path.isfile(opsim_db)
        gen = ObservationMetaDataGenerator(opsim_db, driver='sqlite')
        results = gen.getObservationMetaData(fieldRA=(70.0, 85.0),
                                             telescopeFilter='i')
        self.assertGreater(len(results), 0)
        testCat = PhoSimCatalogSersic2D(bulgeDB, obs_metadata=results[0])
        testCat.phoSimHeaderMap = DefaultPhoSimHeaderMap
        with lsst.utils.tests.getTempFilePath('.txt') as catName:
            with self.assertRaises(RuntimeError):
                testCat.write_catalog(catName)

        if os.path.exists(dbName):
            os.unlink(dbName)

    def test_htmid_query(self):
        """
        Test that queries on HTMID work
        """
        obs_gen = self.gen

        boundLength = 0.7
        htmid  = 33937

        tx = htm.trixelFromHtmid(htmid)
        radius = tx.get_radius()
        ra_c, dec_c = tx.get_center()

        # construct a set of rectangular bounds
        # that will encompass the entire trixel
        # plus extra pointings.  This will be our
        # "control" list
        d_ra_dec = 3.50+2.0*(radius+boundLength)
        ra_min = ra_c - d_ra_dec
        ra_max = ra_c + d_ra_dec
        dec_min = dec_c - d_ra_dec
        dec_max = dec_c + d_ra_dec

        obs_control_list = obs_gen.getObservationMetaData(fieldRA=(ra_min, ra_max),
                                                          fieldDec=(dec_min, dec_max),
                                                          boundLength=boundLength,
                                                          boundType='circle')

        control_id_list = [obs.OpsimMetaData['obsHistID']
                           for obs in obs_control_list]
        obs_control_dict = {}
        for obs in obs_control_list:
            obs_control_dict[obs.OpsimMetaData['obsHistID']] = obs

        obs_htmid_list = obs_gen.obsFromHtmid(htmid,
                                              boundLength=boundLength,
                                              transform_to_degrees=np.degrees)

        htmid_id_list = [obs.OpsimMetaData['obsHistID']
                         for obs in obs_htmid_list]


        obs_htmid_dict = {}
        for obs in obs_htmid_list:
            obs_htmid_dict[obs.OpsimMetaData['obsHistID']] = obs

        self.assertGreater(len(htmid_id_list), 5)
        self.assertGreater(len(control_id_list), len(htmid_id_list))

        # make sure that every pointing selected on htmid is
        # in the control set
        control_id_set = set(control_id_list)
        for hh in htmid_id_list:
            self.assertIn(hh, control_id_set)

        htmid_id_set = set(htmid_id_list)
        for cc in control_id_list:
            obs = obs_control_dict[cc]
            hs = htm.halfSpaceFromRaDec(obs.pointingRA, obs.pointingDec, obs.boundLength)
            if cc not in htmid_id_set:
                # if a control pointing is not in the set selected on htmid,
                # make sure it does not overlap the trixel
                self.assertEqual(hs.contains_trixel(tx), 'outside')
            else:
                # make sure that pointings in both sets overlap the trixel
                # and are identical
                self.assertNotEqual(hs.contains_trixel(tx), 'outside')
                self.assertEqual(obs_control_dict[cc], obs_htmid_dict[cc])


class ObsMetaDataGenMockOpsimTest(unittest.TestCase):
    """
    This class will test the performance of the ObservationMetaDataGenerator
    on a 'mock OpSim' database (i.e. a database of pointings whose Summary
    table contains only a subset of the official OpSim schema)
    """

    @classmethod
    def setUpClass(cls):
        cls.opsim_db_name = tempfile.mktemp(dir=ROOT, prefix='mock_opsim_sqlite-', suffix='.db')

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

        sims_clean_up()

        if os.path.exists(cls.opsim_db_name):
            os.unlink(cls.opsim_db_name)

    def setUp(self):
        self.obs_meta_gen = ObservationMetaDataGenerator(database=self.opsim_db_name)

    def tearDown(self):
        del self.obs_meta_gen

    def testOnNonExistentDatabase(self):
        """
        Test that an exception is raised if you try to connect to an query
        a database that does not exist.
        """

        test_name = 'non_existent.db'
        with self.assertRaises(RuntimeError) as context:
            ObservationMetaDataGenerator(database=test_name,
                                         driver='sqlite')

        self.assertEqual(context.exception.args[0],
                         '%s is not a file' % test_name)

        self.assertFalse(os.path.exists(test_name))

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
        opsim_db_name = tempfile.mktemp(dir=ROOT, prefix='incomplete_mock_opsim_sqlite-', suffix='.db')

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
