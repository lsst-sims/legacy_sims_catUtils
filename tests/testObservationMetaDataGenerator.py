from __future__ import with_statement
import os
import unittest
import numpy
import lsst.utils.tests as utilsTests
from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from lsst.sims.utils import CircleBounds, BoxBounds, altAzPaFromRaDec
from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogSersic2D
from lsst.sims.catUtils.utils import testGalaxyBulgeDBObj
from lsst.sims.catalogs.generation.utils import makePhoSimTestDB
from lsst.utils import getPackageDir


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
        return 1.0/(numpy.cos(numpy.pi/2.0-numpy.radians(alt)))
    elif tag == 'm5':
        return obs.m5[obs.bandpass]
    elif tag == 'skyBrightness':
        return obs.skyBrightness
    elif tag == 'seeing':
        return obs.seeing[obs.bandpass]
    elif tag == 'telescopeFilter':
        return obs.bandpass

    mapping = {'obsHistID': 'Opsim_obshistid',
               'rawSeeing': 'Opsim_rawseeing',
               'dist2Moon': 'Opsim_dist2moon',
               'visitExpTime': 'exptime',
               'moonAlt': 'Opsim_moonalt',
               'sunAlt': 'Opsim_sunalt',
               'expDate': 'SIM_SEED',
               'moonRA': 'Opsim_moonra',
               'moonDec': 'Opsim_moondec',
               'moonPhase': 'Opsim_moonphase'}

    return obs.phoSimMetaData[mapping[tag]]


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

        return numpy.degrees(rec[tag])

    return rec[tag]


class ObservationMetaDataGeneratorTest(unittest.TestCase):

    longMessage = True

    def setUp(self):
        # these are the names of header fields that are handled from the standard ObservationMetaData
        # data, and thus should not be expected to be in the phoSimMetaData
        self.special_field_names = ('pointingRA', 'pointingDec', 'Opsim_altitude', 'Opsim_azimuth',
                                    'Opsim_expmjd', 'airmass', 'Opsim_filter', 'Opsim_rotskypos',
                                    'Unrefracted_RA', 'Unrefracted_Dec')

        dbPath = os.path.join(getPackageDir('sims_data'),
                              'OpSimData/opsimblitz1_1133_sqlite.db')
        self.gen = ObservationMetaDataGenerator(database=dbPath,
                                                driver='sqlite')

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
                  ('fieldRA', (numpy.degrees(1.370916), numpy.degrees(1.40))),
                  ('rawSeeing', (0.728562, 0.9)),
                  ('seeing', (0.7, 0.9)),
                  ('dist2Moon', (numpy.degrees(1.570307), numpy.degrees(1.9))),
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
                self.assertGreater(val, line[1][0], msg=msg)
                self.assertLess(val, line[1][1], msg=msg)

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
                    self.assertGreater(v1, bounds[ix][1][0], msg=msg)
                    self.assertLess(v1, bounds[ix][1][1], msg=msg)
                    self.assertGreater(v2, bounds[jx][1][0], msg=msg)
                    self.assertLess(v2, bounds[jx][1][1], msg=msg)

    def testOpSimQueryOnRanges(self):
        """
        Test that getOpimRecords() returns correct results
        """
        bounds = [('obsHistID', (5973, 7000)),
                  ('fieldRA', (numpy.degrees(1.370916), numpy.degrees(1.40))),
                  ('rawSeeing', (0.728562, 0.9)),
                  ('seeing', (0.7, 0.9)),
                  ('dist2Moon', (numpy.degrees(1.570307), numpy.degrees(1.9))),
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
                self.assertGreater(val, line[1][0], msg=msg)
                self.assertLess(val, line[1][1], msg=msg)

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
                    self.assertGreater(v1, bounds[ix][1][0], msg=msg)
                    self.assertLess(v1, bounds[ix][1][1], msg=msg)
                    self.assertGreater(v2, bounds[jx][1][0], msg=msg)
                    self.assertLess(v2, bounds[jx][1][1], msg=msg)

    def testQueryExactValues(self):
        """
        Test that ObservationMetaData returned by a query demanding an exact value do,
        in fact, adhere to that requirement.
        """
        gen = self.gen

        bounds = [('obsHistID', 5973),
                  ('expDate', 1220779),
                  ('fieldRA', numpy.degrees(1.370916)),
                  ('fieldDec', numpy.degrees(-0.456238)),
                  ('moonRA', numpy.degrees(2.914132)),
                  ('moonDec', numpy.degrees(0.06305)),
                  ('rotSkyPos', numpy.degrees(3.116656)),
                  ('telescopeFilter', 'i'),
                  ('rawSeeing', 0.728562),
                  ('seeing', 0.88911899999999999),
                  ('sunAlt', numpy.degrees(-0.522905)),
                  ('moonAlt', numpy.degrees(0.099096)),
                  ('dist2Moon', numpy.degrees(1.570307)),
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
                  ('fieldRA', numpy.degrees(1.370916)),
                  ('fieldDec', numpy.degrees(-0.456238)),
                  ('moonRA', numpy.degrees(2.914132)),
                  ('moonDec', numpy.degrees(0.06305)),
                  ('rotSkyPos', numpy.degrees(3.116656)),
                  ('telescopeFilter', 'i'),
                  ('rawSeeing', 0.728562),
                  ('seeing', 0.88911899999999999),
                  ('sunAlt', numpy.degrees(-0.522905)),
                  ('moonAlt', numpy.degrees(0.099096)),
                  ('dist2Moon', numpy.degrees(1.570307)),
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

        pointing_list = self.gen.getOpSimRecords(fieldRA=numpy.degrees(1.370916))
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
        results = gen.getObservationMetaData(fieldRA=(numpy.degrees(1.370916), numpy.degrees(1.5348635)),
                                             limit=20)
        self.assertEqual(len(results), 20)

    def testQueryOnFilter(self):
        """
        Test that queries on the filter work.
        """
        gen = self.gen
        results = gen.getObservationMetaData(fieldRA=numpy.degrees(1.370916), telescopeFilter='i')
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
        results = gen.getObservationMetaData(fieldRA=numpy.degrees(1.370916),
                                             telescopeFilter='i',
                                             boundLength=0.9)
        ct = 0
        for obs_metadata in results:
            self.assertTrue(isinstance(obs_metadata.bounds, CircleBounds))

            # include some wiggle room, in case ObservationMetaData needs to
            # adjust the boundLength to accommodate the transformation between
            # ICRS and observed coordinates
            self.assertGreaterEqual(obs_metadata.bounds.radiusdeg, 0.9)
            self.assertLess(obs_metadata.bounds.radiusdeg, 0.95)

            self.assertAlmostEqual(obs_metadata.bounds.RA,
                                   numpy.radians(obs_metadata.pointingRA), 5)
            self.assertAlmostEqual(obs_metadata.bounds.DEC,
                                   numpy.radians(obs_metadata.pointingDec), 5)
            ct += 1

        # Make sure that some ObservationMetaData were tested
        self.assertGreater(ct, 0)

        boundLengthList = [1.2, (1.2, 0.6)]
        for boundLength in boundLengthList:
            results = gen.getObservationMetaData(fieldRA=numpy.degrees(1.370916),
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
                self.assertTrue(isinstance(obs_metadata.bounds, BoxBounds))

                self.assertAlmostEqual(obs_metadata.bounds.RAminDeg, RAdeg-dra, 10)

                self.assertAlmostEqual(obs_metadata.bounds.RAmaxDeg, RAdeg+dra, 10)

                self.assertAlmostEqual(obs_metadata.bounds.DECminDeg, DECdeg-ddec, 10)

                self.assertAlmostEqual(obs_metadata.bounds.DECmaxDeg, DECdeg+ddec, 10)

                self.assertAlmostEqual(obs_metadata.bounds.RA, numpy.radians(obs_metadata.pointingRA), 5)
                self.assertAlmostEqual(obs_metadata.bounds.DEC, numpy.radians(obs_metadata.pointingDec), 5)

                ct += 1

            # Make sure that some ObservationMetaData were tested
            self.assertGreater(ct, 0)

    def testCreationOfPhoSimCatalog(self):
        """
        Make sure that we can create PhoSim input catalogs using the returned
        ObservationMetaData. This test will just make sure that all of the
        expected header entries are there.
        """

        dbName = 'obsMetaDataGeneratorTest.db'
        catName = 'testPhoSimFromObsMetaDataGenerator.txt'
        if os.path.exists(dbName):
            os.unlink(dbName)
        makePhoSimTestDB(filename=dbName)
        bulgeDB = testGalaxyBulgeDBObj(driver='sqlite', database=dbName)
        gen = self.gen
        results = gen.getObservationMetaData(fieldRA=numpy.degrees(1.370916),
                                             telescopeFilter='i')
        testCat = PhoSimCatalogSersic2D(bulgeDB, obs_metadata=results[0])
        testCat.write_catalog(catName)

        with open(catName) as inputFile:
            lines = inputFile.readlines()
            header_entries = []
            for line in lines:
                words = line.split()
                if words[0] == 'object':
                    break
                header_entries.append(words[0])

            for column in self.gen._opsim_to_phosim:
                new_name = self.gen._opsim_to_phosim[column][0]
                self.assertIn(new_name, header_entries)

        if os.path.exists(catName):
            os.unlink(catName)

        if os.path.exists(dbName):
            os.unlink(dbName)


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(ObservationMetaDataGeneratorTest)

    return unittest.TestSuite(suites)


def run(shouldExit=False):
    utilsTests.run(suite(), shouldExit)
if __name__ == "__main__":
    run(True)
