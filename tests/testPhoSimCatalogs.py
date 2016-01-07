from __future__ import with_statement
import os
import unittest
import lsst.utils.tests as utilsTests
import sqlite3
import numpy
import json
import lsst.utils
from collections import OrderedDict
from lsst.sims.catalogs.measures.instance import compound, CompoundInstanceCatalog
from lsst.sims.catUtils.utils import testStarsDBObj, testGalaxyDiskDBObj, \
                                     testGalaxyBulgeDBObj, testGalaxyAgnDBObj
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogSersic2D, PhoSimCatalogPoint, \
                                                         PhoSimCatalogZPoint
from lsst.sims.catalogs.generation.utils import makePhoSimTestDB

class PhoSimCatalogTest(unittest.TestCase):

    def setUp(self):
        self.obs_metadata = makePhoSimTestDB(size=10)
        self.bulgeDB = testGalaxyBulgeDBObj(driver='sqlite', database='PhoSimTestDatabase.db')
        self.diskDB = testGalaxyDiskDBObj(driver='sqlite', database='PhoSimTestDatabase.db')
        self.agnDB = testGalaxyAgnDBObj(driver='sqlite', database='PhoSimTestDatabase.db')
        self.starDB = testStarsDBObj(driver='sqlite', database='PhoSimTestDatabase.db')
        filter_translation={'u':0,'g':1, 'r':2, 'i':3, 'z':4, 'y':5}
        self.control_header = ['Opsim_moondec %.9g\n' % numpy.degrees(self.obs_metadata.phoSimMetaData['Opsim_moondec'][0]),
                              'Opsim_rottelpos %.9g\n' % numpy.degrees(self.obs_metadata.phoSimMetaData['Opsim_rottelpos'][0]),
                              'Unrefracted_Dec %.9g\n' % numpy.degrees(self.obs_metadata.phoSimMetaData['pointingDec'][0]),
                              'Opsim_moonalt %.9g\n' % numpy.degrees(self.obs_metadata.phoSimMetaData['Opsim_moonalt'][0]),
                              'Opsim_rotskypos %.9g\n' % numpy.degrees(self.obs_metadata.phoSimMetaData['Opsim_rotskypos'][0]),
                              'Opsim_moonra %.9g\n' % numpy.degrees(self.obs_metadata.phoSimMetaData['Opsim_moonra'][0]),
                              'Opsim_sunalt %.9g\n' % numpy.degrees(self.obs_metadata.phoSimMetaData['Opsim_sunalt'][0]),
                              'Opsim_expmjd %.9g\n' % self.obs_metadata.phoSimMetaData['Opsim_expmjd'][0],
                              'Unrefracted_Azimuth %.9g\n' % numpy.degrees(self.obs_metadata.phoSimMetaData['Unrefracted_Azimuth'][0]),
                              'Unrefracted_RA %.9g\n' % numpy.degrees(self.obs_metadata.phoSimMetaData['pointingRA'][0]),
                              'Opsim_dist2moon %.9g\n' % numpy.degrees(self.obs_metadata.phoSimMetaData['Opsim_dist2moon'][0]),
                              'Opsim_filter %d\n' % filter_translation[self.obs_metadata.phoSimMetaData['Opsim_filter'][0]],
                              'Unrefracted_Altitude %.9g\n' % numpy.degrees(self.obs_metadata.phoSimMetaData['Unrefracted_Altitude'][0])]


    def tearDown(self):

        if os.path.exists('PhoSimTestDatabase.db'):
            os.unlink('PhoSimTestDatabase.db')


    def verify_catalog(self, file_name):
        """
        Verify that the catalog specified by file_name has the expected header and is not empty
        """
        with open(file_name, 'r') as input_file:
            lines = input_file.readlines()
            for control_line in self.control_header:
                self.assertIn(control_line, lines)

            self.assertGreater(len(lines), len(self.control_header)+3)


    def testCatalog(self):
        """
        This test writes a PhoSim input catalog and compares it, one line at a time
        to a previously written catalog that should be identical.
        """
        testBulge = PhoSimCatalogSersic2D(self.bulgeDB, obs_metadata = self.obs_metadata)
        testDisk = PhoSimCatalogSersic2D(self.diskDB, obs_metadata = self.obs_metadata)
        testAgn = PhoSimCatalogZPoint(self.agnDB, obs_metadata = self.obs_metadata)
        testStar = PhoSimCatalogPoint(self.starDB, obs_metadata = self.obs_metadata)

        catName = os.path.join(lsst.utils.getPackageDir('sims_catUtils'),
                               'tests','scratchSpace','phoSimTestCatalog.txt')

        testBulge.write_catalog(catName)
        testDisk.write_catalog(catName, write_header=False, write_mode='a')
        testAgn.write_catalog(catName, write_header=False, write_mode='a')
        testStar.write_catalog(catName, write_header=False, write_mode='a')

        self.verify_catalog(catName)

        if os.path.exists(catName):
            os.unlink(catName)


    def testCompoundCatalog(self):
        """
        This test writes a PhoSim input catalog and compares it, one line at a time
        to a previously written catalog that should be identical.

        This test uses CompoundInstanceCatalog
        """

        # first, generate the catalog without a CompoundInstanceCatalog
        single_catName = os.path.join(lsst.utils.getPackageDir('sims_catUtils'),
                                      'tests','scratchSpace','phoSimTestCatalog_single.txt')

        testBulge = PhoSimCatalogSersic2D(self.bulgeDB, obs_metadata = self.obs_metadata)
        testDisk = PhoSimCatalogSersic2D(self.diskDB, obs_metadata = self.obs_metadata)
        testAgn = PhoSimCatalogZPoint(self.agnDB, obs_metadata = self.obs_metadata)
        testStar = PhoSimCatalogPoint(self.starDB, obs_metadata = self.obs_metadata)

        testBulge.write_catalog(single_catName)
        testDisk.write_catalog(single_catName, write_header=False, write_mode='a')
        testAgn.write_catalog(single_catName, write_header=False, write_mode='a')
        testStar.write_catalog(single_catName, write_header=False, write_mode='a')

        # now, generate the catalog using CompoundInstanceCatalog
        #
        # because the CompoundCatalogDBObject requires that database
        # connection parameters be set in the input CatalogDBObject
        # daughter class definitions, we have to declare dummy
        # CatalogDBObject daughter classes below

        class dummyDBbase(object):
            driver = 'sqlite'
            database = 'PhoSimTestDatabase.db'

        class dummyBulgeDB(dummyDBbase, testGalaxyBulgeDBObj):
            objid = 'dummy_bulge'

        class dummyDiskDB(dummyDBbase,  testGalaxyDiskDBObj):
            objid = 'dummy_disk'

        class dummyAgnDB(dummyDBbase, testGalaxyAgnDBObj):
            objid = 'dummy_agn'

        class dummyStarDB(dummyDBbase, testStarsDBObj):
            objid = 'dummy_stars'


        compoundCatalog = CompoundInstanceCatalog([PhoSimCatalogSersic2D, PhoSimCatalogSersic2D,
                                                   PhoSimCatalogZPoint, PhoSimCatalogPoint],
                                                  [dummyBulgeDB, dummyDiskDB, dummyAgnDB, dummyStarDB],
                                                  obs_metadata=self.obs_metadata)

        self.assertEqual(len(compoundCatalog._dbObjectGroupList[0]), 3)

        compound_catName = os.path.join(lsst.utils.getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                                        'phoSimTestCatalog_compound.txt')

        compoundCatalog.write_catalog(compound_catName)

        # verify that the compound catalog is what we expect
        self.verify_catalog(compound_catName)

        # verify that the two catalogs are equivalent
        with open(single_catName, 'r') as single_file:
            with open(compound_catName, 'r') as compound_file:
                single_lines = single_file.readlines()
                compound_lines = compound_file.readlines()

                for line in single_lines:
                    self.assertIn(line, compound_lines)

                for line in compound_lines:
                    self.assertIn(line, single_lines)


        if os.path.exists(compound_catName):
            os.unlink(compound_catName)

        if os.path.exists(single_catName):
            os.unlink(single_catName)


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(PhoSimCatalogTest)

    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)
if __name__ == "__main__":
    run(True)
