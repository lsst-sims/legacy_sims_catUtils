#One thing to be aware of: If the logic in Astrometry.py or Photometry.py changes,
#this unittest will fail, even if it is still possible to generate phoSim input
#catalogs.  This test is based on a simple line-by-line comparison of a phoSim
#input catalog with a previously generated catalog that should be identical.
#If we end up changing the logic in Astromtetry or Photometry, we will need to
#re-generate the testData/phoSimControlCatalog.txt to be consistent with the new code.

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
        baseLineFileName = lsst.utils.getPackageDir('sims_catUtils')+'/tests/testData/phoSimControlCatalog.txt'
        self.baseLineFile = open(baseLineFileName,'r')

    def tearDown(self):

        self.baseLineFile.close()

        if os.path.exists('PhoSimTestDatabase.db'):
            os.unlink('PhoSimTestDatabase.db')

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

        testFile = open(catName,'r')
        testLines = testFile.readlines()
        testFile.close()
        controlLines = self.baseLineFile.readlines()
        for line in testLines:
            msg = '%s not in controlLines' % line
            self.assertTrue(line in controlLines, msg=msg)

        for line in controlLines:
            msg = '%s not in testLines'
            self.assertTrue(line in testLines, msg=msg)

        if os.path.exists(catName):
            os.unlink(catName)


    def testCompoundCatalog(self):
        """
        This test writes a PhoSim input catalog and compares it, one line at a time
        to a previously written catalog that should be identical.

        This test uses CompoundInstanceCatalog
        """
        testBulge = PhoSimCatalogSersic2D(self.bulgeDB, obs_metadata = self.obs_metadata)
        testDisk = PhoSimCatalogSersic2D(self.diskDB, obs_metadata = self.obs_metadata)
        testAgn = PhoSimCatalogZPoint(self.agnDB, obs_metadata = self.obs_metadata)
        testStar = PhoSimCatalogPoint(self.starDB, obs_metadata = self.obs_metadata)

        compoundCatalog = CompoundInstanceCatalog([testBulge, testDisk, testAgn, testStar],
                                                   obs_metadata=self.obs_metadata)

        self.assertEqual(len(compoundCatalog._dbObjectGroupList[0]), 3)

        catName = os.path.join(lsst.utils.getPackageDir('sims_catUtils'), 'tests', 'scratchSpace',
                               'phoSimTestCatalog_compound.txt')

        compoundCatalog.write_catalog(catName)

        testFile = open(catName,'r')
        testLines = testFile.readlines()
        testFile.close()
        controlLines = self.baseLineFile.readlines()
        for line in testLines:
            msg = '%s not in controlLines' % line
            self.assertTrue(line in controlLines, msg=msg)

        for line in controlLines:
            msg = '%s not in testLines'
            self.assertTrue(line in testLines, msg=msg)

        if os.path.exists(catName):
            os.unlink(catName)



def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(PhoSimCatalogTest)

    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)
if __name__ == "__main__":
    run(True)
