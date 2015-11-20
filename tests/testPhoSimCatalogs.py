#One thing to be aware of: If the logic in Astrometry.py or Photometry.py changes,
#this unittest will fail, even if it is still possible to generate phoSim input
#catalogs.  This test is based on a simple line-by-line comparison of a phoSim
#input catalog with a previously generated catalog that should be identical.
#If we end up changing the logic in Astromtetry or Photometry, we will need to
#re-generate the testData/phoSimControlCatalog.txt to be consistent with the new code.

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
        self.control_header = ['Opsim_moondec 29.666667\n',
                              'Opsim_rottelpos 180\n',
                              'Unrefracted_Dec -29.666667\n',
                              'Opsim_moonalt -90\n',
                              'Opsim_rotskypos 0\n',
                              'Opsim_moonra 298.825632\n',
                              'Opsim_sunalt -90\n',
                              'Opsim_expmjd 52000\n',
                              'Unrefracted_Azimuth 0\n',
                              'Unrefracted_RA 118.825632\n',
                              'Opsim_dist2moon 179.999998\n',
                              'Opsim_filter 2\n',
                              'Unrefracted_Altitude 1.57079633\n']


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
        testBulge = PhoSimCatalogSersic2D(self.bulgeDB, obs_metadata = self.obs_metadata)
        testDisk = PhoSimCatalogSersic2D(self.diskDB, obs_metadata = self.obs_metadata)
        testAgn = PhoSimCatalogZPoint(self.agnDB, obs_metadata = self.obs_metadata)
        testStar = PhoSimCatalogPoint(self.starDB, obs_metadata = self.obs_metadata)

        # first, generate the catalog without a CompoundInstanceCatalog
        single_catName = os.path.join(lsst.utils.getPackageDir('sims_catUtils'),
                                      'tests','scratchSpace','phoSimTestCatalog_single.txt')

        # now, generate the catalog using CompoundInstanceCatalog
        testBulge.write_catalog(single_catName)
        testDisk.write_catalog(single_catName, write_header=False, write_mode='a')
        testAgn.write_catalog(single_catName, write_header=False, write_mode='a')
        testStar.write_catalog(single_catName, write_header=False, write_mode='a')


        compoundCatalog = CompoundInstanceCatalog([testBulge, testDisk, testAgn, testStar],
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
