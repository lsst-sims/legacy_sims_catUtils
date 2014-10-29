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
from collections import OrderedDict
from lsst.sims.catalogs.measures.instance import compound
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catUtils.baseCatalogModels import GalaxyBulgeObj, GalaxyDiskObj, GalaxyAgnObj, StarObj
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogSersic2D, PhoSimCatalogPoint, \
                                                         PhoSimCatalogZPoint
from lsst.sims.catalogs.generation.utils import makePhoSimTestDB

class SearchReversion(CatalogDBObject):
    """
    This is a mixin which is used below to force the galaxy CatalogDBObjects created for
    this unittest to use the methods defined in the CatalogDBObject class.  This is because
    we are using classes which inherit from GalaxyTileObj but do not actually want to use
    the tiled query routines.

    We also use this mixin for our stellar database object.  This is because StarObj
    implements a query search based on htmid, which the test database for this unit
    test will not have.
    """

    def _get_column_query(self, *args, **kwargs):
        return CatalogDBObject._get_column_query(self,*args, **kwargs)

    def _final_pass(self, *args, **kwargs):
        return CatalogDBObject._final_pass(self,*args, **kwargs)

    def query_columns(self, *args, **kwargs):
        return CatalogDBObject.query_columns(self, *args, **kwargs)
     
class testGalaxyBulge(SearchReversion, GalaxyBulgeObj):
    """
    A class for storing galaxy bulges
    """
    objid = 'phoSimTestBulges'
    objectTypeId = 88

    #The code below makes sure that we can store RA, Dec in degrees
    #in the database but use radians in our calculations.
    #We had to overwrite the original columns list because
    #GalaxyTileObject daughter classes assume that RA and Dec are stored
    #in radians in the database.  This is a side effect of the tiling
    #scheme used to cover the whole sky.

    columns = GalaxyBulgeObj.columns
    _to_remove = []
    for entry in columns:
        if entry[0] == 'raJ2000' or entry[0] == 'decJ2000':
            _to_remove.append(entry)
    for target in _to_remove:
        columns.remove(target)

    columns.append(('raJ2000','ra*PI()/180.'))
    columns.append(('decJ2000','dec*PI()/180.'))

class testGalaxyDisk(SearchReversion, GalaxyDiskObj):
    objid = 'phoSimTestDisks'
    objectTypeId = 89

    #The code below makes sure that we can store RA, Dec in degrees
    #in the database but use radians in our calculations.
    #We had to overwrite the original columns list because
    #GalaxyTileObject daughter classes assume that RA and Dec are stored
    #in radians in the database.  This is a side effect of the tiling
    #scheme used to cover the whole sky.

    columns = GalaxyDiskObj.columns
    _to_remove = []
    for entry in columns:
        if entry[0] == 'raJ2000' or entry[0] == 'decJ2000':
            _to_remove.append(entry)
    for target in _to_remove:
        columns.remove(target)

    columns.append(('raJ2000','ra*PI()/180.'))
    columns.append(('decJ2000','dec*PI()/180.'))

class testGalaxyAgn(SearchReversion, GalaxyAgnObj):
    objid = 'phoSimTestAgn'
    objectTypeId = 90

    #The code below makes sure that we can store RA, Dec in degrees
    #in the database but use radians in our calculations.
    #We had to overwrite the original columns list because
    #GalaxyTileObject daughter classes assume that RA and Dec are stored
    #in radians in the database.  This is a side effect of the tiling
    #scheme used to cover the whole sky.

    columns = GalaxyAgnObj.columns
    _to_remove = []
    for entry in columns:
        if entry[0] == 'raJ2000' or entry[0] == 'decJ2000':
            _to_remove.append(entry)
    for target in _to_remove:
        columns.remove(target)

    columns.append(('raJ2000','ra*PI()/180.'))
    columns.append(('decJ2000','dec*PI()/180.'))

class testStars(SearchReversion, StarObj):
    objid = 'phoSimTestStars'
    objectTypeId = 91

    #The code below removes the definitions of galacticAv and magNorm
    #from this database object.  The definitions of those columns which
    #are implemented in StarObj rely on mathematical functions which
    #are not defined in sqlite.

    columns = StarObj.columns
    _to_remove = []
    for entry in columns:
        if entry[0] == 'galacticAv':
            _to_remove.append(entry)
        elif entry[0] == 'magNorm':
            _to_remove.append(entry)

    for target in _to_remove:
        columns.remove(target)

class PhoSimCatalogTest(unittest.TestCase):

    def setUp(self):
        self.obs_metadata = makePhoSimTestDB(size=10)
        self.bulgeDB = testGalaxyBulge(address='sqlite:///PhoSimTestDatabase.db')
        self.diskDB = testGalaxyDisk(address='sqlite:///PhoSimTestDatabase.db')
        self.agnDB = testGalaxyAgn(address='sqlite:///PhoSimTestDatabase.db')
        self.starDB = testStars(address='sqlite:///PhoSimTestDatabase.db')
        baseLineFileName = os.getenv('SIMS_CATUTILS_DIR')+'/tests/testData/phoSimControlCatalog.txt'
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

        catName = 'phoSimTestCatalog.txt'
        testBulge.write_catalog(catName)
        testDisk.write_catalog(catName, write_header=False, write_mode='a')
        testAgn.write_catalog(catName, write_header=False, write_mode='a')
        testStar.write_catalog(catName, write_header=False, write_mode='a')

        testFile = open(catName,'r')
        testLines = testFile.readlines()
        controlLines = self.baseLineFile.readlines()
        for line in testLines:
            self.assertTrue(line in controlLines)
        testFile.close()

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
