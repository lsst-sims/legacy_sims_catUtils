import os
import numpy
import unittest
import lsst.utils.tests as utilsTests
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.generation.utils import makePhoSimTestDB
from lsst.sims.catUtils.galSimInterface import GalSimGalaxies
from lsst.sims.catUtils.utils import testGalaxyBulge

class testCat(InstanceCatalog):
    column_outputs = ['raJ2000', 'decJ2000']

class GalSimInterfaceTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.dbName = 'galSimTestDB.db'
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)
        
        displacedRA = numpy.array([72.0/3600.0])
        displacedDec = numpy.array([0.0])
        cls.obs_metadata = makePhoSimTestDB(filename=cls.dbName, size=1,
                                            displacedRA=displacedRA, displacedDec=displacedDec)
        cls.connectionString = 'sqlite:///'+cls.dbName

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)

        del cls.dbName
        del cls.connectionString
        del cls.obs_metadata

    def testGalaxyBulge(self):
        catName = 'galaxyBulgeCatalog.sav'
        gals = testGalaxyBulge(address=self.connectionString)
        cat = GalSimGalaxies(gals, obs_metadata = self.obs_metadata)
        cat.write_catalog(catName)
        
        #right now the problem is that, without CATSIM-179, the obs_metadata returned
        #is wrong and the catalog is empty

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(GalSimInterfaceTest)

    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)
if __name__ == "__main__":
    run(True)
