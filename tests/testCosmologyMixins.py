import os
import unittest
import lsst.utils.tests
import numpy as np

from lsst.sims.photUtils import CosmologyObject

from lsst.sims.catalogs.utils import myTestGals, makeGalTestDB

from lsst.sims.catUtils.utils import testGalaxies
from lsst.sims.catUtils.mixins import CosmologyMixin


def setup_module(module):
    lsst.utils.tests.init()


class cosmologicalGalaxyCatalog(testGalaxies, CosmologyMixin):
    catalog_type = __file__ + 'cosmo_galaxy_catalog'
    column_outputs = ['galid', 'lsst_u', 'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y',
                      'uBulge', 'gBulge', 'rBulge', 'iBulge', 'zBulge', 'yBulge',
                      'uDisk', 'gDisk', 'rDisk', 'iDisk', 'zDisk', 'yDisk',
                      'uAgn', 'gAgn', 'rAgn', 'iAgn', 'zAgn', 'yAgn',
                      'redshift', 'cosmologicalDistanceModulus']


class absoluteGalaxyCatalog(testGalaxies):
    catalog_type = __file__ + 'abs_galaxy_catalog'
    column_outputs = ['galid', 'lsst_u', 'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y',
                      'uBulge', 'gBulge', 'rBulge', 'iBulge', 'zBulge', 'yBulge',
                      'uDisk', 'gDisk', 'rDisk', 'iDisk', 'zDisk', 'yDisk',
                      'uAgn', 'gAgn', 'rAgn', 'iAgn', 'zAgn', 'yAgn',
                      'redshift']

    def get_cosmologicalDistanceModulus(self):
        """
        Must set this to zero rather than `None` so that PhotometryGalaxies
        does not apply cosmological dimming
        """
        return np.zeros(len(self.column_by_name('galid')))


class CosmologyMixinUnitTest(unittest.TestCase):
    """
    This class will test to make sure that our example CosmologyMixin
    (defined in lsst/sims/photUtils/examples/CosmologyMixin.py)
    can produce a catalog
    """

    @classmethod
    def setUpClass(cls):
        cls.dbName = 'cosmologyTestDB.db'
        cls.dbSize = 100
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)
        makeGalTestDB(size=cls.dbSize, seedVal=1, filename=cls.dbName)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)

        del cls.dbName
        del cls.dbSize

    def setUp(self):
        self.catName = 'cosmologyCatalog.txt'
        if os.path.exists(self.catName):
            os.unlink(self.catName)

    def tearDown(self):
        if os.path.exists(self.catName):
            os.unlink(self.catName)
        del self.catName

    def testCosmologyCatalog(self):
        """
        Does a catalog get written?
        """
        dbObj = myTestGals(database=self.dbName)
        cat = cosmologicalGalaxyCatalog(dbObj)
        cat.write_catalog(self.catName)

    def testCatalogDistanceModulus(self):
        """
        Does cosmologicalDistanceModulus get properly applied
        """
        dbObj = myTestGals(database=self.dbName)
        cosmoCat = cosmologicalGalaxyCatalog(dbObj)
        controlCat = absoluteGalaxyCatalog(dbObj)
        cosmoIter = cosmoCat.iter_catalog(chunk_size=self.dbSize)
        controlIter = controlCat.iter_catalog(chunk_size=self.dbSize)

        cosmology = CosmologyObject()

        for (cosmoRow, controlRow) in zip(cosmoIter, controlIter):
            modulus = cosmology.distanceModulus(controlRow[25])
            self.assertEqual(cosmoRow[0], controlRow[0])
            self.assertEqual(cosmoRow[25], controlRow[25])
            self.assertEqual(cosmoRow[26], modulus)
            for i in range(1, 25):
                self.assertAlmostEqual(cosmoRow[i], controlRow[i] + modulus, 6)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
