import os
import unittest
import lsst.utils.tests as utilsTests
import numpy
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catUtils.utils import testStarsDBObj, testGalaxyDiskDBObj, \
                                     testGalaxyBulgeDBObj, testGalaxyAgnDBObj
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogSersic2D, PhoSimCatalogPoint, \
                                                         PhoSimCatalogZPoint
from lsst.sims.catalogs.generation.utils import makePhoSimTestDB
from lsst.sims.catalogs.generation.db import ObservationMetaData
from lsst.sims.photUtils import VariabilityStars, VariabilityGalaxies
from lsst.sims.photUtils.utils import TestVariabilityMixin
from lsst.sims.coordUtils import AstrometryStars, AstrometryGalaxies


class PhoSimPointVariable(PhoSimCatalogPoint, VariabilityStars, TestVariabilityMixin):
    pass


class PhoSimZPointVariable(PhoSimCatalogZPoint, VariabilityStars, TestVariabilityMixin):
    pass


class ControlCatalog(InstanceCatalog):
    transformations = {'raPhoSim':numpy.degrees, 'decPhoSim':numpy.degrees}


class AgnControlCatalog(ControlCatalog, VariabilityGalaxies, TestVariabilityMixin, AstrometryGalaxies):
    column_outputs = ['raPhoSim', 'decPhoSim', 'magNorm', 'delta_rAgn']


class BulgeControlCatalog(ControlCatalog, AstrometryGalaxies):
    column_outputs = ['raPhoSim', 'decPhoSim', 'magNorm']


class DiskControlCatalog(ControlCatalog, AstrometryGalaxies):
    column_outputs = ['raPhoSim', 'decPhoSim', 'magNorm']


class StarControlCatalog(ControlCatalog, AstrometryStars, VariabilityStars, TestVariabilityMixin):
    column_outputs = ['raPhoSim', 'decPhoSim', 'magNorm', 'delta_lsst_r']


class PhoSimVariabilityTest(unittest.TestCase):

    @classmethod
    def setUp(cls):
        cls.dbName = 'PhoSimVariabilityDatabase.db'
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)

        cls.obs_metadata = makePhoSimTestDB(size=10, filename=cls.dbName)

        connection = 'sqlite:///' + cls.dbName

        cls.bulgeDB = testGalaxyBulgeDBObj(address=connection)
        cls.diskDB = testGalaxyDiskDBObj(address=connection)
        cls.agnDB = testGalaxyAgnDBObj(address=connection)
        cls.starDB = testStarsDBObj(address=connection)

    @classmethod
    def tearDownClass(cls):
        if os.path.exists(cls.dbName):
            os.unlink(cls.dbName)

    def testAgn(self):
        baseline = AgnControlCatalog(self.agnDB, obs_metadata=self.obs_metadata)
        test = PhoSimZPointVariable(self.agnDB, obs_metadata=self.obs_metadata)
        
        for bb, tt in zip(baseline.iter_catalog(), test.iter_catalog()):
            self.assertEqual(tt[2], bb[0])
            self.assertEqual(tt[3], bb[1])
            self.assertAlmostEqual(bb[2] + bb[3], tt[4], 10)
            self.assertTrue(numpy.abs(bb[3]) > 0.0)

    def testStars(self):
        baseline = StarControlCatalog(self.starDB, obs_metadata=self.obs_metadata)
        test = PhoSimPointVariable(self.starDB, obs_metadata=self.obs_metadata)
        
        for bb, tt in zip(baseline.iter_catalog(), test.iter_catalog()):
            self.assertEqual(tt[2], bb[0])
            self.assertEqual(tt[3], bb[1])
            self.assertAlmostEqual(bb[2] + bb[3], tt[4], 10)
            self.assertTrue(numpy.abs(bb[3]) > 0.0)        

    def testBulges(self):
        baseline = BulgeControlCatalog(self.bulgeDB, obs_metadata=self.obs_metadata)
        test = PhoSimCatalogSersic2D(self.bulgeDB, obs_metadata=self.obs_metadata)

        for bb, tt in zip(baseline.iter_catalog(), test.iter_catalog()):
            self.assertEqual(tt[2], bb[0])
            self.assertEqual(tt[3], bb[1])
            self.assertAlmostEqual(bb[2], tt[4], 10)


    def testDisks(self):
        baseline = DiskControlCatalog(self.diskDB, obs_metadata=self.obs_metadata)
        test = PhoSimCatalogSersic2D(self.diskDB, obs_metadata=self.obs_metadata)

        for bb, tt in zip(baseline.iter_catalog(), test.iter_catalog()):
            self.assertEqual(tt[2], bb[0])
            self.assertEqual(tt[3], bb[1])
            self.assertAlmostEqual(bb[2], tt[4], 10)


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(PhoSimVariabilityTest)

    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)
if __name__ == "__main__":
    run(True)
