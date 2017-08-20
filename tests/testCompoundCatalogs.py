import unittest
import os
import shutil
import tempfile
import lsst.utils.tests
from lsst.utils import getPackageDir

from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import (GalaxyBulgeObj, GalaxyDiskObj,
                                                  GalaxyAgnObj, GalaxyTileCompoundObj,
                                                  StarObj)

from lsst.sims.catalogs.definitions import InstanceCatalog, CompoundInstanceCatalog

_testCompoundCatalogs_is_connected = True
try:
    _example_db = GalaxyBulgeObj()
except:
    _testCompoundCatalogs_is_connected = False


def setup_module(module):
    lsst.utils.tests.init()


class BulgeDiskCatalog(InstanceCatalog):
    catalog_type = __file__ + 'bulge_disk_catalog'
    cannot_be_null = ['sedFilename']
    column_outputs = ['galtileid', 'raJ2000', 'decJ2000',
                      'componentra', 'componentdec',
                      'magNorm', 'sedFilename',
                      'majorAxis', 'minorAxis',
                      'positionAngle',
                      'halfLightRadius',
                      'internalExtinctionModel',
                      'internalAv', 'internalRv']


class AgnCatalog(InstanceCatalog):
    catalog_type = __file__ + 'agn_catalog'
    cannot_be_null = ['sedFilename']
    column_outputs = ['galtileid', 'raJ2000', 'decJ2000',
                      'componentra', 'componentdec',
                      'magNorm', 'sedFilename',
                      'variabilityParameters']


class StarCatalog(InstanceCatalog):
    catalog_type = __file__ + 'star_catalog'
    cannot_be_null = ['sedFilename']
    column_outputs = ['id', 'raJ2000', 'decJ2000',
                      'glon', 'glat', 'magNorm',
                      'properMotionRa', 'properMotionDec',
                      'parallax', 'galacticAv', 'radialVelocity',
                      'variabilityParameters', 'sedFilename']


class CompoundCatalogTest(unittest.TestCase):

    def setUp(self):
        self.baseDir = tempfile.mkdtemp(dir=ROOT, prefix='compoundCatalogTest-')

    def tearDown(self):
        if os.path.exists(self.baseDir):
            shutil.rmtree(self.baseDir)

    @unittest.skipIf(not _testCompoundCatalogs_is_connected,
                     "We are not connected to fatboy")
    def testGalaxyCatalog(self):
        """
        Test GalaxyTileCompoundObj by creating a catalog of galaxy bulges, disks,
        and agns using both the 'old fashioned way' (one catalog at a time), and
        using CompoundInstanceCatalog
        """
        controlFileName = os.path.join(self.baseDir, 'gal_compound_control.txt')
        testFileName = os.path.join(self.baseDir, 'gal_compound_test.txt')

        if os.path.exists(controlFileName):
            os.unlink(controlFileName)
        if os.path.exists(testFileName):
            os.unlink(testFileName)

        obs = ObservationMetaData(pointingRA=25.0, pointingDec=-45.0,
                                  boundType='circle', boundLength=0.05)

        dbBulge = GalaxyBulgeObj()
        dbDisk = GalaxyDiskObj()
        dbAgn = GalaxyAgnObj()

        catBulge = BulgeDiskCatalog(dbBulge, obs_metadata=obs)
        catDisk = BulgeDiskCatalog(dbDisk, obs_metadata=obs)
        catAgn = AgnCatalog(dbAgn, obs_metadata=obs)

        catBulge.write_catalog(controlFileName, write_header=False, chunk_size=10000)
        catDisk.write_catalog(controlFileName, write_mode='a', write_header=False, chunk_size=10000)
        catAgn.write_catalog(controlFileName, write_mode='a', write_header=False, chunk_size=10000)

        totalCat = CompoundInstanceCatalog([BulgeDiskCatalog, BulgeDiskCatalog, AgnCatalog],
                                           [GalaxyDiskObj, GalaxyBulgeObj, GalaxyAgnObj],
                                           obs_metadata=obs,
                                           compoundDBclass=GalaxyTileCompoundObj)

        totalCat.write_catalog(testFileName, write_header=False, chunk_size=10000)

        with open(controlFileName, 'r') as controlFile:
            control = controlFile.readlines()
        controlFile.close()

        with open(testFileName, 'r') as testFile:
            test = testFile.readlines()

        for line in control:
            self.assertIn(line, test)

        for line in test:
            self.assertIn(line, control)

    @unittest.skipIf(not _testCompoundCatalogs_is_connected,
                     "We are not connected to fatboy")
    def testGalaxyAndStarCatalog(self):
        """
        Test GalaxyTileCompoundObj by creating a catalog of galaxy bulges, disks,
        agns, and stars using both the 'old fashioned way' (one catalog at a time), and
        using CompoundInstanceCatalog
        """
        controlFileName = os.path.join(self.baseDir, 'galStar_compound_control.txt')
        testFileName = os.path.join(self.baseDir, 'galStar_compound_test.txt')

        if os.path.exists(controlFileName):
            os.unlink(controlFileName)
        if os.path.exists(testFileName):
            os.unlink(testFileName)

        obs = ObservationMetaData(pointingRA=25.0, pointingDec=-45.0,
                                  boundType='circle', boundLength=0.05)

        dbBulge = GalaxyBulgeObj()
        dbDisk = GalaxyDiskObj()
        dbAgn = GalaxyAgnObj()
        dbStar = StarObj()

        catBulge = BulgeDiskCatalog(dbBulge, obs_metadata=obs)
        catDisk = BulgeDiskCatalog(dbDisk, obs_metadata=obs)
        catAgn = AgnCatalog(dbAgn, obs_metadata=obs)
        catStar = StarCatalog(dbStar, obs_metadata=obs)

        catBulge.write_catalog(controlFileName, write_header=False, chunk_size=10000)
        catDisk.write_catalog(controlFileName, write_mode='a', write_header=False, chunk_size=10000)
        catAgn.write_catalog(controlFileName, write_mode='a', write_header=False, chunk_size=10000)
        catStar.write_catalog(controlFileName, write_mode='a', write_header=False, chunk_size=10000)

        totalCat = CompoundInstanceCatalog([BulgeDiskCatalog, BulgeDiskCatalog, StarCatalog, AgnCatalog],
                                           [GalaxyBulgeObj, GalaxyDiskObj, StarObj, GalaxyAgnObj],
                                           obs_metadata=obs,
                                           compoundDBclass=GalaxyTileCompoundObj)

        totalCat.write_catalog(testFileName, write_header=False, chunk_size=10000)

        with open(controlFileName, 'r') as controlFile:
            control = controlFile.readlines()

        with open(testFileName, 'r') as testFile:
            test = testFile.readlines()

        for line in control:
            self.assertIn(line, test)

        for line in test:
            self.assertIn(line, control)


class MemoryTestClass(lsst.utils.tests.MemoryTestCase):
    pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
