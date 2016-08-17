import unittest
import os
import lsst.utils.tests as utilsTests
from lsst.utils import getPackageDir

from lsst.sims.utils import ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import GalaxyBulgeObj, GalaxyDiskObj, \
                                                 GalaxyAgnObj, GalaxyTileCompoundObj, \
                                                 StarObj

from lsst.sims.catalogs.definitions import InstanceCatalog, CompoundInstanceCatalog

_is_connected = True
try:
    _example_db = GalaxyBulgeObj()
except:
    _is_connected = False


class BulgeDiskCatalog(InstanceCatalog):
    cannot_be_null = ['sedFilename']
    column_outputs = ['galtileid', 'raJ2000', 'decJ2000',
                     'componentra', 'componentdec',
                     'magNorm', 'sedFilename',
                     'majorAxis', 'minorAxis',
                     'positionAngle',
                     'halfLightRadius',
                     'internalExtinctionModel',
                     'internalAv', 'internalRv',
                     ]


class AgnCatalog(InstanceCatalog):
    cannot_be_null = ['sedFilename']
    column_outputs = ['galtileid', 'raJ2000', 'decJ2000',
                      'componentra','componentdec',
                      'magNorm', 'sedFilename',
                      'variabilityParameters',
                      ]


class StarCatalog(InstanceCatalog):
    cannot_be_null = ['sedFilename']
    column_outputs = ['id', 'raJ2000', 'decJ2000',
                      'glon', 'glat', 'magNorm',
                      'properMotionRa', 'properMotionDec',
                      'parallax', 'galacticAv', 'radialVelocity',
                      'variabilityParameters', 'sedFilename']

class CompoundCatalogTest(unittest.TestCase):

    def setUp(self):
        self.baseDir = os.path.join(getPackageDir('sims_catUtils'), \
                                    'tests', 'scratchSpace')

    @unittest.skipIf(not _is_connected, "We are not connected to fatboy")
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

        controlFile = open(controlFileName, 'r')
        control = controlFile.readlines()
        controlFile.close()

        testFile = open(testFileName, 'r')
        test = testFile.readlines()
        testFile.close()

        for line in control:
            self.assertTrue(line in test)

        for line in test:
            self.assertTrue(line in control)

        if os.path.exists(controlFileName):
            os.unlink(controlFileName)

        if os.path.exists(testFileName):
            os.unlink(testFileName)

    @unittest.skipIf(not _is_connected, "We are not connected to fatboy")
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

        controlFile = open(controlFileName, 'r')
        control = controlFile.readlines()
        controlFile.close()

        testFile = open(testFileName, 'r')
        test = testFile.readlines()
        testFile.close()

        for line in control:
            self.assertTrue(line in test)

        for line in test:
            self.assertTrue(line in control)

        if os.path.exists(controlFileName):
            os.unlink(controlFileName)

        if os.path.exists(testFileName):
            os.unlink(testFileName)


def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(CompoundCatalogTest)

    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)
if __name__ == "__main__":
    run(True)
