from __future__ import with_statement
import os
import copy
import numpy
import unittest
import eups
import galsim
import lsst.utils.tests as utilsTests
from lsst.sims.photUtils import Bandpass
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.generation.utils import makePhoSimTestDB
from lsst.sims.catUtils.galSimInterface import GalSimGalaxies, GalSimStars, GalSimAgn, \
                                               ExampleGaussianPSF, ExampleOpticalPSF
from lsst.sims.catUtils.utils import calcADUwrapper, testGalaxyBulgeDBObj, testGalaxyDiskDBObj, \
                                     testGalaxyAgnDBObj, testStarsDBObj
import lsst.afw.image as afwImage

class testGalaxyCatalog(GalSimGalaxies):
    column_outputs = copy.deepcopy(GalSimGalaxies.column_outputs)
    column_outputs.remove('fitsFiles')
    column_outputs.append('magNorm')
    column_outputs.append('redshift')
    column_outputs.append('internalAv')
    column_outputs.append('internalRv')
    column_outputs.append('galacticAv')
    column_outputs.append('galacticRv')
    column_outputs.append('fitsFiles')

class testStarCatalog(GalSimStars):
    column_outputs = copy.deepcopy(GalSimStars.column_outputs)
    column_outputs.remove('fitsFiles')
    column_outputs.append('magNorm')
    column_outputs.append('redshift')
    column_outputs.append('internalAv')
    column_outputs.append('internalRv')
    column_outputs.append('galacticAv')
    column_outputs.append('galacticRv')
    column_outputs.append('fitsFiles')
    
    PSF = ExampleOpticalPSF()

class testAgnCatalog(GalSimAgn):
    column_outputs = copy.deepcopy(GalSimAgn.column_outputs)
    column_outputs.remove('fitsFiles')
    column_outputs.append('magNorm')
    column_outputs.append('redshift')
    column_outputs.append('internalAv')
    column_outputs.append('internalRv')
    column_outputs.append('galacticAv')
    column_outputs.append('galacticRv')
    column_outputs.append('fitsFiles')
    
    PSF = ExampleGaussianPSF()

class psfCatalog(testGalaxyCatalog):
    PSF = ExampleGaussianPSF()
    
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
    
    def catalogTester(self, catName=None, catalog=None, nameRoot=None):

        catalog.write_images(nameRoot=nameRoot)
        
        galsimCounts = {}
        controlCounts = {}
        
        for name in catalog.galSimInterpreter.detectorImages:
            if nameRoot is not None:
                name = nameRoot+'_'+name
            im = afwImage.ImageF(name)
            imArr = im.getArray()
            galsimCounts[name[-6]] = imArr.sum()
            controlCounts[name[-6]] = 0.0
            os.unlink(name)
            
        with open(catName, 'r') as testFile:
            lines = testFile.readlines()
            for line in lines:
                if line[0] != '#':
                    gg = line.split(';')
                    sedName = gg[5]
                    magNorm = float(gg[11])
                    redshift = float(gg[12])
                    internalAv = float(gg[13])
                    internalRv = float(gg[14])
                    galacticAv = float(gg[15])
                    galacticRv = float(gg[16])
        

                    for filterName in controlCounts:
                        bandPassName=os.path.join(eups.productDir('throughputs'),'baseline',('total_'+filterName+'.dat'))
                        bandpass = Bandpass()
                        bandpass.readThroughput(bandPassName)
                        controlCounts[filterName] += calcADUwrapper(sedName=sedName, bandpass=bandpass,
                                                                    redshift=redshift, magNorm=magNorm,
                                                                    internalAv=internalAv, internalRv=internalRv, 
                                                                    galacticAv=galacticAv, galacticRv=galacticRv)
               
            drawnFilters = 0
            for ff in controlCounts:
                if controlCounts[ff] > 1000.0:
                    drawnFilters += 1
                    self.assertTrue(numpy.abs(controlCounts[ff] - galsimCounts[ff]) < 0.05*controlCounts[ff])
                    
            self.assertTrue(drawnFilters>4)

    def testGalaxyBulges(self):
        catName = 'testBulgeCat.sav'
        gals = testGalaxyBulgeDBObj(address=self.connectionString)
        cat = testGalaxyCatalog(gals, obs_metadata = self.obs_metadata)
        cat.write_catalog(catName)
        self.catalogTester(catName=catName, catalog=cat, nameRoot='bulge')
        if os.path.exists(catName):
            os.unlink(catName)
       
    def testGalaxyDisks(self):
        catName = 'testDiskCat.sav'
        gals = testGalaxyDiskDBObj(address=self.connectionString)
        cat = testGalaxyCatalog(gals, obs_metadata = self.obs_metadata)
        cat.write_catalog(catName)
        self.catalogTester(catName=catName, catalog=cat, nameRoot='disk')
        if os.path.exists(catName):
            os.unlink(catName)
    
    def testStars(self):
        catName = 'testStarCat.sav'
        stars = testStarsDBObj(address=self.connectionString)
        cat = testStarCatalog(stars, obs_metadata = self.obs_metadata)
        cat.write_catalog(catName)
        self.catalogTester(catName=catName, catalog=cat, nameRoot='stars')
        if os.path.exists(catName):
            os.unlink(catName)

    def testAgns(self):
        catName = 'testAgnCat.sav'
        agn = testGalaxyAgnDBObj(address=self.connectionString)
        cat = testAgnCatalog(agn, obs_metadata = self.obs_metadata)
        cat.write_catalog(catName)
        self.catalogTester(catName=catName, catalog=cat, nameRoot='agn')
        if os.path.exists(catName):
            os.unlink(catName)
  
    
    def testPSFimages(self):
        catName = 'testPSFcat.sav'
        gals = testGalaxyBulgeDBObj(address=self.connectionString)
        cat = psfCatalog(gals, obs_metadata = self.obs_metadata)
        cat.write_catalog(catName)
        self.catalogTester(catName=catName, catalog=cat, nameRoot='psf')
        if os.path.exists(catName):
            os.unlink(catName)

    def testMultipleImages(self):
        dbName = 'galSimTestMultipleDB.db'
        if os.path.exists(dbName):
            os.unlink(dbName)
        
        displacedRA = numpy.array([72.0/3600.0, 50.0/3600.0, 75.0/3600.0])
        displacedDec = numpy.array([0.0, 15.0/3600.0, -15.0/3600.0])
        obs_metadata = makePhoSimTestDB(filename=dbName, size=1,
                                            displacedRA=displacedRA, displacedDec=displacedDec)
        connectionString = 'sqlite:///'+dbName
        gals = testGalaxyBulgeDBObj(address=connectionString)
        cat = testGalaxyCatalog(gals, obs_metadata=obs_metadata)
        catName = 'multipleCatalog.sav'
        cat.write_catalog(catName)
        self.catalogTester(catName=catName, catalog=cat, nameRoot='multiple')
        if os.path.exists(catName):
            os.unlink(catName)
        
        stars = testStarsDBObj(address=connectionString)
        cat = testStarCatalog(stars, obs_metadata=obs_metadata)
        catName = 'multipleStarCatalog.sav'
        cat.write_catalog(catName)
        self.catalogTester(catName=catName, catalog=cat, nameRoot='multipleStars')
        if os.path.exists(catName):
            os.unlink(catName)
        
        if os.path.exists(dbName):
            os.unlink(dbName)

    def testCompoundFitsFiles(self):
        dbName1 = 'galSimTestCompound1DB.db'
        if os.path.exists(dbName1):
            os.unlink(dbName1)
        
        displacedRA = numpy.array([72.0/3600.0, 50.0/3600.0, 75.0/3600.0])
        displacedDec = numpy.array([0.0, 15.0/3600.0, -15.0/3600.0])
        obs_metadata1 = makePhoSimTestDB(filename=dbName1, size=1,
                                            displacedRA=displacedRA, displacedDec=displacedDec)
        connectionString1 = 'sqlite:///'+dbName1
        
        dbName2 = 'galSimTestCompound2DB.db'
        if os.path.exists(dbName2):
            os.unlink(dbName2)
        
        displacedRA = numpy.array([50.0/3600.0, 60.0/3600.0, 65.0/3600.0])
        displacedDec = numpy.array([-5.0/3600.0, 10.0/3600.0, 10.0/3600.0])
        obs_metadata2 = makePhoSimTestDB(filename=dbName2, size=1,
                                            displacedRA=displacedRA, displacedDec=displacedDec)
        connectionString2 = 'sqlite:///'+dbName2
        
        
        gals = testGalaxyBulgeDBObj(address=connectionString1)
        cat1 = testGalaxyCatalog(gals, obs_metadata=obs_metadata1)
        catName = 'compoundCatalog.sav'
        cat1.write_catalog(catName)

        stars = testStarsDBObj(address=connectionString2)
        cat2 = testStarCatalog(stars, obs_metadata=obs_metadata2)
        cat2.copyGalSimInterpreter(cat1)
        cat2.write_catalog(catName, write_header=False, write_mode='a')
        self.catalogTester(catName=catName, catalog=cat2, nameRoot='compound')
        
        if os.path.exists(dbName1):
            os.unlink(dbName1)
        if os.path.exists(dbName2):
            os.unlink(dbName2)
        if os.path.exists(catName):
            os.unlink(catName)
    

def suite():
    utilsTests.init()
    suites = []
    suites += unittest.makeSuite(GalSimInterfaceTest)

    return unittest.TestSuite(suites)

def run(shouldExit = False):
    utilsTests.run(suite(), shouldExit)
if __name__ == "__main__":
    run(True)
