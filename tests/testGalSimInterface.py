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
    """
    Wraps the GalSimGalaxies class.  Adds columns to the output
    so that we can read the InstanceCatalog back in and verify that
    GalSim put the correct number of ADU in each FITS file.
    """
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
    """
    Wraps the GalSimStars class.  Adds columns to the output
    so that we can read the InstanceCatalog back in and verify that
    GalSim put the correct number of ADU in each FITS file.
    """
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
    """
    Wraps the GalSimAgn class.  Adds columns to the output
    so that we can read the InstanceCatalog back in and verify that
    GalSim put the correct number of ADU in each FITS file.
    """
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
    """
    Adds a PSF to testGalaxyCatalog
    """
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
        """
        Reads in a GalSim Instance Catalog.  Writes the images from that catalog.
        Then reads those images back in.  Uses AFW to calculate the number of counts
        in each FITS image.  Reads in the InstanceCatalog associated with those images.
        Uses sims_photUtils code to calculate the ADU for each object on the FITS images.
        Verifies that the two independent calculations of counts agree (to within a tolerance,
        since the GalSim images are generated in a pseudo-random way).

        @param [in] catName is the name of the InstanceCatalog that has been written to disk

        @paranm [in] catalog is the actual InstanceCatalog instantiation

        @param [in] nameRoot is a string appended to the names of the FITS files being written
        """

        #write the fits files
        catalog.write_images(nameRoot=nameRoot)

        #a dictionary of ADU for each FITS file as calculated by GalSim
        #(indexed on the name of the FITS file)
        galsimCounts = {}

        #a dictionary of ADU for each FITS file as calculated by sims_photUtils
        #(indexed on the name of the FITS file)
        controlCounts = {}

        #a list of bandpasses over which we are integraging
        listOfFilters = []

        #read in the names of all of the written fits files directly from the
        #InstanceCatalog's GalSimInterpreter
        #Use AFW to read in the FITS files and calculate the ADU
        for name in catalog.galSimInterpreter.detectorImages:
            if nameRoot is not None:
                name = nameRoot+'_'+name
            im = afwImage.ImageF(name)
            imArr = im.getArray()
            galsimCounts[name] = imArr.sum()
            controlCounts[name] = 0.0

            if name[-6] not in listOfFilters:
                listOfFilters.append(name[-6])

            os.unlink(name)

        #Read in the InstanceCatalog.  For each object in the catalog, use sims_photUtils
        #to calculate the ADU.  Keep track of how many ADU should be in each FITS file.
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
                    listOfFileNames = gg[17].split('//')
                    alreadyWritten = []

                    for name in listOfFileNames:

                        #guard against objects being written on one
                        #chip more than once
                        msg = '%s was written on %s more than once' % (sedName, name)
                        self.assertTrue(name not in alreadyWritten, msg=msg)
                        alreadyWritten.append(name)

                        #loop over all of the detectors on which an object fell
                        #(this is not a terribly great idea, since our conservative implementation
                        #of GalSimInterpreter._doesObjectImpingeOnDetector means that some detectors
                        #will be listed here even though the object does not illumine them)
                        for filterName in listOfFilters:
                            chipName = name.replace(':','_')
                            chipName = chipName.replace(' ','_')
                            chipName = chipName.replace(',','_')
                            chipName = chipName.strip()

                            fullName = nameRoot+'_'+chipName+'_'+filterName+'.fits'

                            bandpassName=os.path.join(eups.productDir('throughputs'),'baseline',('total_'+filterName+'.dat'))
                            bandpass = Bandpass()
                            bandpass.readThroughput(bandpassName)
                            controlCounts[fullName] += calcADUwrapper(sedName=sedName, bandpass=bandpass,
                                                                        redshift=redshift, magNorm=magNorm,
                                                                        internalAv=internalAv, internalRv=internalRv,
                                                                        galacticAv=galacticAv, galacticRv=galacticRv)

            drawnDetectors = 0
            unDrawnDetectors = 0
            for ff in controlCounts:
                if controlCounts[ff] > 1000.0 and galsimCounts[ff] > 0.001:
                    #because, for really dim images, there could be enough statistical imprecision in the GalSim drawing routine
                    #to violate the condition below
                    drawnDetectors += 1
                    self.assertTrue(numpy.abs(controlCounts[ff] - galsimCounts[ff]) < 0.05*controlCounts[ff])
                elif galsimCounts[ff] > 0.001:
                    unDrawnDetectors += 1

            #to make sure we did not neglect more than one detector
            self.assertTrue(unDrawnDetectors<2)
            self.assertTrue(drawnDetectors>0)

    def testGalaxyBulges(self):
        """
        Test that GalSimInterpreter puts the right number of counts on images of galaxy bulges
        """
        catName = 'testBulgeCat.sav'
        gals = testGalaxyBulgeDBObj(address=self.connectionString)
        cat = testGalaxyCatalog(gals, obs_metadata = self.obs_metadata)
        cat.write_catalog(catName)
        self.catalogTester(catName=catName, catalog=cat, nameRoot='bulge')
        if os.path.exists(catName):
            os.unlink(catName)

    def testGalaxyDisks(self):
        """
        Test that GalSimInterpreter puts the right number of counts on images of galaxy disks
        """
        catName = 'testDiskCat.sav'
        gals = testGalaxyDiskDBObj(address=self.connectionString)
        cat = testGalaxyCatalog(gals, obs_metadata = self.obs_metadata)
        cat.write_catalog(catName)
        self.catalogTester(catName=catName, catalog=cat, nameRoot='disk')
        if os.path.exists(catName):
            os.unlink(catName)

    def testStars(self):
        """
        Test that GalSimInterpreter puts the right number of counts on images of stars
        """
        catName = 'testStarCat.sav'
        stars = testStarsDBObj(address=self.connectionString)
        cat = testStarCatalog(stars, obs_metadata = self.obs_metadata)
        cat.write_catalog(catName)
        self.catalogTester(catName=catName, catalog=cat, nameRoot='stars')
        if os.path.exists(catName):
            os.unlink(catName)

    def testAgns(self):
        """
        Test that GalSimInterpreter puts the right number of counts on images of AGN
        """
        catName = 'testAgnCat.sav'
        agn = testGalaxyAgnDBObj(address=self.connectionString)
        cat = testAgnCatalog(agn, obs_metadata = self.obs_metadata)
        cat.write_catalog(catName)
        self.catalogTester(catName=catName, catalog=cat, nameRoot='agn')
        if os.path.exists(catName):
            os.unlink(catName)


    def testPSFimages(self):
        """
        Test that GalSimInterpreter puts the right number of counts on images of Galaxy bulges convolved
        with a PSF
        """
        catName = 'testPSFcat.sav'
        gals = testGalaxyBulgeDBObj(address=self.connectionString)
        cat = psfCatalog(gals, obs_metadata = self.obs_metadata)
        cat.write_catalog(catName)
        self.catalogTester(catName=catName, catalog=cat, nameRoot='psf')
        if os.path.exists(catName):
            os.unlink(catName)

    def testMultipleImages(self):
        """
        Test that GalSimInterpreter puts the right number of counts on images of multiple objects
        """
        dbName = 'galSimTestMultipleDB.db'
        if os.path.exists(dbName):
            os.unlink(dbName)

        displacedRA = numpy.array([72.0/3600.0, 55.0/3600.0, 75.0/3600.0])
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
        """
        Test that GalSimInterpreter puts the right number of counts on images containgin different types of objects
        """
        dbName1 = 'galSimTestCompound1DB.db'
        if os.path.exists(dbName1):
            os.unlink(dbName1)

        displacedRA = numpy.array([72.0/3600.0, 55.0/3600.0, 75.0/3600.0])
        displacedDec = numpy.array([0.0, 15.0/3600.0, -15.0/3600.0])
        obs_metadata1 = makePhoSimTestDB(filename=dbName1, size=1,
                                            displacedRA=displacedRA, displacedDec=displacedDec)
        connectionString1 = 'sqlite:///'+dbName1

        dbName2 = 'galSimTestCompound2DB.db'
        if os.path.exists(dbName2):
            os.unlink(dbName2)

        displacedRA = numpy.array([55.0/3600.0, 60.0/3600.0, 62.0/3600.0])
        displacedDec = numpy.array([-3.0/3600.0, 10.0/3600.0, 10.0/3600.0])
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
