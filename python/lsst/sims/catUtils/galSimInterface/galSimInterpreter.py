"""
This script contains all of the methods needed to turn a GalSimCatalog
(see exampleCatalogExamples/galSimCatalogExamples.py) into FITS images

There is the class GalSimInterpreter, which can be imported into other scripts
(assuming you are running a stack-enabled Python) and a main function which uses
it to create FITS images.
"""

import os
import numpy
import galsim
from lsst.sims.catalogs.generation.db import radiansToArcsec

__all__ = ["GalSimInterpreter", "GalSimDetector", "radiansToArcsec"]

class GalSimDetector(object):
    """
    This class stores information about individual detectors for use by the GalSimInterpreter
    """

    def __init__(self, name=None, xCenter=None, yCenter=None,
                 xMin=None, xMax=None, yMin=None, yMax=None,
                 plateScale=None):
        """
        param [in] name is a string denoting the name of the detector (this should be the
        same name that will be returned by the astrometry method findChipName()

        param [in] xCenter is the x pupil coordinate of the center of the detector in arcseconds

        param [in] yCenter is the y pupil coordinate of the cneter of the detector in arcseconds

        param [in] xMin, xMax, yMin, yMax are the corresponding minimum and maximum values of the
        pupil coordinates on this detector in arcseconds

        param [in] plateScale in arcseconds per pixel on this detector
        """

        self.name = name
        self.xCenter = xCenter
        self.yCenter = yCenter
        self.xMin = xMin
        self.xMax = xMax
        self.yMin = yMin
        self.yMax = yMax
        self.plateScale = plateScale

class GalSimInterpreter(object):
    """
    This class will read in a GalSimCatalog from a file and write
    FITS images for each of the detectors specified in the catalog.

    The basic work flow is:

    myInterpreter = GalSimInterpreter()

    myInterpreter.readCatalog('myCatalog.dat')

    myInterpreter.drawCatalog(fileNameRoot='myFitsImages', bandpass='/path/to/bandpass/file.dat')

    This will produce a series of FITS files named myFitsImages_R_0_0_S_0_0.fits, etc., one for each
    detector in the catalog.
    """

    def __init__(self, detectors=None):

        #in case we want to draw images using the Fourier transform
        self.bigfft = galsim.GSParams(maximum_fft_size=10000)

        self.data = None
        
        if detectors is None:
            raise RuntimeError("Will not create images; you passed no detectors to the GalSimInterpreter")
        
        self.detectors = detectors

        #this is a list of which detectors which objects fall on
        #(in case an object is near the edge of a detector and some
        #of the light from the object falls on that detector)
        self.chipsImpinged = None

        self.detectorImages = {}
        self.bandPasses = {}


    def _doesObjectImpingeOnChip(self, xPupil=None, yPupil=None, halfLightRadius=None,
                                 minorAxis=None, majorAxis=None, detector=None):
        """
        Compare an object to a detector and determine whether or not that object will cast any
        light on that detector (in case the object is near the edge of a detector and will cast some
        incidental light onto it).

        param [in] obj the row from self.data describing the astronomical object

        param [in] detector a GalSimDetector object describing the detector

        returns True if the object does illuminate the detector; False if not
        """

        #12 November 2014
        #The method below is not terribly clever.  It triples the half light radius of the object,
        #maps it to an ellipse with the same major-to-minor axis ratio as the actual object,
        #then assumes that the object is a circle with this (enlarged) half light radius.  Any detector
        #that is within this radius of the object is considered illuminated by the object.
        #This was meant to be a very safe estimate (choosing detectors that are farther from the object
        #than they probably need to be).  Something more clever and targeted is welcome

        if xPupil >= detector.xMin and \
           xPupil <= detector.xMax:

            isBetweenX = True
        else:
            isBetweenX = False

        if yPupil >= detector.yMin and \
           yPupil <= detector.yMax:

            isBetweenY = True
        else:
            isBetweenY = False

        if isBetweenX and isBetweenY:
            #the object falls on the detector directly
            return True

        radius = 3.0*halfLightRadius
        ratio = minorAxis/majorAxis
        distance = radius/numpy.sqrt(ratio)

        #check if light from the object bleed across any of the detector's boundaries
        if isBetweenY:
            if xPupil <= detector.xMin and detector.xMin - xPupil < distance:
                return True

            if xPupil >= detector.xMax and xPupil - detector.xMax < distance:
                return True

        if isBetweenX:
            if yPupil <= detector.yMin and detector.yMin - yPupil < distance:
                return True

            if yPupil >= detector.yMax and yPupil - detector.yMax < distance:
                return True

        #see if light from the object bleeds through any of the detector's corners
        for xx in [detector.xMin, detector.xMax]:
            for yy in [detector.yMin, detector.yMax]:
                testDistance = numpy.sqrt(numpy.power(xx - xPupil,2) + \
                           numpy.power(yy - yPupil,2))

                if testDistance < distance:
                    return True

        return False


    def findAllChips(self, xPupil=None, yPupil=None, halfLightRadius=None, minorAxis=None, majorAxis=None):
        outputString = ''
        outputList = []
        for dd in self.detectors:
            if self._doesObjectImpingeOnChip(xPupil=radiansToArcsec(xPupil), yPupil=radiansToArcsec(yPupil),
                                             halfLightRadius=radiansToArcsec(halfLightRadius),
                                             minorAxis=radiansToArcsec(minorAxis), majorAxis=radiansToArcsec(majorAxis),
                                             detector=dd):
            
                if outputString != '':
                    outputString += '//'
                outputString += dd.name
                outputList.append(dd)
        return outputString, outputList

    def drawCatalog(self, fileNameRoot=None, bandPass=None):
        """
        Loop through the detectors read in by readCatalog and create FITS images for each

        param [in] fileNameRoot is a string containing the root of the names of the resulting FITS
        files.  The output files will be named fileNameRoot_detectorName.fits, i.e.
        fileNameRoot_R_0_0_S_1_2.fits, etc.

        param [in] bandPass is a string containing the name of the file which contains the filter
        bandpass over which you want to integrate the flux.  This file should be two columns: wavelength
        in nanometers and throughput from 0.0 to 1.0
        """

        #read in the bandpass
        bp = galsim.Bandpass(bandPass)
        for detector in self.detectors:
            self.drawImage(fileNameRoot=fileNameRoot, bandPass=bp, detector=detector)

    def _getImageName(self, fileNameRoot='', detector=None):
        #format the name of the detector to add to the name of the FITS file
        detectorName = detector.name
        detectorName = detectorName.replace(',','_')
        detectorName = detectorName.replace(':','_')
        detectorName = detectorName.replace(' ','_')

        fileName = fileNameRoot+detectorName+'.fits'
        return fileName

    def initializeImage(self, fileName='', detector=None):
        """
        Draw the FITS file associated with a specific detector

        param [in] fileNameRoot is a string containing the root of the names of the resulting fits
        files.  The output files will be named fileNameRoot_detectorName.fits, i.e.
        fileNameRoot_R_0_0_S_1_2.fits, etc.

        param [in] bandPass is galsim.bandpass object denoting the bandpass over which to integrate flux

        param [in] detector is a GalSimDetector object indicating the detector for which to draw the FITS file
        """

        #set the size of the image
        nx = int((detector.xMax - detector.xMin)/detector.plateScale)
        ny = int((detector.yMax - detector.yMin)/detector.plateScale)
        image = galsim.Image(nx,ny)
        
        self.detectorImages[fileName] = image

    def drawObject(self, galSimType=None, detectorList=None, fileNameRoot='', bandPassName=None, **kwargs):
        if bandPassName is None:
            return
        
        if bandPassName in self.bandPasses:
            bandPass = self.bandPasses[bandPassName]
        else:
            bandPass = galsim.Bandpass(bandPassName)
            self.bandPasses[bandPassName] = bandPass
        
        #draw all of the objects who illumine this detector
        fileNameList = []
        for dd in detectorList:
            fileName = self._getImageName(fileNameRoot=fileNameRoot, detector=dd)
            fileNameList.append(fileName)
            if fileName not in self.detectorImages:    
                self.initializeImage(fileName=fileName, detector=dd)
        
        if galSimType == 'galaxy':
            self.drawGalaxy(fileNameList=fileNameList, detectorList=detectorList, 
                                bandPass=bandPass, **kwargs)
        else:
            print "Apologies: the GalSimInterpreter does not yet have a method to draw "
            print objectParams['galSimType']
            print " objects\n"

    def drawGalaxy(self, fileNameList=None, detectorList=None, bandPass=None, sindex=None, minorAxis=None,
                   majorAxis=None, positionAngle=None, halfLightRadius=None, 
                   x_pupil=None, y_pupil=None, sed=None):
        """
        Draw the image of a galaxy.

        param [in] entry is an element from self.data containing the information on an astronomical object

        param [in] image is a galsim Image object into which we will draw this object

        param [in] detector is a GalSimDetector object denoting the detectror with which the image is associated

        param [in] bandPass is a galsim.bandpass object denoting the bandpass over which to integrate flux
        """
        
        if sed is None:
            return
        
        hlr = radiansToArcsec(halfLightRadius)
        minor = radiansToArcsec(minorAxis)
        major = radiansToArcsec(majorAxis)
        xp = radiansToArcsec(x_pupil)
        yp = radiansToArcsec(y_pupil)

        #create a Sersic profile
        centeredObj = galsim.Sersic(n=float(sindex), half_light_radius=float(hlr))

        #turn the Sersic profile into an ellipse
        centeredObj = centeredObj.shear(q=minor/major, beta=positionAngle*galsim.radians)
        
        for (detector, fileName) in zip(detectorList, fileNameList):
            
            image = self.detectorImages[fileName]
            
            #by default, galsim draws objects at the center of the image;
            #this will shift the object into its correct position
            dx=xp-detector.xCenter
            dy=yp-detector.yCenter
            obj = centeredObj.shift(dx, dy)

            #declare a spectrum() function for galSim to use
            #spectrum = galsim.SED(objectParams['galSimSedName'],
            #                         flux_type='flambda')

            spectrum = galsim.SED(spec = lambda ll: numpy.interp(ll, sed.wavelen, sed.flambda),
                                  flux_type='flambda')

            #convolve the Sersic profile with the spectrum
            obj = obj*spectrum

            #add this object to the image, integrating over the bandPass
            #Note: by specifying method='real_space', we are asking GalSim
            #to directly integrate over the pixels.  Other options are to use
            #method='fft' and integrate using a Fourier transform (though this method can
            #be finnicky for large images) or method='phot' in which case photons are drawn
            #from the SED and shot at the chip (a la phoSim)
            image = obj.drawImage(bandpass=bandPass, scale=detector.plateScale, image=image,
                                  add_to_image=True, method='real_space')

    def writeImages(self):
        for imageName in self.detectorImages:
            self.detectorImages[imageName].write(file_name=imageName)
            


def main():
    """
    This method reads in the galSim_example.txt catalog created by
    galSimCatalogGenerator.py and uses it to draw FITS files for each
    of the detectors defined in that catalog.

    See the documentation at the top of galSimCatalogGenerator.py for an
    explanation of why the two scripts must be separate.

    Run this script using the python for which you have GalSim installed.

    Be sure to run it in a shell in which all of the LSST stack environment
    variables have been set, otherwise you will not have access to all of the
    photUtils functionality needed by the GalSimInterpreter.
    """

    try:
        #specify a bandpass through which to observe the galaxies
        bandPass = os.path.join(os.getenv('THROUGHPUTS_DIR'),'baseline','total_g.dat')
    except AttributeError:
        print "You must set the environment variable THROUGPUTS_DIR to point"
        print "to wherever the LSST throughputs are stored on your machine\n"
        print "probably something like $LSST_HOME/yourOS/throughputs/version/"
        exit()

    gs = GalSimInterpreter()

    #read in our catalog of galaxy bulges
    gs.readCatalog('galSim_example.txt')

    #write the images to files of the name galsimTest_detectorName.fits
    name = 'galsimTest_'
    gs.drawCatalog(bandPass=bandPass, fileNameRoot=name)

if __name__ == "__main__":
    main()
