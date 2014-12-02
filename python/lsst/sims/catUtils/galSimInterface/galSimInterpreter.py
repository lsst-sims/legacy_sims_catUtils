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
                 plateScale=None, fileNameRoot=''):
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
        self.fileName = self._getFileName(fileNameRoot=fileNameRoot)


    def _getFileName(self, fileNameRoot=''):
        #format the name of the detector to add to the name of the FITS file
        detectorName = self.name
        detectorName = detectorName.replace(',','_')
        detectorName = detectorName.replace(':','_')
        detectorName = detectorName.replace(' ','_')

        name = fileNameRoot+detectorName
        return name

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

    def __init__(self, detectors=None, bandPassNames=None, bandPassFiles=None):

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
        
        self.setBandPasses(bandPassNames=bandPassNames, bandPassFiles=bandPassFiles)

    def setBandPasses(self, bandPassNames=None, bandPassFiles=None):
        for bpn, bpf in zip(bandPassNames, bandPassFiles):
            bp = galsim.Bandpass(bpf)
            self.bandPasses[bpn] = bp

    def _getFileName(self, detector=None, bandPassName=None):
        return detector.fileName+'_'+bandPassName+'.fits'

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

    def initializeImage(self, detector=None):
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
        
        return image

    def drawObject(self, galSimType=None, detectorList=None, fileNameRoot='', sed=None, x_pupil=None,
                   y_pupil=None, **kwargs):
        
        if sed is None:
            return
        
        for dd in detectorList:
            for bandPassName in self.bandPasses:
                name = self._getFileName(detector=dd, bandPassName=bandPassName)
                if name not in self.detectorImages:    
                    self.detectorImages[name] = self.initializeImage(detector=dd)
        
        xp = radiansToArcsec(x_pupil)
        yp = radiansToArcsec(y_pupil)
        
        if galSimType == 'galaxy':
            centeredObj = self.drawGalaxy(**kwargs)
        else:
            print "Apologies: the GalSimInterpreter does not yet have a method to draw "
            print objectParams['galSimType']
            print " objects\n"
            return
            
        for detector in detectorList:
            
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
            
            for bandPassName in self.bandPasses:
                image = self.detectorImages[self._getFileName(detector=detector, bandPassName=bandPassName)]
                bandPass = self.bandPasses[bandPassName]

                #add this object to the image, integrating over the bandPass
                #Note: by specifying method='real_space', we are asking GalSim
                #to directly integrate over the pixels.  Other options are to use
                #method='fft' and integrate using a Fourier transform (though this method can
                #be finnicky for large images) or method='phot' in which case photons are drawn
                #from the SED and shot at the chip (a la phoSim)
                image = obj.drawImage(bandpass=bandPass, scale=detector.plateScale, image=image,
                                      add_to_image=True, method='real_space')

    def drawGalaxy(self, sindex=None, minorAxis=None,
                   majorAxis=None, positionAngle=None, halfLightRadius=None):
        """
        Draw the image of a galaxy.

        param [in] entry is an element from self.data containing the information on an astronomical object

        param [in] image is a galsim Image object into which we will draw this object

        param [in] detector is a GalSimDetector object denoting the detectror with which the image is associated

        param [in] bandPass is a galsim.bandpass object denoting the bandpass over which to integrate flux
        """
        
        hlr = radiansToArcsec(halfLightRadius)
        minor = radiansToArcsec(minorAxis)
        major = radiansToArcsec(majorAxis)

        #create a Sersic profile
        centeredObj = galsim.Sersic(n=float(sindex), half_light_radius=float(hlr))

        #turn the Sersic profile into an ellipse
        centeredObj = centeredObj.shear(q=minor/major, beta=positionAngle*galsim.radians)
        
        return centeredObj


    def writeImages(self):
        for imageName in self.detectorImages:
            self.detectorImages[imageName].write(file_name=imageName)
