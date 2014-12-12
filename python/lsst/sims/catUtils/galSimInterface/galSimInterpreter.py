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

__all__ = ["GalSimInterpreter", "GalSimDetector", "PSFbase", "ExampleGalSimPSF"]

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

class PSFbase(object):
        
        wavelength_dependent = False
    
        def applyPSF(self, x_pupil=None, y_pupil=None, obj=None, **kwargs):
            """
            Apply the PSF to a GalSim GSObject
        
            In theory, this method will accept the x and y pupil coordinates in arc seconds as well
            as a GalSim GSObject.  The method will calculate the PSF parameters based on x_pupil
            and y_pupil, construct a Galsim GSObject corresponding to the PSF function, and convolve
            the PSF with the baseObject, returning the rsult of the convolution.
        
            Must have the option to return the raw psf in the case of point sources.
        
            This example uses a Gaussian PSF.
            """
            psf = self._getPSF(x_pupil=x_pupil, y_pupil=y_pupil, **kwargs)
            if obj is not None:
                obj = galsim.Convolve(obj, psf)
                return obj
            else:
                return psf

class ExampleGalSimPSF(PSFbase):

    def _getPSF(self, x_pupil=None, y_pupil=None, **kwargs):
        """
        Apply the PSF to a GalSim GSObject
        
        In theory, this method will accept the x and y pupil coordinates in arc seconds as well
        as a GalSim GSObject.  The method will calculate the PSF parameters based on x_pupil
        and y_pupil, construct a Galsim GSObject corresponding to the PSF function, and convolve
        the PSF with the baseObject, returning the rsult of the convolution.
        
        Must have the option to return the raw psf in the case of point sources.
        
        This example uses a Gaussian PSF.
        """
        psf = galsim.Gaussian(sigma=0.14)
        psf = psf.shear(q=0.05, beta=numpy.pi*0.25*galsim.radians)
        return psf

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

    def __init__(self, detectors=None, bandPassNames=None, bandPassFiles=None,
                 gain=2.3):

        self.PSF = None
        self.gain = gain

        if detectors is None:
            raise RuntimeError("Will not create images; you passed no detectors to the GalSimInterpreter")
        
        self.detectors = detectors

        #this is a list of which detectors which objects fall on
        #(in case an object is near the edge of a detector and some
        #of the light from the object falls on that detector)
        self.chipsImpinged = None

        self.detectorObjects = {}
        self.bandPasses = {}
        
        self.setBandPasses(bandPassNames=bandPassNames, bandPassFiles=bandPassFiles)

    def setBandPasses(self, bandPassNames=None, bandPassFiles=None):
        for bpn, bpf in zip(bandPassNames, bandPassFiles):
            bp = galsim.Bandpass(bpf)
            self.bandPasses[bpn] = bp

    def setPSF(self, PSF=None):
        self.PSF=PSF

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

        if halfLightRadius==0.0 or minorAxis==0.0 or majorAxis==0.0:
            #I am not sure in the case of point sources how to deal with this,
            #since there is not general PSF formalism with a defined size.
            #For the moment, I will do a very conservative test (letting in more objects
            #than is probably necessary for each chip).  I will allow anything that is
            #within 2 arcseconds of the detector's boundaries (on the assumption that no
            #reasonably PSF would smear a source more than 2 arcseconds)
            
            if xPupil < detector.xMin - 2.0 or xPupil > detector.xMax + 2.0:
                return False
            
            if yPupil < detector.yMin -2.0 or yPupil > detector.yMax + 2.0:
                return False
                
            return True
            
        else:
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
        
        if outputString == '':
            outputString = None
        
        return outputString, outputList

    def blankImage(self, detector=None):
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
        
        if sed is None or len(detectorList) == 0:
            return
        
        for dd in detectorList:
            for bandPassName in self.bandPasses:
                name = self._getFileName(detector=dd, bandPassName=bandPassName)
                if name not in self.detectorObjects:
                    self.detectorObjects[name] = []
        
        centeredObj = None
        xp = radiansToArcsec(x_pupil)
        yp = radiansToArcsec(y_pupil)
        spectrum = galsim.SED(spec = lambda ll: numpy.interp(ll, sed.wavelen, sed.flambda),
                              flux_type='flambda')
        
        for bandPassName in self.bandPasses:
            if centeredObj is None or (self.PSF is not None and self.PSF.wavelength_dependent):
                if galSimType == 'sersic':
                    centeredObj = self.drawGalaxy(x_pupil=xp, y_pupil=yp,
                                                  bandpass=self.bandPasses[bandPassName], **kwargs)
                elif galSimType == 'pointSource':
                    centeredObj = self.drawPointSource(x_pupil=xp, y_pupil=yp,
                                                       bandpass=self.bandPasses[bandPassName])
                else:
                    print "Apologies: the GalSimInterpreter does not yet have a method to draw "
                    print objectParams['galSimType']
                    print " objects\n"
                    return
            
            for detector in detectorList:
            
                name = self._getFileName(detector=detector, bandPassName=bandPassName)
            
                #by default, galsim draws objects at the center of the image;
                #this will shift the object into its correct position
                dx=xp-detector.xCenter
                dy=yp-detector.yCenter
                
                #12 December 2014
                #tests indicate that this will not overwrite centeredObj but will, in fact,
                #create a new object (which is what we want)
                obj = centeredObj.shift(dx, dy)
            
                #convolve the Sersic profile with the spectrum
                obj = obj*spectrum
                self.detectorObjects[name].append(obj)


    def drawPointSource(self, x_pupil=None, y_pupil=None, bandpass=None):
    
        if self.PSF is None:
            raise RuntimeError("Cannot draw a point source in GalSim without a PSF")
        
        return self.PSF.applyPSF(x_pupil=x_pupil, y_pupil=y_pupil, bandpass=bandpass)

    def drawGalaxy(self, x_pupil=None, y_pupil=None, sindex=None, minorAxis=None,
                   majorAxis=None, positionAngle=None, halfLightRadius=None, bandpass=None):
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
        if self.PSF is not None:
            centeredObj = self.PSF.applyPSF(x_pupil=x_pupil, y_pupil=y_pupil, obj=centeredObj,
                                            bandpass=bandpass)
            
        return centeredObj

    def compressObjectList(self):
        for name in self.detectorObjects:
            obj = galsim.ChromaticSum(self.detectorObjects[name])
            self.detectorObjects[name] = [obj]

    def writeImages(self, isTest=False, nameRoot=None):
        for detector in self.detectors:
            for bandPassName in self.bandPasses:
                name = self._getFileName(detector=detector, bandPassName=bandPassName)
                if name in self.detectorObjects:
                    obj = galsim.ChromaticSum(self.detectorObjects[name])
                    image = self.blankImage(detector=detector)
                    if isTest:
                        image = obj.drawImage(bandpass=self.bandPasses[bandPassName], scale=detector.plateScale,
                                              image=image, add_to_image=True, method='real_space', gain=self.gain)
                    else:
                        image = obj.drawImage(bandpass=self.bandPasses[bandPassName], scale=detector.plateScale,
                                              image=image, add_to_image=True, method='phot', gain=self.gain)
                    
                    if nameRoot is not None:
                        name = nameRoot + '_' + name
                        
                    image.write(file_name=name)
