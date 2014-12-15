"""
This file defines the following classes:

GalSimInterpreter -- a class which takes objects passed by a GalSim Instance Catalog
(see galSimCatalogs.py) and uses GalSim to write them to FITS images.

GalSimDetector -- a class which stored information about a detector in a way that
GalSimInterpreter expects.

PSFbase and various examples of PSFs -- classes which wrap the PSF implementations
from GalSim in a way that GalSimInterpreter expects.
"""

import os
import numpy
import galsim
from lsst.sims.catalogs.generation.db import radiansToArcsec

__all__ = ["GalSimInterpreter", "GalSimDetector", "PSFbase",
           "ExampleGaussianPSF", "ExampleOpticalPSF"]

class GalSimDetector(object):
    """
    This class stores information about individual detectors for use by the GalSimInterpreter
    """

    def __init__(self, name=None, xCenter=None, yCenter=None,
                 xMin=None, xMax=None, yMin=None, yMax=None,
                 plateScale=None):
        """
        param [in] name is a string denoting the name of the detector (this should be the
        same name that will be returned by the astrometry method findChipName())

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
        self.fileName = self._getFileName()


    def _getFileName(self):
        #format the name of the detector to add to the name of the FITS file
        detectorName = self.name
        detectorName = detectorName.replace(',','_')
        detectorName = detectorName.replace(':','_')
        detectorName = detectorName.replace(' ','_')

        name = detectorName
        return name

class PSFbase(object):
        """
        This is the base class for wrappers of GalSim's PSF classes.  To apply a PSF to GalSim images
        using the GalSim Instance Catalog and GalSim Interpreter, the user must define a daughter
        class of this class and instantiate it as the member variable self.PSF in the GalSim Instance Catalog.
        
        Any Daughter class of this class must have:
        
        1) a boolean member variable wavelength_dependent which tells the GalSimInterpreter whether or not
        it needs to worry about the PSF changing with wavelength
        
        2) a member method _getPSF which accepts the coordinates x_pupil and y_pupil in arcseconds as kwargs
        and (optionally) the kwarg bandpass, which is a galsim bandpass object.  This method will instantiate
        a psf object at those coordinates (and, if relevant, the effective wavelength) of the bandpass, and return
        it.
        
        The method applyPSF is defined in this class and should not be overwritten.  It handles the task of actually
        convolving the PSF returned by _getPSF.
        """
        
        wavelength_dependent = False
        
        def _getPSF(self, x_pupil=None, y_pupil=None, bandpass=None):
            """
            If it had been implemented, this would return a GalSim PSF instantiation at the
            coordinates and wavelength specified and returned it to applyPSF.  As it is, this
            class has not been implemented and is left to the user to implement in Daughter
            classes of PSFbase.
            
            @param [in] x_pupil the x coordinate on the pupil in arc seconds
            
            @param [in] y_pupil the y coordinate on the pupil in arc seconds
            
            @param [in] bandpass is an instantiation of the GalSim bandpass class which contains
            data defining the bandpass in question (in case the PSF is wavelength dependent)
            """
            
            raise NotImplementedError("There is not _getPSF for PSFbase; define a daughter class and define your own")
        
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

class ExampleGaussianPSF(PSFbase):

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

class ExampleOpticalPSF(PSFbase):
    
    wavelength_dependent = True
    
    def _getPSF(self, x_pupil=None, y_pupil=None, **kwargs):
        eff = kwargs['bandpass'].effective_wavelength
        psf = galsim.OpticalPSF(lam_over_diam=radiansToArcsec(eff*1.0e-9/8.0), astig1=1.0,
                                astig2=2.0)
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

        self.detectorImages = {}
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
                    if xPupil>100.0:
                        print 'returning because of xMin'
                    return True

                if xPupil >= detector.xMax and xPupil - detector.xMax < distance:
                    if xPupil>100.0:
                        print 'returning because of xMax'
                    return True

            if isBetweenX:
                if yPupil <= detector.yMin and detector.yMin - yPupil < distance:
                    if xPupil>100.0:
                        print 'returning because of yMin'
                    return True

                if yPupil >= detector.yMax and yPupil - detector.yMax < distance:
                    if xPupil>100.0:
                        print 'returning because of yMax'
                    return True

            #see if light from the object bleeds through any of the detector's corners
            for xx in [detector.xMin, detector.xMax]:
                for yy in [detector.yMin, detector.yMax]:
                    testDistance = numpy.sqrt(numpy.power(xx - xPupil,2) + \
                               numpy.power(yy - yPupil,2))

                    if testDistance < distance:
                        if xPupil>100.0:
                            print 'returning because corner ',testDistance,distance
                            print xx,yy
                            print xPupil, yPupil
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
        image = galsim.Image(nx, ny, scale=detector.plateScale)
        
        return image

    def drawObject(self, galSimType=None, detectorList=None, fileNameRoot='', sed=None, x_pupil=None,
                   y_pupil=None, **kwargs):
        
        if sed is None or len(detectorList) == 0:
            return
        
        for detector in detectorList:
            for bandPassName in self.bandPasses:
                name = self._getFileName(detector=detector, bandPassName=bandPassName)
                if name not in self.detectorImages:
                    self.detectorImages[name] = self.blankImage(detector=detector)
        
        centeredObj = None
        xp = radiansToArcsec(x_pupil)
        yp = radiansToArcsec(y_pupil)
        spectrum = galsim.SED(spec = lambda ll: numpy.interp(ll, sed.wavelen, sed.flambda),
                              flux_type='flambda')
        
        for bandPassName in self.bandPasses:
            if centeredObj is None or (self.PSF is not None and self.PSF.wavelength_dependent):
                if galSimType == 'sersic':
                    centeredObj = self.drawSersic(x_pupil=xp, y_pupil=yp,
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
                
                dx = xp - detector.xCenter
                dy = yp - detector.yCenter
                obj = centeredObj.shift(dx, dy)
                
                #convolve the Sersic profile with the spectrum
                obj = obj*spectrum
                localImage = self.blankImage(detector=detector)
                localImage = obj.drawImage(bandpass=self.bandPasses[bandPassName], scale=detector.plateScale,
                                           method='phot', gain=self.gain, image=localImage)
                
                self.detectorImages[name] += localImage

    def drawPointSource(self, x_pupil=None, y_pupil=None, bandpass=None):
    
        if self.PSF is None:
            raise RuntimeError("Cannot draw a point source in GalSim without a PSF")
        
        return self.PSF.applyPSF(x_pupil=x_pupil, y_pupil=y_pupil, bandpass=bandpass)

    def drawSersic(self, x_pupil=None, y_pupil=None, sindex=None, minorAxis=None,
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

    def writeImages(self, nameRoot=None):
        for name in self.detectorImages:
            if nameRoot is not None:
                fileName = nameRoot+'_'+name
            else:
                fileName = name
            self.detectorImages[name].write(file_name=fileName)
