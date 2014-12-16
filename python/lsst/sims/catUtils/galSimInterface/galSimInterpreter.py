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
           "ExampleGaussianPSF", "ExampleOpticalPSF", "ExampleCCDNoise"]

class GalSimDetector(object):
    """
    This class stores information about individual detectors for use by the GalSimInterpreter
    """

    def __init__(self, name=None, xCenter=None, yCenter=None,
                 xMin=None, xMax=None, yMin=None, yMax=None,
                 plateScale=None, electronsPerADU=1.71234e3, readNoise=52.1237):
        """
        param [in] name is a string denoting the name of the detector (this should be the
        same name that will be returned by the astrometry method findChipName())

        param [in] xCenter is the x pupil coordinate of the center of the detector in arcseconds

        param [in] yCenter is the y pupil coordinate of the cneter of the detector in arcseconds

        param [in] xMin, xMax, yMin, yMax are the corresponding minimum and maximum values of the
        pupil coordinates on this detector in arcseconds

        param [in] plateScale in arcseconds per pixel on this detector
        
        param [in] electronsPerADU default value is taken from afw/testUtils.py line 76
        
        param [in] readNoise per pixel in electrons
        default value is taken from afw/testUtils.py line 77

        This class will generate its own internal variable self.fileName which is
        the name of the detector as it will appear in the output FITS files
        """

        self.name = name
        self.xCenter = xCenter
        self.yCenter = yCenter
        self.xMin = xMin
        self.xMax = xMax
        self.yMin = yMin
        self.yMax = yMax
        self.plateScale = plateScale
        self.electronsPerADU = electronsPerADU
        self.readNoise = readNoise
        self.fileName = self._getFileName()


    def _getFileName(self):
        """
        Format the name of the detector to add to the name of the FITS file
        """
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

    Consult GalSim's documentation to see what kinds of PSFs are available.

    See the classes ExampleGaussianPSF and ExampleOpticalPSF below for example implementations.

    See galSimCompoundGenerator.py and galSimStarGenerator.py for example usages.
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

        This method accepts the x and y pupil coordinates in arc seconds as well
        as a GalSim GSObject.  The method calculates the PSF parameters based on x_pupil
        and y_pupil, constructs a Galsim GSObject corresponding to the PSF function, and convolves
        the PSF with the GSObject, returning the result of the convolution.

        In the case of point sources, this object returns the raw PSF, rather than attempting
        a convolution (since there is nothing to convolve with).

        @param [in] x_pupil the x pupil coordinate in arc seconds

        @param [in] y_pupil the y pupil coordinate in arc seconds

        @param [in] obj is a GalSim GSObject (an astronomical object) with which
        to convolve the PSF (optional)

        **kwargs is there so that a bandpass can also be passed in and sent to _getPSF
        """

        #use the user-defined _getPSF method to calculate the PSF at these specific
        #coordinates and (optionally) wavelength
        psf = self._getPSF(x_pupil=x_pupil, y_pupil=y_pupil, **kwargs)

        if obj is not None:
            #if we are dealing with an extended object, convolve it with the psf
            obj = galsim.Convolve(obj, psf)
            return obj
        else:
            #if there is no object (i.e. if this is a point source), just return the PSF
            return psf

class ExampleGaussianPSF(PSFbase):
    """
    This is an example implementation of a wavelength- and position-independent
    Gaussian PSF.  See the documentation in PSFbase to learn how it is used.
    """

    wavelength_dependent = False

    def _getPSF(self, x_pupil=None, y_pupil=None, **kwargs):
        """
        Return a Gaussian PSF to be convolved with sources.

        @param [in] x_pupil the x coordinate on the pupil in arc seconds

        @param [in] y_pupil the y coordinate on the pupil in arc seconds

        @param [in] bandpass is an instantiation of the GalSim bandpass class which contains
        data defining the bandpass in question (in case the PSF is wavelength dependent)
        """
        psf = galsim.Gaussian(sigma=0.14)
        psf = psf.shear(q=0.05, beta=numpy.pi*0.25*galsim.radians)
        return psf

class ExampleOpticalPSF(PSFbase):
    """
    This is an example implementation of a position-independent version of GalSims OpticalPSF class.
    See documentation for PSFbase to learn how it is used.
    """

    wavelength_dependent = True

    def _getPSF(self, x_pupil=None, y_pupil=None, **kwargs):
        """
        Return an OpticalPSF to be convolved with sources.

        @param [in] x_pupil the x coordinate on the pupil in arc seconds

        @param [in] y_pupil the y coordinate on the pupil in arc seconds

        @param [in] bandpass is an instantiation of the GalSim bandpass class which contains
        data defining the bandpass in question (in case the PSF is wavelength dependent)
        """

        eff = kwargs['bandpass'].effective_wavelength
        psf = galsim.OpticalPSF(lam_over_diam=radiansToArcsec(eff*1.0e-9/8.0), astig1=1.0,
                                astig2=2.0)
        return psf

class ExampleCCDNoise(object):
    """
    This class wraps the GalSim class CCDNoise.  It is meant to be assigned as
    the self.noise member variable in a GalSim InstanceCatalog.  To instantiate
    a different noise model, write a class like this one that defines a method
    getNoiseModel() which accepts as its argument an ObservationMetaData
    instantiation (see 
    
    sims_catalogs_generation/python/lsst/sims/catalogs/generation/db/ObservationMetaData.py
    
    for definition) and a GalSimDetector instantiation and returns an instantiation of 
    a GalSim noise model
    """
    
    def __init__(self, seed=None):
        if seed is None:
            self.randomNumbers = galsim.UniformDeviate()
        else:
            self.randomNumbers = galsim.UniformDeviate(seed)
    
    def getNoiseModel(self, obs_metadata=None, detector=None):
        skyLevel = 10.0 #this is obviously nonsense; GalSim wants electrons-per-pixel; Peter thinks we are storing
                        #sky brightness in magnitudes per square arc-second; once we have the new interface to
                        #OpSim written, we will need to use that, plus filter information, to convert 
                        #between the two (which, in principle, can be done)

        gain = detector.electronsPerADU
        readNoise = detector.readNoise
        
        return galsim.CCDNoise(self.randomNumbers, sky_level=skyLevel, gain=gain, read_noise=readNoise)


class GalSimInterpreter(object):
    """
    This is the class which actually takes the objects contained in the GalSim Instance Catalog and converts them
    into FITS images.
    """

    def __init__(self, detectors=None, bandPassNames=None, bandPassFiles=None,
                 gain=2.3):

        """
        @param [in] detectors is a list of GalSimDetectors for which we are drawing FITS images

        @param [in] bandPassNames is a list of the form ['u', 'g', 'r', 'i', 'z', 'y'] denoting
        the bandpasses for which we are drawing FITS images.

        @param [in] bandPassFiles is a list of paths to the bandpass data files corresponding to
        bandPassNames

        @param gain is the number of photons per ADU for our detectors.  Note: this requires all
        detectors to have the same gain.  Maybe that is inappropriate...
        """

        self.PSF = None
        self.gain = gain

        if detectors is None:
            raise RuntimeError("Will not create images; you passed no detectors to the GalSimInterpreter")

        self.detectors = detectors

        self.detectorImages = {} #this dict will contain the FITS images (as GalSim images)
        self.bandPasses = {} #this dict will contain the GalSim bandpass instantiations corresponding to the input bandpasses

        self.setBandPasses(bandPassNames=bandPassNames, bandPassFiles=bandPassFiles)

    def setBandPasses(self, bandPassNames=None, bandPassFiles=None):
        """
        Read in files containing bandpass data and store them in a dict of GalSim bandpass instantiations.

        @param [in] bandPassNames is a list of the names by which the bandpasses are to be referred
        i.e. ['u', 'g', 'r', 'i', 'z', 'y']

        @param [in] bandPassFiles is a list of paths to the files containing data for the bandpasses

        The bandpasses will be stored in the member variable self.bandPasses, which is a dict
        """
        for bpn, bpf in zip(bandPassNames, bandPassFiles):
            bp = galsim.Bandpass(bpf)
            self.bandPasses[bpn] = bp

    def setPSF(self, PSF=None):
        """
        Set the PSF wrapper for this GalSimInterpreter

        @param [in] PSF is an instantiation of a class which inherits from PSFbase and defines _getPSF()
        """
        self.PSF=PSF

    def _getFileName(self, detector=None, bandPassName=None):
        """
        Given a detector and a bandpass name, return the name of the FITS file to be written

        @param [in] detector is an instantiation of GalSimDetector

        @param [in] bandPassName is a string i.e. 'u' denoting the filter being drawn

        The resulting filename will be detectorName_bandPassName.fits
        """
        return detector.fileName+'_'+bandPassName+'.fits'

    def _doesObjectImpingeOnDetector(self, xPupil=None, yPupil=None, halfLightRadius=None,
                                 minorAxis=None, majorAxis=None, detector=None):
        """
        Compare an object to a detector and determine whether or not that object will cast any
        light on that detector (in case the object is near the edge of a detector and will cast some
        incidental light onto it).

        @param [in] xPupil the x pupil coordinate of the object in arc seconds

        @param [in] yPupil the y pupil coordinate of the object in arc seconds

        @param [in] detector an instantiation of GalSimDetector.  This is the detector against
        which we will compare the object.

        @param [in] halfLightRadius the half light radius of the object in arc seconds (optional)

        @param [in] minorAxis the semi-minor axis of the object in arc seconds (optional)

        @param [in] majorAxis the semi-major axis of the object in arc seconds (optional)

        returns True if the object does illuminate the detector; False if not.

        Because it is possible for an extended source to cast light on a detector even if
        its center does not fall on the detector, this method will return True if the object's
        center is within 3 x halfLightRadius x sqrt(majorAxis/minorAxis) of the detector's bounds.

        In the case of a point source (halfLightRadius, minorAxis, majorAxis == 0.0),
        this method will return True if the object is within 2 arc seconds of the detector's
        bounds (on the assumption that no PSF will smear the object outmore than 2 arc seconds)
        """

        if halfLightRadius==0.0 or minorAxis==0.0 or majorAxis==0.0:
            #I am not sure in the case of point sources how to deal with this,
            #since there is not general PSF formalism with a defined size.
            #For the moment, I will do a very conservative test (letting in more objects
            #than is probably necessary for each detector).  I will allow anything that is
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


    def findAllDetectors(self, xPupil=None, yPupil=None, halfLightRadius=None, minorAxis=None, majorAxis=None):
        """
        For a given object, find all of the detectors on which it casts light.

        @params [in] xPupil the x pupil coordinate of the object in radians

        @params [in] yPupil the y pupil coordinate of the object in radiasn

        @params [in] halfLightRadius the half light radius of the object in radians

        @params [in] minorAxis the semi-minor axis of the object in radians

        @params [in] majorAxis the semi-major axis of the object in radians

        returns a string listing the names of all the detectors on which the object casts light. The names
        are separated by a '//'
        """
        outputString = ''
        outputList = []
        for dd in self.detectors:
            if self._doesObjectImpingeOnDetector(xPupil=radiansToArcsec(xPupil), yPupil=radiansToArcsec(yPupil),
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
        Draw a blank image associated with a specific detector.  The image will have the correct size
        for the given detector.

        param [in] detector is an instantiation of GalSimDetector
        """

        #set the size of the image
        nx = int((detector.xMax - detector.xMin)/detector.plateScale)
        ny = int((detector.yMax - detector.yMin)/detector.plateScale)
        image = galsim.Image(nx, ny, scale=detector.plateScale)

        return image

    def drawObject(self, galSimType=None, detectorList=None, sed=None, x_pupil=None,
                   y_pupil=None, halfLightRadius=None, minorAxis=None, majorAxis=None,
                   positionAngle=None, sindex=None):
        """
        Draw an object on all of the relevant FITS files.

        @param [in] galSimType is a string, either 'pointSource' or 'sersic' denoting the shape of the object

        @param [in] detectorList is a list of GalSimDetectors on which to draw the object

        @param [in] sed is the SED of the object (an instantiation of the Sed class defined in
        sims_photUtils/../../Sed.py

        @param [in] x_pupil is the x pupil coordinate of the object in radians

        @param [in] y_pupil is the y pupil coordinate of the object in radians

        @param [in] halfLightRadius is the halfLightRadius of the object in radians

        @param [in] minorAxis is the semi-minor axis of the object in radians

        @param [in] majorAxis is the semi-major axis of the object in radians

        @param [in] positionAngle is the position angle of the object in radians

        @param [in] sindex is the sersic index of the object
        """

        if sed is None or len(detectorList) == 0:
            #there is nothing to draw
            return

        #go through the list of detector/bandpass combinations and initialize
        #all of the FITS files we will need (if they have not already been initialized)
        for detector in detectorList:
            for bandPassName in self.bandPasses:
                name = self._getFileName(detector=detector, bandPassName=bandPassName)
                if name not in self.detectorImages:
                    self.detectorImages[name] = self.blankImage(detector=detector)

        centeredObj = None
        xp = radiansToArcsec(x_pupil)
        yp = radiansToArcsec(y_pupil)
        hlr = radiansToArcsec(halfLightRadius)
        spectrum = galsim.SED(spec = lambda ll: numpy.interp(ll, sed.wavelen, sed.flambda),
                              flux_type='flambda')

        for bandPassName in self.bandPasses:

            #create a new object if one has not already been created or if the PSF is wavelength
            #dependent (in which case, each filter is going to need its own initialized object)
            if centeredObj is None or (self.PSF is not None and self.PSF.wavelength_dependent):
                if galSimType == 'sersic':
                    centeredObj = self.drawSersic(x_pupil=xp, y_pupil=yp,
                                                  bandpass=self.bandPasses[bandPassName],
                                                  sindex=sindex, halfLightRadius=hlr,
                                                  positionAngle=positionAngle,
                                                  minorAxis=minorAxis, majorAxis=majorAxis)
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

                #convolve the object's shape profile with the spectrum
                obj = obj*spectrum
                localImage = self.blankImage(detector=detector)
                localImage = obj.drawImage(bandpass=self.bandPasses[bandPassName], scale=detector.plateScale,
                                           method='phot', gain=self.gain, image=localImage)

                self.detectorImages[name] += localImage

    def drawPointSource(self, x_pupil=None, y_pupil=None, bandpass=None):
        """
        Draw an image of a point source.

        @param [in] x_pupil is the x pupil coordinate of the object in arc seconds

        @param [in] y_pupil is the y pupil coordinate of the objec tin arc seconds

        @param [in] bandpass is an instantiation of the galsim.Bandpass class characterizing
        the bandpass over which we are integrating (in case the PSF is wavelength dependent)
        """

        if self.PSF is None:
            raise RuntimeError("Cannot draw a point source in GalSim without a PSF")

        return self.PSF.applyPSF(x_pupil=x_pupil, y_pupil=y_pupil, bandpass=bandpass)

    def drawSersic(self, x_pupil=None, y_pupil=None, sindex=None, minorAxis=None,
                   majorAxis=None, positionAngle=None, halfLightRadius=None, bandpass=None):
        """
        Draw the image of a Sersci profile.

        @param [in] x_pupil is the x pupil coordinate of the object in arc seconds

        @param [in] y_pupil is the y pupil coordinate of the object in arc seconds

        @param [in] sindex is the Sersic index of the object

        @param [in] minorAxis is the semi-minor axis of the object in any units (we only care
        about the ratio of the semi-minor to semi-major axes)

        @param [in] majorAxis is the semi-major axis of the object in the same units
        as minorAxis

        @param [in] halfLightRadius is the half light radius of the object in arc seconds

        @param [in] bandpass is an instantiation of the galsim.Bandpass class characterizing
        the bandpass over which we are integrating (in case the PSF is wavelength dependent)
        """

        #create a Sersic profile
        centeredObj = galsim.Sersic(n=float(sindex), half_light_radius=float(halfLightRadius))

        #turn the Sersic profile into an ellipse
        centeredObj = centeredObj.shear(q=minorAxis/majorAxis, beta=positionAngle*galsim.radians)
        if self.PSF is not None:
            centeredObj = self.PSF.applyPSF(x_pupil=x_pupil, y_pupil=y_pupil, obj=centeredObj,
                                            bandpass=bandpass)

        return centeredObj

    def addNoise(self, noiseWrapper=None, obs_metadata=None):
        """
        Adds a GalSim noise model to the images being stored by this
        GalSimInterpreter
        
        @param [in] noiseWrapper is a wrapper for GalSim's noise models
        (see ExampleCCDNoise for an example)
        
        @param [in] obs_metadata is an ObservationMetaData instantiation
        """
        
        for detector in detectorList:
            for bandPassName in self.bandPasses:
                name = self._getFileName(detector=detector, bandPassName=bandPassName)
                if name in self.detectorImages:
                    noiseModel = noiseWrapper.getNoiseModel(obs_metadata=obs_metadata, detector=detector)
                    self.detectorImages[name].addNoise(noiseModel)

    def writeImages(self, nameRoot=None):
        """
        Write the FITS files to disk.

        @param [in] nameRoot is a string that will be prepended to the names of the output
        FITS files.  The files will be named like

        nameRoot_detectorName_bandpassName.fits

        myImages_R_0_0_S_1_1_y.fits is an example of an image for an LSST-like camera with
        nameRoot = 'myImages'
        """
        for name in self.detectorImages:
            if nameRoot is not None:
                fileName = nameRoot+'_'+name
            else:
                fileName = name
            self.detectorImages[name].write(file_name=fileName)
