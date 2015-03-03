"""
This file defines the following classes:

GalSimInterpreter -- a class which takes objects passed by a GalSim Instance Catalog
(see galSimCatalogs.py) and uses GalSim to write them to FITS images.

GalSimDetector -- a class which stored information about a detector in a way that
GalSimInterpreter expects.
"""

import os
import numpy
import galsim
from lsst.sims.utils import radiansToArcsec

__all__ = ["GalSimInterpreter", "GalSimDetector"]

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
        self.localFluxes = []

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

    def _doesObjectImpingeOnDetector(self, xPupil=None, yPupil=None, detector=None, imgScale=None, nonZeroPixels=None):
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
        this method will return True if the object is within 120 arc seconds of the detector's
        bounds (on the assumption that no PSF will smear the object outmore than 2 arc seconds)
        """
        if detector is None:
            return False
       
        for (ix, iy) in zip(nonZeroPixels[0], nonZeroPixels[1]):
            xx = xPupil + ix*imgScale
            yy = yPupil + iy*imgScale
            if xx<detector.xMax and xx>detector.xMin and yy<detector.yMax and yy>detector.yMin:
                return True

        return False

    def findAllDetectors(self, galSimType=None, x_pupil=None,
                   y_pupil=None, halfLightRadius=None, minorAxis=None, majorAxis=None,
                   positionAngle=None, sindex=None):

        outputString = ''
        outputList = []
        centeredObj = None
        testScale = 0.1
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

            centeredImage = centeredObj.drawImage(scale=testScale, method='phot', n_photons=1000)
            xmax = testScale * (centeredImage.getXMax()/2) + xp
            xmin = testScale * (-1*centeredImage.getXMax()/2) + xp
            ymax = testScale * (centeredImage.getYMax()/2) + yp
            ymin = testScale *(-1*centeredImage.getYMin()/2) + yp
            
            viableDetectors = []
            for dd in self.detectors:
                xOverLaps = False
                if xmax > dd.xMin and xmax < dd.xMax:
                    xOverLaps = True
                elif xmin > dd.xMin and xmin < dd.xMax:
                    xOverLaps = True
                elif xmin < dd.xMin and xmax > dd.xMax:
                    xOverLaps = True
                
                yOverLaps = False
                if ymax > dd.yMin and ymax < dd.yMax:
                    yOverLaps = True
                elif ymin > dd.yMin and ymin < dd.yMax:
                    yOverLaps = True
                elif ymin < dd.yMin and ymax > dd.yMax:
                    yOverLaps = True
                
                if xOverLaps and yOverLaps and dd not in outputList:
                    viableDetectors.append(dd)

            
            if len(viableDetectors)>0:
                maxPixel = centeredImage(centeredImage.getXMax()/2, centeredImage.getYMax()/2)
                activePixels = numpy.where(centeredImage.array>maxPixel*0.001)
            
                xmin = testScale * (activePixels[0].min() - centeredImage.getXMax()/2) + xp
                xmax = testScale * (activePixels[0].max() - centeredImage.getXMax()/2) + xp
                ymin = testScale * (activePixels[1].min() - centeredImage.getYMax()/2) + yp
                ymax = testScale * (activePixels[1].max() - centeredImage.getYMax()/2) + yp
            
                for dd in viableDetectors:
                    xOverLaps = False
                    if xmax > dd.xMin and xmax < dd.xMax:
                        xOverLaps = True
                    elif xmin > dd.xMin and xmin < dd.xMax:
                        xOverLaps = True
                    elif xmin < dd.xMin and xmax > dd.xMax:
                        xOverLaps = True
                
                    yOverLaps = False
                    if ymax > dd.yMin and ymax < dd.yMax:
                        yOverLaps = True
                    elif ymin > dd.yMin and ymin < dd.yMax:
                        yOverLaps = True
                    elif ymin < dd.yMin and ymax > dd.yMax:
                        yOverLaps = True
                
                    if xOverLaps and yOverLaps:
                        if self._doesObjectImpingeOnDetector(xPupil=xp - centeredImage.getXMax()*testScale/2.0,
                                                 yPupil=yp - centeredImage.getYMax()*testScale/2.0,
                                                 detector=dd, imgScale=centeredImage.scale,
                                                 nonZeroPixels=activePixels):

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
                
                self.localFluxes.append(localImage.array.sum())
                
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

        for detector in self.detectors:
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
            total = 0.0
            for ix in range(self.detectorImages[name].getXMax()):
                for iy in range(self.detectorImages[name].getYMax()):
                    total += self.detectorImages[name](ix+1,iy+1)
            print name,': ',total

        self.localFluxes.sort()
        for f in self.localFluxes:
            print f
