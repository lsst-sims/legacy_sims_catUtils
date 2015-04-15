"""
This file defines the model classes that wrap noise models from
galsim into the CatSim interface
"""

import numpy
import galsim
from lsst.sims.photUtils import expectedSkyCountsForM5, PhotometricDefaults

__all__ = ["ExampleCCDNoise"]

class NoiseAndBackgroundBase(object):
    """
    This is the base class for all wrappers of sky background
    and noise in the GalSim interface.  Daughters of this class
    are meant to be included as the noise_and_background
    class member variable of GalSim InstanceCatalog classes.
    To implement a new noise model, users should write a new class
    that inherits from this one.  That new class should only define
    a method getNoiseModel() that takes as arguments skyLevel, gain,
    and readNoise.  See the docstring for getNoiseModel() for further
    details.
    """

    def __init__(self, seed=None):
        """
        @param [in] seed is an (optional) int that will seed the
        random number generator used by the noise model
        """

        if seed is None:
            self.randomNumbers = galsim.UniformDeviate()
        else:
            self.randomNumbers = galsim.UniformDeviate(seed)


    def getNoiseModel(self, skyLevel=0.0, readNoise=None, gain=None):
        """
        This method returns the noise model implemented for this wrapper
        class.

        @param [in] skyLevel is the number of electrons per pixel due
        to the sky background.  However, this value should only be non-zero
        if the sky background has been subtracted from the image.  The
        purpose of this parameter is to provide an extra background value
        when calculating the level of Poisson noise in each pixel.  If the
        sky background is already present in the image, then the noise model
        will just set the noise level based on the intensity in each pixel
        and there is no need to add an additional skyLevel.  If the sky
        background is still included in the image, set skyLevel equal to zero.

        @param [in] readNoise is the read noise in electrons per pixel

        @param [in] gain is the number of electrons per ADU

        @param [out] returns an instantiation of a GalSim noise class, as
        specified by the particular wrapper class to which this method belongs.
        """

        raise NotImplementedError("There is no noise model for NoiseAndBackgroundBase")


    def addNoiseAndBackground(self, image, bandpass=None, m5=None,
                              addBackground=True,
                              addNoise = True,
                              expTime=PhotometricDefaults.exptime,
                              nexp=1,
                              readnoise=PhotometricDefaults.rdnoise,
                              darkcurrent=PhotometricDefaults.darkcurrent,
                              othernoise=PhotometricDefaults.othernoise,
                              seeing=PhotometricDefaults.seeing['r'],
                              platescale=PhotometricDefaults.platescale,
                              gain=PhotometricDefaults.gain,
                              effarea=PhotometricDefaults.effarea):
        """
        This method actually adds the sky background and noise to an image.

        @param [in] image is the GalSim image object to which the background
        and noise are being added.

        @param [in] bandpass is a CatSim bandpass object (not a GalSim bandpass
        object) characterizing the filter through which the image is being taken.

        @param [in] addBackground = True if you want to actually add the sky background
        to the image; False otherwise (default True).

        @param [in] addNoise = True if you want to add noise to the image; False if
        otherwise (default True).

        @param [in] expTime is the exposure time in seconds.

        @param [in] nexp is the number of exposures (default 1).

        @param [in] readnoise is the readnoise in electrons per pixel.

        @param [in] darkcurrent

        @param [in] othernoise additional systematic noise.

        @param [in] seeing is the seeing in arcseconds.

        @param [in] platescale is the number of arcseconds per pixel.

        @param [in] gain is the number of electrons per ADU.

        @param [in] effarea is the effective area of the telescope's primary mirror.

        @param [out] the input image with the background and noise model added to it.
        """


        #calculate the sky background to be added to each pixel
        skyCounts = expectedSkyCountsForM5(m5, bandpass,
                                           expTime=expTime, nexp=nexp,
                                           readnoise=readnoise, darkcurrent=darkcurrent,
                                           othernoise=othernoise, seeing=seeing,
                                           platescale=platescale, gain=gain,
                                           effarea=effarea)

        image = image.copy()

        if addBackground:
            image += skyCounts
            skyLevel = 0.0 #if we are adding the skyCounts to the image,there is no need
                           #to pass a skyLevel parameter to the noise model.  skyLevel is
                           #just used to calculate the level of Poisson noise.  If the
                           #sky background is included in the image, the Poisson noise
                           #will be calculated from the actuall image brightness.
        else:
            skyLevel = skyCounts/gain

        if addNoise:
            noiseModel = self.getNoiseModel(skyLevel=skyLevel, gain=gain, readNoise=readnoise)
            image.addNoise(noiseModel)

        return image


class ExampleCCDNoise(NoiseAndBackgroundBase):
    """
    This class wraps the GalSim class CCDNoise.  It is meant to be assigned as
    the self.noise_and_background member variable in a GalSim InstanceCatalog.
    To instantiatea different noise model, write a class like this one that
    defines a method getNoiseModel() which accepts as its arguments skyLevel,
    readNoise, and gain and returns an instantiation of a GalSim noise model
    """

    def getNoiseModel(self, skyLevel=0.0, gain=PhotometricDefaults.gain,
                      readNoise=PhotometricDefaults.rdnoise):

        """
        This method returns the noise model implemented for this wrapper
        class.

        @param [in] skyLevel is the number of electrons per pixel due
        to the sky background.  However, this value should only be non-zero
        if the sky background has been subtracted from the image.  The
        purpose of this parameter is to provide an extra background value
        when calculating the level of Poisson noise in each pixel.  If the
        sky background is already present in the image, then the noise model
        will just set the noise level based on the intensity in each pixel
        and there is no need to add an additional skyLevel.  If the sky
        background is still included in the image, set skyLevel equal to zero.

        @param [in] readNoise is the read noise in electrons per pixel

        @param [in] gain is the number of electrons per ADU

        @param [out] returns an instantiation of the GalSim CCDNoise class
        """

        return galsim.CCDNoise(self.randomNumbers, sky_level=skyLevel, gain=gain, read_noise=readNoise)


