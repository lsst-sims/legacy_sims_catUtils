"""
This file defines the model classes that wrap noise models from
galsim into the CatSim interface
"""

import numpy
import galsim
from lsst.sims.photUtils import expectedSkyCountsForM5, PhotometricDefaults

__all__ = ["ExampleCCDNoise"]

class NoiseAndBackgroundBase(object):

    def __init__(self, seed=None):
        """
        @param [in] seed is an (optional) int that will seed the
        random number generator used by the noise model
        """

        if seed is None:
            self.randomNumbers = galsim.UniformDeviate()
        else:
            self.randomNumbers = galsim.UniformDeviate(seed)

        self.addedBackground = False


    def getNoiseModel(self, *args, **kwargs):
        raise NotImplementedError("There is no noise model for NoiseAndBackgroundBase")


    def addNoiseAndBackground(self, image, bandpass=None, m5=None,
                              addBackground=True,
                              addNoise = True,
                              expTime=PhotometricDefaults.exptime,
                              nexp=PhotometricDefaults.nexp,
                              readnoise=PhotometricDefaults.rdnoise,
                              darkcurrent=PhotometricDefaults.darkcurrent,
                              othernoise=PhotometricDefaults.othernoise,
                              seeing=PhotometricDefaults.seeing['r'],
                              platescale=PhotometricDefaults.platescale,
                              gain=PhotometricDefaults.gain,
                              effarea=PhotometricDefaults.effarea):
    
        skyCounts = expectedSkyCountsForM5(m5, bandpass,
                                           expTime=expTime, nexp=nexp,
                                           readnoise=readnoise, darkcurrent=darkcurrent,
                                           othernoise=othernoise, seeing=seeing,
                                           platescale=platescale, gain=gain,
                                           effarea=effarea)

        image = image.copy()

        if addBackground:
            image += skyCounts
            skyLevel = 0.0
        else:
            skyLevel = skyCounts

        if addNoise:
            noiseModel = self.getNoiseModel(skyLevel=skyLevel, gain=gain, readNoise=readnoise)
            image.addNoise(noiseModel)

        return image
        

class ExampleCCDNoise(NoiseAndBackgroundBase):
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

    def getNoiseModel(self, skyLevel=0.0, gain=PhotometricDefaults.gain,
                      readNoise=PhotometricDefaults.rdnoise):

        return galsim.CCDNoise(self.randomNumbers, sky_level=skyLevel, gain=gain, read_noise=readNoise)


