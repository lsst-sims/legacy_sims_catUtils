"""
This file defines the model classes that wrap noise models from
galsim into the CatSim interface
"""

import numpy
import galsim
from lsst.sims.utils import radiansToArcsec

__all__ = ["ExampleCCDNoise"]

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


