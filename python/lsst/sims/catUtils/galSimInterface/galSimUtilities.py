"""
This file defines the model classes that wrap PSFs and noise models from
galsim into the CatSim interface
"""

import numpy
import galsim
from lsst.sims.catalogs.generation.db import radiansToArcsec

__all__ = ["PSFbase", "ExampleGaussianPSF", "ExampleOpticalPSF",
           "ExampleCCDNoise"]

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


