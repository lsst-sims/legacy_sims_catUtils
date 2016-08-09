import numpy
import astropy.cosmology as cosmology
import astropy.units as units
from lsst.sims.catalogs.decorators import cached
from lsst.sims.photUtils import CosmologyObject
flatnessthresh = 1.0e-12

__all__ = ["CosmologyMixin"]

class CosmologyMixin(object):
    """
    This class is designed to operate as a mixin for InstanceCatalog classes.
    It provides a member variable self.cosmology which is an instantiation
    of the CosmologyObject class.  self.cosmology defaults to the
    Milliennium Simulation cosmology.  This mixin also provides a method
    self.setCosmology() which will allow the user to customize self.cosmology
    and a getter for the column cosmologicalDistanceModulus that reflects the
    effect of the luminosity distance on a galaxy's component magnitudes.

    NOTE: one should only include this mixin in catalogs whose magNorm is
    normalized to an absolute magnitude of some sort (i.e. the magnitude if
    the galaxy was at redshift=0).  The magNorms for galaxies stored on the
    University of Washington LSSTCATSIM database do not fit this criterion.
    magNorms on the University of Washington LSSTCATSIM database include the
    effects of cosmological distance modulus.
    """

    cosmology = CosmologyObject()

    def setCosmology(self, H0=73.0, Om0=0.25, Ok0=None, w0=None, wa=None):
        """
        This method customizes the member variable self.cosmology by re-instantiating
        the CosmologyObject class

        param [in] H0 is the Hubble parameter today in km/s/Mpc

        param [in] Om0 is the density paramter (fraction of critical) associated with matter today

        param [in] Ode0 is the density paratmer associated with dark energy today

        param [in] w0 is the w0 parameter associated with the equation of state of dark energy
        w = w0 + wa z/(1+z)

        param [in] wa is the wa parameter usesd to set the equation of state of dark energy
        w = w0 + wa z/(1+z)
        """
        self.cosmology = CosmologyObject(H0=H0, Om0=Om0, Ok0=Ok0, w0=w0, wa=wa)

    @cached
    def get_cosmologicalDistanceModulus(self):
        """
        getter for cosmologicalDistanceModulus (the effect of the luminosity
        distance on a galaxy's component magnitudes)
        """
        redshift = self.column_by_name("redshift")

        if len(redshift) == 0:
            #newer versions of astropy do not appreciate being passed an
            #empty numpy array of redshifts; avoid nasty exceptions by
            #just returning an empty numpy array if we got an empty numpy array
            return numpy.array([])

        return self.cosmology.distanceModulus(redshift)
