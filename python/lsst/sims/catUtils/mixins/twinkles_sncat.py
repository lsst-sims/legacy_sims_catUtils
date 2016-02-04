"""
Mixin to InstanceCatalog class to give SN catalogs in catsim
"""
import numpy as np

from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.measures.instance import compound
from lsst.sims.photUtils import BandpassDict
from lsst.sims.catUtils.mixins import CosmologyMixin
from lsst.sims.catUtils.mixins import PhotometryBase
import lsst.sims.photUtils.PhotometricParameters as PhotometricParameters

from lsst.sims.catUtils.exampleCatalogDefinitions.phoSimCatalogExamples import\
        PhoSimCatalogSN
import astropy

from lsst.sims.catUtils.supernovae import SNObject
from lsst.sims.catUtils.supernovae import SNUniverse
from lsst.sims.catUtils.mixins import SNFunctionality


__all__ = ['TwinkSNCat']
cosmo = CosmologyMixin()



class TwinkSNCat ( SNFunctionality,  InstanceCatalog, CosmologyMixin, SNUniverse):

    """
    `lsst.sims.catalogs.measures.instance.InstanceCatalog` class with SN
    characterized by the  following attributes

    Attributes
    ----------
    column_outputs :
    suppressHighzSN :
    maxTimeSNVisible :
    maxz :
    variables :
    override_formats :
    cannot_be_null :
    mjdobs :
    badvalues position :
    3-tuple of floats (ra, dec, redshift), velocity : 3 tuple of floats
        velocity wrt host galaxy in Km/s, the supernova model (eg. SALT2)
    and parameters of the supernova model that predict the SED.
    """

    # t_0, c, x_1, x_0 are parameters characterizing a SALT
    # based SN model as defined in sncosmo
    column_outputs = ['snid', 'snra', 'sndec', 'z', 't0', 'c', 'x1', 'x0']

    surveyStartDate = 59580. # For Kraken_1042
    suppressHighzSN = True
    maxTimeSNVisible = 100.
    maxz = 1.2
    # Flux variables are convenient to display in exponential format to avoid
    # having them cut off
    variables = ['flux_u', 'flux_g', 'flux_r', 'flux_i', 'flux_z', 'flux_y']
    variables += ['flux', 'flux_err', 'mag_err']

    override_formats = {'snra': '%8e', 'sndec': '%8e', 'c': '%8e',
                        'x0': '%8e'}
    for var in variables:
        override_formats[var] = '%8e'
    # You can add parameters like fluxes and magnitudes by adding the following
    # variables to the list
    # 'flux_u', 'flux_g', 'flux_r', 'flux_i', 'flux_z', 'flux_y' ,
    # 'mag_u', 'mag_g', 'mag_r', 'mag_i', 'mag_z', 'mag_y']
    cannot_be_null = ['x0', 'z', 't0']

    @astropy.utils.lazyproperty
    def mjdobs(self):
        '''
        The time of observation for the catalog, which is set to be equal
        to obs_metadata.mjd
        '''
        return self.obs_metadata.mjd.TAI

    @astropy.utils.lazyproperty
    def badvalues(self):
        '''
        The representation of bad values in this catalog is numpy.nan
        '''
        return np.nan

    @property
    def suppressDimSN(self):
        """
        Boolean to decide whether to output observations of SN that are too dim
        should be represented in the catalog or not. By default set to True
        """
        if not hasattr(self, '_suppressDimSN'):
            suppressDimSN_default = True
            self._suppressDimSN = suppressDimSN_default
        return self._suppressDimSN

    @suppressDimSN.setter
    def suppressDimSN(self, suppressDimSN):
        """
        set the value of suppressDimSN of the catalog Parameters
        Parameters
        ----------
        supressDimSN : Boolean, mandatory
            Value to set suppressDimSN to
        """
        self._suppressDimSN = suppressDimSN
        return self._suppressDimSN 

    @astropy.utils.lazyproperty
    def photometricparameters(self, expTime=15., nexp=2):
        lsstPhotometricParameters = PhotometricParameters(exptime=expTime,
                                                          nexp=nexp)
        return lsstPhotometricParameters

    @astropy.utils.lazyproperty
    def lsstBandpassDict(self):
        return BandpassDict.loadTotalBandpassesFromFiles()

    @astropy.utils.lazyproperty
    def observedIndices(self):
        bandPassNames = self.obs_metadata.bandpass
        return [self.lsstBandpassDict.keys().index(x) for x in bandPassNames]

    def get_snid(self):
        # Not necessarily unique if the same galaxy hosts two SN
        # Use refIdCol to access the relevant id column of the dbobj
        # Should revert to galTileID for galaxyTiled catalogDBObj and
        # id for galaxyObj catalogDBObj
        # (email from Scott)
        return self.column_by_name('Tsnid')

    @property
    def numobjs(self):
        return len(self.column_by_name('Tsnid'))

    @compound('snra', 'sndec', 'z', 'vra', 'vdec', 'vr')
    def get_angularCoordinates(self):
        '''
        Obtain the coordinates and velocity of the SN from the host galaxy

        Returns
        -------
        `np.ndarray` of coordinara, dec, z, vra, vdec, and vr

        '''
        snra, sndec, snz = self.column_by_name('raJ2000'),\
            self.column_by_name('decJ2000'),\
            self.column_by_name('Tredshift')
        snvra = np.zeros(self.numobjs)
        snvdec = np.zeros(self.numobjs)
        snvr = np.zeros(self.numobjs)

        return (snra, sndec, snz, snvra, snvdec, snvr)

    @compound('c', 'x1', 'x0', 't0')
    def get_snparams(self):

        c, x1, x0 = self.column_by_name('Tc'), \
                    self.column_by_name('Tx1'),\
                    self.column_by_name('Tx0')
        t0 = self.column_by_name('Tt0') + self.surveyStartDate
        if self.suppressDimSN :
            print('badvalues ', self.badvalues)
            print('mjd ', self.mjdobs)
            print('maxTime', self.maxTimeSNVisible)
            print('number of cases ',
                    len(t0[np.abs(t0 - self.mjdobs) > self.maxTimeSNVisible]))
            t0 = np.where(np.abs(t0 - self.mjdobs) > self.maxTimeSNVisible,
                      self.badvalues, t0)

        return (c, x1, x0, t0)


    def get_time(self):

        return np.repeat(self.mjdobs, self.numobjs)

    def get_band(self):
        bandname = self.obs_metadata.bandpass
        return np.repeat(bandname, self.numobjs)

