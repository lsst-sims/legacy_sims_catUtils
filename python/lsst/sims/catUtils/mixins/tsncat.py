"""
Mixin to InstanceCatalog class to give SN catalogs in catsim
"""
import numpy as np
from copy import deepcopy

from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.measures.instance import compound
from lsst.sims.photUtils import BandpassDict
from lsst.sims.catUtils.mixins import CosmologyMixin
from lsst.sims.catUtils.mixins import PhotometryBase
import lsst.sims.photUtils.PhotometricParameters as PhotometricParameters

import astropy

from lsst.sims.catUtils.supernovae import SNObject
from lsst.sims.catUtils.supernovae import SNUniverse


__all__ = ['TSNIaCatalog']
cosmo = CosmologyMixin()


class TSNIaCatalog (InstanceCatalog, CosmologyMixin, SNUniverse):

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
    column_outputs = ['Tsnid', 'snra', 'sndec', 'Tredshift', 'Tt0', 'Tc', 'Tx1', 'Tx0']
    #column_outputs = ['snid', 'snra', 'sndec', 'z', 't0', 'c', 'x1', 'x0']
    # The t0 value stored in the database is in terms of MJD - survey start
    # date. survey start date must be stored so that times correspond to
    # ObservationMetaData times
    # This can be reset by the user
    surveyStartDate = 59580. # For Kraken_1042


    # You can add parameters like fluxes and magnitudes by adding the following
    # variables to the list
    # 'flux_u', 'flux_g', 'flux_r', 'flux_i', 'flux_z', 'flux_y' ,
    # 'mag_u', 'mag_g', 'mag_r', 'mag_i', 'mag_z', 'mag_y']

    suppressHighzSN = True
    maxTimeSNVisible = 100.
    maxz = 1.2
    # Flux variables are convenient to display in exponential format to avoid
    # having them cut off
    variables = ['flux_u', 'flux_g', 'flux_r', 'flux_i', 'flux_z', 'flux_y']
    variables += ['flux', 'flux_err', 'mag_err']

    override_formats = {'snra': '%8e', 'sndec': '%8e', 'Tc': '%8e',
                        'Tx0': '%8e'}
    for var in variables:
        override_formats[var] = '%8e'

    cannot_be_null = ['Tx0', 'Tredshift', 'Tt0']

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

    @property
    def numobjs(self):
        return len(self.column_by_name('snid'))
        
    @compound('Tt0', 'TTx1')
    def get_snparams(self):
        tt0 = deepcopy(self.column_by_name('Tt0'))
        tt0 += self.surveyStartDate
        ttx1 = self.column_by_name('Tx1')

        return(tt0, ttx1) 

    def get_time(self):
        return np.repeat(self.mjdobs, self.numobjs)

    def get_band(self):
        bandname = self.obs_metadata.bandpass
        return np.repeat(bandname, self.numobjs)

    @compound('flux', 'mag', 'flux_err', 'mag_err')
    def get_snbrightness(self):

        c, x1, x0, t0, _z, ra, dec = self.column_by_name('Tc'),\
            self.column_by_name('Tx1'),\
            self.column_by_name('Tx0'),\
            self.column_by_name('Tt0'),\
            self.column_by_name('Tredshift'),\
            self.column_by_name('raJ2000'),\
            self.column_by_name('decJ2000')

        SNobject = SNObject()
        bandname = self.obs_metadata.bandpass
        print('Printing the value of obs_metadat.bandpass', bandname)
        if isinstance(bandname, list):
            raise ValueError('bandname expected to be string, but is list\n')
        bandpass = self.lsstBandpassDict[bandname]

        # Initialize return array
        vals = np.zeros(shape=(self.numobjs, 4))

        for i, _ in enumerate(vals):
            SNobject.set(z=_z[i], c=c[i], x1=x1[i], t0=t0[i], x0=x0[i])
            SNobject.setCoords(ra=ra[i], dec=dec[i])
            SNobject.mwEBVfromMaps()

            # Calculate fluxes
            fluxinMaggies = SNobject.catsimBandFlux(time=self.mjdobs,
                                                    bandpassobject=bandpass)
            mag = SNobject.catsimBandMag(time=self.mjdobs,
                                         fluxinMaggies=fluxinMaggies,
                                         bandpassobject=bandpass)
            vals[i, 0] = fluxinMaggies
            vals[i, 1] = mag
            flux_err = SNobject.catsimBandFluxError(time=self.mjdobs,
                                                    bandpassobject=bandpass,
                                                    m5=self.obs_metadata.m5[
                                                        bandname],
                                                    photParams=None,
                                                    fluxinMaggies=fluxinMaggies,
                                                    magnitude=mag)

            mag_err = SNobject.catsimBandMagError(time=self.mjdobs,
                                                  bandpassobject=bandpass,
                                                  m5=self.obs_metadata.m5[
                                                      bandname],
                                                  photParams=None,
                                                  magnitude=mag)
            vals[i, 2] = flux_err
            vals[i, 3] = mag_err
        return (vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3])

    @compound('flux_u', 'flux_g', 'flux_r', 'flux_i', 'flux_z', 'flux_y',
              'mag_u', 'mag_g', 'mag_r', 'mag_i', 'mag_z', 'mag_y',
              'adu_u', 'adu_g', 'adu_r', 'adu_i', 'adu_z', 'adu_y', 'mwebv')
    def get_snfluxes(self):

        c, x1, x0, t0, _z, ra, dec = self.column_by_name('Tc'),\
            self.column_by_name('Tx1'),\
            self.column_by_name('Tx0'),\
            self.column_by_name('Tt0'),\
            self.column_by_name('Tredshift'),\
            self.column_by_name('raJ2000'),\
            self.column_by_name('decJ2000')

        SNobject = SNObject()
        # Initialize return array
        vals = np.zeros(shape=(self.numobjs, 19))
        for i, _ in enumerate(vals):
            SNobject.set(z=_z[i], c=c[i], x1=x1[i], t0=t0[i], x0=x0[i])
            SNobject.setCoords(ra=ra[i], dec=dec[i])
            SNobject.mwEBVfromMaps()
            # Calculate fluxes
            vals[i, :6] = SNobject.catsimManyBandFluxes(time=self.mjdobs,
                                                        bandpassDict=self.lsstBandpassDict,
                                                        observedBandPassInd=None)
            # Calculate magnitudes
            vals[i, 6:12] = SNobject.catsimManyBandMags(time=self.mjdobs,
                                                        bandpassDict=self.lsstBandpassDict,
                                                        observedBandPassInd=None)

            vals[i, 12:18] = SNobject.catsimManyBandADUs(time=self.mjdobs,
                                                bandpassDict=self.lsstBandpassDict,
                                                photParams=self.photometricparameters)
            vals[i, 18] = SNobject.ebvofMW
        return (vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3],
                vals[:, 4], vals[:, 5], vals[:, 6], vals[:, 7],
                vals[:, 8], vals[:, 9], vals[:, 10], vals[:, 11],
                vals[:, 12], vals[:, 13], vals[:, 14], vals[:, 15],
                vals[:, 16], vals[:, 17], vals[:, 18])
