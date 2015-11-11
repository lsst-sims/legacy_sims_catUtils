#!/usr/bin/env python

# from cStringIO import StringIO
import numpy as np

from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.measures.instance import compound
from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.utils import ObservationMetaData

from lsst.sims.photUtils import Sed
from lsst.sims.photUtils import BandpassDict
from lsst.sims.catUtils.mixins import CosmologyMixin
from lsst.sims.catUtils.mixins import PhotometryBase
import lsst.sims.photUtils.PhotometricParameters as PhotometricParameters
from lsst.sims.photUtils.SignalToNoise import calcSNR_m5

import astropy
import sncosmo

from .snObject import SNObject
from .snUniversalRules import SNUniverse

import sqlite3

wavelenstep = 0.1
cosmo = CosmologyMixin()


class SNIaCatalog (InstanceCatalog, CosmologyMixin, SNUniverse):

    """
    `lsst.sims.catalogs.measures.instance.InstanceCatalog` class with SN
    characterized by the  following attributes

    Attributes
    ----------
    position : 3-tuple of floats
              (ra, dec, redshift),
    velocity : 3 tuple of floats
              velocity wrt host galaxy in Km/s,
    the supernova model (eg. SALT2)
    and parameters of the supernova model that predict the SED.
    """

    # t_0, c, x_1, x_0 are parameters characterizing a SALT
    # based SN model as defined in sncosmo
    column_outputs = ['snid', 'snra', 'sndec', 'z', 't0', 'c', 'x1',
                      'x0']

    # You can add parameters like fluxes and magnitudes by adding the following
    # variables to the list
    # 'flux_u', 'flux_g', 'flux_r', 'flux_i', 'flux_z', 'flux_y' ,
    # 'mag_u', 'mag_g', 'mag_r', 'mag_i', 'mag_z', 'mag_y']

    # Flux variables are convenient to display in exponential format to avoid
    # having them cut off
    variables = ['flux_u', 'flux_g', 'flux_r', 'flux_i', 'flux_z', 'flux_y']
    variables += ['flux', 'flux_err', 'mag_err']

    override_formats = {'snra': '%8e', 'sndec': '%8e', 'c': '%8e',
            'x0': '%8e'} #, 'flux': '%8e', 'mag': '%8e', 'flux_err': '%8e'}
    for var in variables:
        override_formats[var] = '%8e'

    cannot_be_null = ['x0', 'z', 't0']

    @astropy.utils.lazyproperty
    def mjdobs(self):
        '''
        The time of observation for the catalog, which is set to be equal
        to obs_metadata.mjd
        '''
        return self.obs_metadata.mjd

    @astropy.utils.lazyproperty
    def badvalues(self):
        '''
        The representation of bad values in this catalog is numpy.nan
        '''
        return np.nan

    @property
    def suppressDimSN(self):
        '''
        Boolean to decide whether to output observations of SN that are too dim
        should be represented in the catalog or not. By default set to True
        '''
        # self._suppressDimSN = True
        if not hasattr(self, '_suppressDimSN'):
            suppressDimSN_default = True
            self._suppressDimSN = suppressDimSN_default
        return self._suppressDimSN

    @suppressDimSN.setter
    def suppressDimSN(self, suppressDimSN):
        '''
        set the value of suppressDimSN of the catalog 

        Parameters
        ----------
        value : Boolean, mandatory
            Value to set suppressDimSN to 
        '''
        # if suppressDimSN is None:
        #    self._suppressDimSN = True
        # else:
        self._suppressDimSN = suppressDimSN
        return self._suppressDimSN

    @property
    def suppressHighzSN(self):
        '''
        Boolean to decide whether to output information corresponding to SN 
        at redshift above self.maxz 
        '''

        return True

    @property
    def midSurveyTime(self):
        '''
        The time at the middle of the survey: ie. at the 5 year period.


        .. note: Changing this should not change the statistical
        properties of the survey, but will change the exact SN we find.
        '''
        if not hasattr(self, '_midSurveyTime'):
            midSurveyTime_default = 570000.0
            self._midSurveyTime = midSurveyTime_default
        return self._midSurveyTime

    @midSurveyTime.setter
    def midSurveyTime(self, mymidSurveyTime):
        '''
        set the value of suppressDimSN of the catalog 

        Parameters
        ----------
        value : Boolean, mandatory
            Value to set suppressDimSN to 
        '''
        # if suppressDimSN is None:
        #    self._suppressDimSN = True
        # else:
        self._midSurveyTime = mymidSurveyTime
        return self._midSurveyTime

    @property
    def maxTimeSNVisible(self):
        '''
        The catalog will provide values for SN flux (even if zero according to
        model) for peak \pm maxTimeVisibilityforSN 
        '''
        return 100.

    @property
    def maxz(self):
        return 1.2

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

    @astropy.utils.lazyproperty
    def lsstpbase(self):
        pbase = PhotometryBase()
        return pbase

    def get_snid(self):
        # Not necessarily unique if the same galaxy hosts two SN
        return self.column_by_name('id')

    @property
    def numobjs(self):
        return len(self.column_by_name('id'))

    @compound('snra', 'sndec', 'z', 'vra', 'vdec', 'vr')
    def get_angularCoordinates(self):
        '''
        Obtain the coordinates and velocity of the SN from the host galaxy

        Returns
        -------
        `np.ndarray` of coordinara, dec, z, vra, vdec, and vr

        '''
        hostra , hostdec, hostz = self.column_by_name('raJ2000'),\
            self.column_by_name('decJ2000'),\
            self.column_by_name('redshift')
        snra, sndec, snz, snvra, snvdec, snvr = self.SNCoordinatesFromHost(
            hostra, hostdec, hostz)

        return ([snra, sndec, snz, snvra, sndec, snvr])

    @compound('c', 'x1', 'x0', 't0')
    def get_snparams(self):
        hostz, hostid, hostmu = self.column_by_name('redshift'),\
            self.column_by_name('snid'),\
            self.column_by_name('cosmologicalDistanceModulus')

        vals = self.SNparamDistfromHost(hostz, hostid, hostmu)

        return (vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3])

    def get_SNsed(self):
        """
        returns a list of SN seds in `lsst.sims.photUtils.Sed` observed within
        the spatio-temporal range specified by obs_metadata

        """
        c, x1, x0, t0, _z, ra, dec = self.column_by_name('c'),\
            self.column_by_name('x1'),\
            self.column_by_name('x0'),\
            self.column_by_name('t0'),\
            self.column_by_name('redshift'),\
            self.column_by_name('raJ2000'),\
            self.column_by_name('decJ2000')

        SNobject = SNObject()

        sedlist = []
        for i in range(self.numobjs):
            SNobject.set(z=_z[i], c=c[i], x1=x1[i], t0=t0[i], x0=x0[i])
            SNobject.setCoords(ra=ra[i], dec=dec[i])
            SNobject.mwEBVfromMaps()
            sed = SNobject.SNObjectSED(time=self.obs_metadata.mjd,
                                       bandpass=lsstBandpassDict,
                                       applyExitinction=True)
            sedlist.append(sed)

        return sedlist

    def get_time(self):
        
        return np.repeat(self.mjdobs, self.numobjs)

    def get_band(self):
        bandname = self.obs_metadata.bandpass
        return np.repeat(bandname, self.numobjs)

    def get_fluxerr(self):

        bandname = self.obs_metadata.bandpass
        bandpass = self.lsstBandpassDict[bandname]
        vals = np.zeros(self.numobjs)
        flux = self.column_by_name('flux')
        mag = self.column_by_name('mag')
        photParams = PhotometricParameters()
        m5 = self.obs_metadata.m5[bandname]
        m5 = np.asarray([[m5]])

        SNR, gamma = calcSNR_m5(magnitudes=[mag],
                                bandpasses=[bandpass],
                                m5=m5,
                                photParams=photParams)

        return flux

    @compound('flux', 'mag', 'flux_err', 'mag_err')
    def get_snbrightness(self):

        c, x1, x0, t0, _z , _id, ra, dec = self.column_by_name('c'),\
            self.column_by_name('x1'),\
            self.column_by_name('x0'),\
            self.column_by_name('t0'),\
            self.column_by_name('redshift'),\
            self.column_by_name('snid'),\
            self.column_by_name('raJ2000'),\
            self.column_by_name('decJ2000')

        SNobject = SNObject()
        bandname = self.obs_metadata.bandpass
        bandpass = self.lsstBandpassDict[bandname]

        # Initialize return array
        vals = np.zeros(shape=(self.numobjs, 4))

        # vals = np.zeros(self.numobjs)
        for i, v in enumerate(vals):
            arr = [_z[i], c[i], x1[i], t0[i], x0[i]]
            testnan = lambda x: x is np.nan
            SNobject.set(z=_z[i], c=c[i], x1=x1[i], t0=t0[i], x0=x0[i])
            SNobject.setCoords(ra=ra[i], dec=dec[i])
            SNobject.mwEBVfromMaps()

            # Calculate fluxes
            flux = SNobject.catsimBandFluxes(time=self.mjdobs,
                                             bandpassobject=bandpass)
            mag = SNobject.catsimBandMags(time=self.mjdobs,
                                          bandpassobject=bandpass)
            vals[i, 0] = flux
            vals[i, 1] = mag
            flux_err = SNobject.catsimBandFluxError(time=self.mjdobs,
                                                     bandpassobject=bandpass,
                                                     m5=self.obs_metadata.m5[bandname],
                                                     photParams=None,
                                                     magnitude=mag)
            mag_err = SNobject.catsimBandMagError(time=self.mjdobs,
                                                     bandpassobject=bandpass,
                                                     m5=self.obs_metadata.m5[bandname],
                                                     photParams=None,
                                                     magnitude=mag)
            vals[i, 2] = flux_err
            vals[i, 3] = mag_err
        return (vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3])


    @compound('flux_u', 'flux_g', 'flux_r', 'flux_i', 'flux_z', 'flux_y',
              'mag_u', 'mag_g', 'mag_r', 'mag_i', 'mag_z', 'mag_y',
              'adu_u', 'adu_g', 'adu_r', 'adu_i', 'adu_z', 'adu_y')
    def get_snfluxes(self):

        c, x1, x0, t0, _z , _id, ra, dec = self.column_by_name('c'),\
            self.column_by_name('x1'),\
            self.column_by_name('x0'),\
            self.column_by_name('t0'),\
            self.column_by_name('redshift'),\
            self.column_by_name('snid'),\
            self.column_by_name('raJ2000'),\
            self.column_by_name('decJ2000')

        SNobject = SNObject()
        # Initialize return array
        vals = np.zeros(shape=(self.numobjs, 18))
        for i, v in enumerate(vals):
            arr = [_z[i], c[i], x1[i], t0[i], x0[i]]
            testnan = lambda x: x is np.nan
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

            vals[i, 12:] = SNobject.catsimADU(time=self.obs_metadata,
                                              bandpassDict=self.lsstBandpassDict,
                                              photParams=self.photometricparameters)
        return (vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3],
                vals[:, 4], vals[:, 5], vals[:, 6], vals[:, 7],
                vals[:, 8], vals[:, 9], vals[:, 10], vals[:, 11],
                vals[:, 12], vals[:, 13], vals[:, 14], vals[:, 15],
                vals[:, 16], vals[:, 17])
