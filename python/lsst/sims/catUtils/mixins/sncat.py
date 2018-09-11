"""
Mixins for the InstanceCatalog class to provide SN catalogs in catsim. There
are three classes here:
    - SNFunctionality: provides common functions required by all SN instance
        catalogs. It does not make sense to instantiate this class, but rather
        it should be used as a mixin alongside another class.
    - SNIaCatalog: Dynamically created catalogs by sampling user specified
        distributions of SN parameters on the fly based on host galaxies in the
        catsim database.
    - FrozenSNCat: catalogs that are 'frozen' on the catsim database. For a
        user to use one of these catalogs, such a catalog would have to be
        uploaded to catsim. Examples of such catalogs that are in the catsim
        database are the tables `TwinkSN` and `TwinkSNKraken`
"""
from builtins import str
from builtins import range
import numpy as np
import numbers

from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catalogs.decorators import compound
from lsst.sims.photUtils import (BandpassDict, Bandpass)
from lsst.sims.catUtils.mixins import CosmologyMixin
import lsst.sims.photUtils.PhotometricParameters as PhotometricParameters
from lsst.sims.catUtils.supernovae import SNObject
from lsst.sims.catUtils.supernovae import SNUniverse
from lsst.sims.catUtils.mixins import EBVmixin
from lsst.sims.utils import _galacticFromEquatorial
import astropy


__all__ = ['SNIaCatalog', 'SNFunctionality', 'FrozenSNCat']
cosmo = CosmologyMixin()


class SNFunctionality(InstanceCatalog, EBVmixin, CosmologyMixin, SNUniverse):
    """
    SNFunctionality is a mixin that provides functionality of getting fluxes
    and magnitudes for SN defined by parameters of `~sims_catUtils.SNObject` as
    defined in `~sims_catUtils/python/lsst/sims/catUtils/supernovae/SNObject`


    This class is not meant to be used by itself, as it does not have any way
    of obtaining its attributes, but as a mixin to classes like SNIaCatalog
    which define these attributes.
    """

    # Write the location of SED file (for example for PhoSim)
    writeSedFile = False
    # prefix to use for SED File name
    sn_sedfile_prefix = ''

    # t_0, c, x_1, x_0 are parameters characterizing a SALT
    # based SN model as defined in sncosmo
    column_outputs = ['snid', 'snra', 'sndec', 'z', 't0', 'c', 'x1', 'x0']

    _lsstmwebv = None
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

    _sn_object_cache = None

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

    @compound('TsedFilepath', 'magNorm')
    def get_phosimVars(self):
        """
        Obtain variables sedFilepath to be used to obtain unique filenames
        for each SED for phoSim and MagNorm which is also used. Note that aside
        from acting as a getter, this also writes spectra to 
        `self.sn_sedfile_prefix`snid_mjd_band.dat for each observation of
        interest
        """
        # construct the unique filename
        # method: snid_mjd(to 4 places of decimal)_bandpassname
        mjd = "_{:0.4f}_".format(self.mjdobs)
        mjd += self.obs_metadata.bandpass + '.dat'
        fnames = np.array([self.sn_sedfile_prefix + str(int(elem)) + mjd
                           if isinstance(elem, numbers.Number) else
                           self.sn_sedfile_prefix + str(elem) + mjd
                           for elem in self.column_by_name('snid')], dtype='str')

        c, x1, x0, t0, z = self.column_by_name('c'),\
            self.column_by_name('x1'),\
            self.column_by_name('x0'),\
            self.column_by_name('t0'),\
            self.column_by_name('redshift')

        bp = Bandpass()
        bp.imsimBandpass()

        magNorms = np.zeros(len(fnames))

        snobject = SNObject()
        snobject.rectifySED = True
        for i in range(len(self.column_by_name('snid'))):
            # if t0 is nan, this was set by the catalog for dim SN, or SN
            #   outside redshift range, We will not provide a SED file for these
            if np.isnan(t0[i]):
                magNorms[i] = np.nan
                fnames[i] = None

            else:
                snobject.set(c=c[i], x1=x1[i], x0=x0[i], t0=t0[i],
                             z=z[i])
                if snobject.modelOutSideTemporalRange == 'zero':
                    if self.mjdobs > snobject.maxtime() or self.mjdobs < snobject.mintime():
                        magNorms[i] = np.nan
                        fnames[i] = None

                # SED in rest frame
                sed = snobject.SNObjectSourceSED(time=self.mjdobs)
                try:
                    magNorms[i] = sed.calcMag(bandpass=bp)
                except:
                    # sed.flambda = 1.0e-20
                    magNorms[i] = 1000.  # sed.calcMag(bandpass=bp)

                if self.writeSedFile:
                    sed.writeSED(fnames[i])

        return (fnames, magNorms)

    def get_snid(self):
        # Not necessarily unique if the same galaxy hosts two SN
        # Use refIdCol to access the relevant id column of the dbobj
        # Should revert to galTileID for galaxyTiled catalogDBObj and
        # id for galaxyObj catalogDBObj
        # (email from Scott)
        return self.column_by_name(self.refIdCol)

    def load_SNsed(self):
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

        raDeg = np.degrees(ra)
        decDeg = np.degrees(dec)

        sedlist = []
        for i in range(self.numobjs):
            SNobject.set(z=_z[i], c=c[i], x1=x1[i], t0=t0[i], x0=x0[i])
            SNobject.setCoords(ra=raDeg[i], dec=decDeg[i])
            SNobject.mwEBVfromMaps()
            sed = SNobject.SNObjectSED(time=self.mjdobs,
                                       bandpass=self.lsstBandpassDict,
                                       applyExitinction=True)
            sedlist.append(sed)

        return sedlist

    @property
    def numobjs(self):
        return len(self.column_by_name('id'))

    def get_time(self):
        """
        mjd at SALT2 'peak'
        """
        return np.repeat(self.mjdobs, self.numobjs)

    def get_band(self):
        bandname = self.obs_metadata.bandpass
        return np.repeat(bandname, self.numobjs)

    @compound('flux', 'mag', 'flux_err', 'mag_err', 'adu')
    def get_snbrightness(self):
        """
        getters for brightness related parameters of sn
        """
        if self._sn_object_cache is None or len(self._sn_object_cache) > 1000000:
            self._sn_object_cache = {}

        c, x1, x0, t0, _z, ra, dec = self.column_by_name('c'),\
            self.column_by_name('x1'),\
            self.column_by_name('x0'),\
            self.column_by_name('t0'),\
            self.column_by_name('redshift'),\
            self.column_by_name('raJ2000'),\
            self.column_by_name('decJ2000')

        raDeg = np.degrees(ra)
        decDeg = np.degrees(dec)

        ebv = self.column_by_name('EBV')
        id_list = self.column_by_name('snid')

        bandname = self.obs_metadata.bandpass
        if isinstance(bandname, list):
            raise ValueError('bandname expected to be string, but is list\n')
        bandpass = self.lsstBandpassDict[bandname]

        # Initialize return array so that it contains the values you would get
        # if you passed through a t0=self.badvalues supernova
        vals = np.array([[0.0]*len(t0), [np.inf]*len(t0),
                        [np.nan]*len(t0), [np.inf]*len(t0),
                        [0.0]*len(t0)]).transpose()

        for i in np.where(np.logical_and(np.isfinite(t0),
                                         np.abs(self.mjdobs - t0) < self.maxTimeSNVisible))[0]:

            if id_list[i] in self._sn_object_cache:
                SNobject = self._sn_object_cache[id_list[i]]
            else:
                SNobject = SNObject()
                SNobject.set(z=_z[i], c=c[i], x1=x1[i], t0=t0[i], x0=x0[i])
                SNobject.setCoords(ra=raDeg[i], dec=decDeg[i])
                SNobject.set_MWebv(ebv[i])
                self._sn_object_cache[id_list[i]] = SNobject

            if self.mjdobs <= SNobject.maxtime() and self.mjdobs >= SNobject.mintime():

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
                                                        photParams=self.photometricparameters,
                                                        fluxinMaggies=fluxinMaggies,
                                                        magnitude=mag)

                mag_err = SNobject.catsimBandMagError(time=self.mjdobs,
                                                      bandpassobject=bandpass,
                                                      m5=self.obs_metadata.m5[
                                                          bandname],
                                                      photParams=self.photometricparameters,
                                                      magnitude=mag)
                sed = SNobject.SNObjectSED(time=self.mjdobs,
                                           bandpass=self.lsstBandpassDict,
                                           applyExtinction=True)
                adu = sed.calcADU(bandpass, photParams=self.photometricparameters)
                vals[i, 2] = flux_err
                vals[i, 3] = mag_err
                vals[i, 4] = adu

        return (vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3], vals[:, 4])

    @compound('flux_u', 'flux_g', 'flux_r', 'flux_i', 'flux_z', 'flux_y',
              'mag_u', 'mag_g', 'mag_r', 'mag_i', 'mag_z', 'mag_y',
              'adu_u', 'adu_g', 'adu_r', 'adu_i', 'adu_z', 'adu_y', 'mwebv')
    def get_snfluxes(self):

        c, x1, x0, t0, _z, ra, dec = self.column_by_name('c'),\
            self.column_by_name('x1'),\
            self.column_by_name('x0'),\
            self.column_by_name('t0'),\
            self.column_by_name('redshift'),\
            self.column_by_name('raJ2000'),\
            self.column_by_name('decJ2000')

        raDeg = np.degrees(ra)
        decDeg = np.degrees(dec)

        snobject = SNObject()
        # Initialize return array
        vals = np.zeros(shape=(self.numobjs, 19))
        for i, _ in enumerate(vals):
            snobject.set(z=_z[i], c=c[i], x1=x1[i], t0=t0[i], x0=x0[i])
            snobject.setCoords(ra=raDeg[i], dec=decDeg[i])
            snobject.mwEBVfromMaps()
            # Calculate fluxes
            vals[i, :6] = snobject.catsimManyBandFluxes(time=self.mjdobs,
                                                        bandpassDict=self.lsstBandpassDict,
                                                        observedBandPassInd=None)
            # Calculate magnitudes
            vals[i, 6:12] = snobject.catsimManyBandMags(time=self.mjdobs,
                                                        bandpassDict=self.lsstBandpassDict,
                                                        observedBandPassInd=None)

            vals[i, 12:18] = snobject.catsimManyBandADUs(time=self.mjdobs,
                                                         bandpassDict=self.lsstBandpassDict,
                                                         photParams=self.photometricparameters)
            vals[i, 18] = snobject.ebvofMW
        return (vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3],
                vals[:, 4], vals[:, 5], vals[:, 6], vals[:, 7],
                vals[:, 8], vals[:, 9], vals[:, 10], vals[:, 11],
                vals[:, 12], vals[:, 13], vals[:, 14], vals[:, 15],
                vals[:, 16], vals[:, 17], vals[:, 18])

    #def get_EBV(self):
    #    return self.column_by_name('EBV')


class SNIaCatalog (SNFunctionality):

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

    @compound('snra', 'sndec', 'z', 'vra', 'vdec', 'vr')
    def get_angularCoordinates(self):
        '''
        Obtain the coordinates and velocity of the SN from the host galaxy

        Returns
        -------
        `np.ndarray` of coordinara, dec, z, vra, vdec, and vr

        '''
        hostra, hostdec, hostz = self.column_by_name('raJ2000'),\
            self.column_by_name('decJ2000'),\
            self.column_by_name('redshift')
        snra, sndec, snz, snvra, snvdec, snvr = self.SNCoordinatesFromHost(
            hostra, hostdec, hostz)

        return ([snra, sndec, snz, snvra, snvdec, snvr])

    @compound('glon', 'glat')
    def get_galacticCoords(self):
        return _galacticFromEquatorial(self.column_by_name('snra'), self.column_by_name('sndec'))

    @compound('c', 'x1', 'x0', 't0')
    def get_snparams(self):
        hostz, hostid, hostmu = self.column_by_name('redshift'),\
            self.column_by_name('snid'),\
            self.column_by_name('cosmologicalDistanceModulus')

        vals = self.SNparamDistFromHost(hostz, hostid, hostmu)
        return (vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3])


class FrozenSNCat(SNFunctionality):

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

    surveyStartDate = 59580.0   # For Kraken_1042 / Minion_1016

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
        """
        getter for SN parameters (SALT2)
        """

        c, x1, x0 = self.column_by_name('Tc'), \
            self.column_by_name('Tx1'),\
            self.column_by_name('Tx0')
        t0 = self.column_by_name('Tt0') + self.surveyStartDate
        if self.suppressDimSN:
            t0 = np.where(np.abs(t0 - self.mjdobs) > self.maxTimeSNVisible,
                          self.badvalues, t0)

        return (c, x1, x0, t0)

