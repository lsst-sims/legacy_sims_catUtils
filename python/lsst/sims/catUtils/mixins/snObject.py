"""
Class describing the SN object itself. The SN object derives from
SNCosmo.Model and provides a sed of SNIa based on the SALT2 model-like
 model often called 'salt2-extended'. This model is extended to longer
 wavelength ranges compared to models in Guy10 or Betoule14, at the cost of
 larger model variance. SNObject has additional attributes such as ra, dec. and
additional methods to calculate band magnitudes using the LSST software stack
after applying MW extinction:
 -  calc_mags which use the magnitude calculations in LSST stack
 -  extinction which use the extinction calculations in LSST stack

"""
import numpy as np

from lsst.sims.photUtils.Sed import Sed
from lsst.sims.photUtils.EBV import EBVbase
from lsst.sims.photUtils.BandpassDict import BandpassDict
from lsst.sims.photUtils.SignalToNoise import calcSNR_m5, calcMagError_m5
from lsst.sims.photUtils.PhotometricParameters import PhotometricParameters

import sncosmo

__all__ = ['SNObject']


class SNObject(sncosmo.Model):

    """
    Extension of the SNCosmo `TimeSeriesModel` to include more parameters and
    use methods in the catsim stack. We constrain ourselves to the use of a
    specific SALT model for the Supernova (Salt2-Extended), and set its MW
    extinction to be 0, since we will use the LSST software to calculate
    extinction.

    Parameters
    ----------
    ra : float
         ra of the SN in degrees
    dec : float
        dec of the SN in degrees


    Attributes
    ----------
    _ra : float or None
        ra of the SN in radians

    _dec : float or None
        dec of the SN in radians

    skycoord : `np.ndarray' of size 2 or None
        np.array([[ra], [dec]]), which are in radians

    ebvofMW : float or None
        mwebv value calculated from the self.skycoord if not None, or set to
        a value using self.set_MWebv. If neither of these are done, this value
        will be None, leading to exceptions in extinction calculation.
        Therefore, the value must be set explicitly to 0. to get unextincted
        quantities.

    Methods
    -------

    Examples
    --------
    >>> SNObject  = SNObject(ra=30., dec=60.)
    >>> SNObject._ra
    >>> 0.5235987755982988
    >>> SNObject._dec
    >>> 1.0471975511965976
    """

    def __init__(self, ra=None, dec=None, source='salt2-extended'):
        """
        Parameters
        ----------
        ra : float
            ra of the SN in degrees
        dec : float
            dec of the SN in degrees
        """

        dust = sncosmo.OD94Dust()
        sncosmo.Model.__init__(self, source=source, effects=[dust, dust],
                               effect_names=['host', 'mw'],
                               effect_frames=['rest', 'obs'])

        # Current implementation of Model has a default value of mwebv = 0.
        # ie. no extinction, but this is not part of the API, so should not
        # depend on it, set explicitly in order to unextincted SED from
        # SNCosmo. We will use catsim extinction from `lsst.sims.photUtils`.

        self.ModelSource = source
        self.set(mwebv=0.)

        # ra and dec if passed are assumed to be in degrees and converted into
        # radians.
        self._ra = ra
        self._dec = dec

        # NB: More lines of code to support the possibility that ra, dec are
        # not provided at instantiation, and default to None
        self._hascoords = True
        if self._dec is None:
            self._hascoords = False
        if self._ra is None:
            self._hascoords = False

        # Satisfied that coordinates provided are floats
        if self._hascoords:
            self.setCoords(ra, dec)

        # For values of ra, dec, the E(B-V) is calculated directly
        # from DustMaps
        self.lsstmwebv = EBVbase()
        self.ebvofMW = None
        if self._hascoords:
            self.mwEBVfromMaps()
        return

    @property
    def SNstate(self):
        """
        Dictionary summarizing the property of the state. Can be used to
        serialize to disk

        Returns : Dictionary with values of parameters of the model.
        """
        statedict = dict()

        # SNCosmo Parameters
        statedict['ModelSource'] = self.source.name
        for param_name in self.param_names:
            statedict[param_name] = self.get(param_name)

        # New Attributes
        # statedict['lsstmwebv'] = self.lsstmwebv
        statedict['_ra'] = self._ra
        statedict['_dec'] = self._dec
        statedict['MWE(B-V)'] = self.ebvofMW

        return statedict

    @staticmethod
    def equivsncosmoParamDict(SNstate, SNCosmoModel):
        """
        return a dictionary that contains the parameters of SNCosmoModel
        that are contained in SNstate

        Parameters
        ----------
        SNstate : `SNObject.SNstate`, mandatory
            Dictionary defining the state of a SNObject
        SNCosmoModel : A `sncosmo.Model` instance, mandatory

        Returns
        -------
        sncosmoParams: Dictionary of sncosmo parameters

        """
        sncosmoParams = dict()
        for param in SNstate.keys():
            if param in SNCosmoModel.param_names:
                sncosmoParams[param] = SNstate[param]
        sncosmoParams['mwebv'] = SNstate['MWE(B-V)']
        return sncosmoParams

    @staticmethod
    def sncosmoParamDict(SNstate, SNCosmoModel):
        """
        return a dictionary that contains the parameters of SNCosmoModel
        that are contained in SNstate

        Parameters
        ----------
        SNstate : `SNObject.SNstate`, mandatory
            Dictionary defining the state of a SNObject
        SNCosmoModel : A `sncosmo.Model` instance, mandatory

        Returns
        -------
        sncosmoParams: Dictionary of sncosmo parameters

        """
        sncosmoParams = dict()
        for param in SNstate.keys():
            if param in SNCosmoModel.param_names:
                sncosmoParams[param] = SNstate[param]
        return sncosmoParams

    @classmethod
    def fromSNState(cls, snState):
        """
        creates an instance of SNObject with a state described by snstate.

        Parameters
        ----------
        snState: Dictionary summarizing the state of SNObject

        Returns
        -------
        Instance of SNObject class with attributes set by snstate

        Example
        -------

        """
        # Separate into SNCosmo parameters and SNObject parameters
        dust = sncosmo.OD94Dust()
        sncosmoModel = sncosmo.Model(source=snState['ModelSource'],
                                     effects=[dust, dust],
                                     effect_names=['host', 'mw'],
                                     effect_frames=['rest', 'obs'])

        sncosmoParams = cls.sncosmoParamDict(snState, sncosmoModel)

        # Now create the class
        cls = SNObject(source=snState['ModelSource'])

        # Set the SNObject coordinate properties
        # Have to be careful to not convert `None` type objects to degrees
        setdec, setra = False, False
        if snState['_ra'] is not None:
            ra = np.degrees(snState['_ra'])
            setra = True
        if snState['_dec'] is not None:
            dec = np.degrees(snState['_dec'])
            setdec = True
        if setdec and setra:
            cls.setCoords(ra, dec)

        # Set the SNcosmo parameters
        cls.set(**sncosmoParams)

        # Set the ebvofMW by hand
        cls.ebvofMW = snState['MWE(B-V)']

        return cls

    def equivalentSNCosmoModel(self):
        """
        returns an SNCosmo Model which is equivalent to SNObject
        """
        snState = self.SNstate
        dust = sncosmo.OD94Dust()
        sncosmoModel = sncosmo.Model(source=snState['ModelSource'],
                                     effects=[dust, dust],
                                     effect_names=['host', 'mw'],
                                     effect_frames=['rest', 'obs'])

        sncosmoParams = self.sncosmoParamDict(snState, sncosmoModel)
        sncosmoParams['mwebv'] = snState['MWE(B-V)']
        sncosmoModel.set(**sncosmoParams)
        return sncosmoModel

    def summary(self):
        '''
        summarizes the current state of the SNObject class in a returned
        string.

        Parameters
        ----------
        None

        Returns
        -------
        Summary State in string format

        Examples
        --------
        >>> t = SNObject()
        >>> print t.summary()
        '''
        state = '  SNObject Summary      \n'

        state += 'Model = ' + '\n'
        state += 'z = ' + str(self.get('z')) + '\n'
        state += 'c = ' + str(self.get('c')) + '\n'
        state += 'x1 = ' + str(self.get('x1')) + '\n'
        state += 'x0 = ' + str(self.get('x0')) + '\n'
        state += 't0 = ' + str(self.get('t0')) + '\n'
        state += 'ra = ' + str(self._ra) + ' in radians \n'
        state += 'dec = ' + str(self._dec) + ' in radians \n'
        state += 'MW E(B-V) = ' + str(self.ebvofMW) + '\n'

        return state

    def setCoords(self, ra, dec):
        """
        set the ra and dec coordinate of SNObject to values in radians
        corresponding to the given values in degrees

        Parameters
        ----------
        ra: float, mandatory
            the ra in degrees
        dec: float, mandatory
            dec in degrees

        Returns
        -------
        None

        Examples
        --------
        >>> t = SNObject()
        >>> t.setCoords(ra=30., dec=90.)
        >>> t._ra
        >>> 0.5235987755982988
        >>> t._dec
        >>> 1.0471975511965976
        """

        if ra is None or dec is None:
            raise ValueError('Why try to set coordinates without full'
                             'coordiantes?\n')
        self._ra = np.radians(ra)
        self._dec = np.radians(dec)
        self.skycoord = np.array([[self._ra], [self._dec]])

        self._hascoords = True

        return

    def set_MWebv(self, value):
        """
        if mwebv value is known, this can be used to set the attribute
        ebvofMW of the SNObject class to the value (float).

        Parameters
        ----------
        value: float, mandatory
               value of mw extinction parameter E(B-V) in mags to be used in
               applying extinction to the SNObject spectrum

        Returns
        -------
        None

        Examples
        --------
        >>> t = SNObject()
        >>> t.set_MWebv(0.)
        >>> 0.

        .. note:: For a large set of SN, one may use fast `np.ndarray` valued
                  functions to obtain an array of such values, and then set
                  the values from such an array.
        """
        self.ebvofMW = value
        return

    def mwEBVfromMaps(self):
        """
        set the attribute ebvofMW of the class from the ra and dec
        of the SN. If the ra or dec attribute of the class is None,
        set this attribute to None.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Examples
        --------
        >>> t = SNObject()
        >>> t.setCoords(ra=30., dec=60.)
        >>> t.mwEBVfromMaps()
        >>> t.ebvofMW
        >>> 0.977767825127

        .. note:: This function must be run after the class has attributes ra
                  and dec set. In case it is run before this, the mwebv value
                  will be set to None.

        """

        if not self._hascoords:
            raise ValueError('Cannot Calculate EBV from dust maps if ra or dec'
                             'is `None`')
        self.ebvofMW = self.lsstmwebv.calculateEbv(
            equatorialCoordinates=self.skycoord)[0]
        return

    def SNObjectSED(self, time, wavelen=None, bandpass=None,
                    applyExtinction=True):
        '''
        return a `lsst.sims.photUtils.sed` object from the SN model at the
        requested time and wavelengths with or without extinction from MW
        according to the SED extinction methods. The wavelengths may be
        obtained from a `lsst.sims.Bandpass` object or a `lsst.sims.BandpassDict`
        object instead. (Currently, these have the same wavelengths). See notes
        for details on handling of exceptions.

        If the sed is requested at times outside the validity range of the
        model, the flux density is returned as 0. If the time is within the
        range of validity of the model, but the wavelength range requested
        is outside the range, then the returned fluxes are np.nan outside
        the range, and the model fluxes inside

        Parameters
        ----------
        time: float
            time of observation
        wavelen: `np.ndarray` of floats, optional, defaults to None
            array containing wavelengths in nm
        bandpass: `lsst.sims.photUtils.Bandpass` object or
            `lsst.sims.photUtils.BandpassDict`, optional, defaults to `None`.
            Using the dict assumes that the wavelength sampling and range
            is the same for all elements of the dict.

            if provided, overrides wavelen input and the SED is
            obtained at the wavelength values native to bandpass
            object.


        Returns
        -------
        `sims_photutils.sed` object containing the wavelengths and SED
        values from the SN at time time in units of ergs/cm^2/sec/nm


        .. note: If both wavelen and bandpassobject are `None` then exception,
                 will be raised.
        Examples
        --------
        >>> sed = SN.SNObjectSED(time=0., wavelen=wavenm)
        '''

        if wavelen is None and bandpass is None:
            raise ValueError('A non None input to either wavelen or\
                              bandpassobject must be provided')

        # if bandpassobject present, it overrides wavelen
        if bandpass is not None:
            if isinstance(bandpass, BandpassDict):
                firstfilter = bandpass.keys()[0]
                bp = bandpass[firstfilter]
            else:
                bp = bandpass
            # remember this is in nm
            wavelen = bp.wavelen

        flambda = np.zeros(len(wavelen))

        # self.mintime() and self.maxtime() are properties describing
        # the ranges of SNCosmo.Model in time.

        # Set SED to 0 beyond the model phase range, will change this if
        # SNCosmo includes a more sensible decay later.
        if (time > self.mintime()) & (time < self.maxtime()):

            # If SNCosmo is requested a SED value beyond the model range
            # it will crash. Try to prevent that by returning np.nan for
            # such wavelengths. This will still not help band flux calculations
            # but helps us get past this stage.

            flambda = flambda * np.nan

            # Convert to Ang
            wave = wavelen * 10.0
            mask1 = wave > self.minwave()
            mask2 = wave < self.maxwave()
            mask = mask1 & mask2
            wave = wave[mask]

            # flux density dE/dlambda returned from SNCosmo in
            # ergs/cm^2/sec/Ang, convert to ergs/cm^2/sec/nm

            flambda[mask] = self.flux(time=time, wave=wave)
            flambda[mask] = flambda[mask] * 10.0

        SEDfromSNcosmo = Sed(wavelen=wavelen, flambda=flambda)

        if not applyExtinction:
            return SEDfromSNcosmo

        # Apply LSST extinction
        ax, bx = SEDfromSNcosmo.setupCCMab()

        if self.ebvofMW is None:
            raise ValueError('ebvofMW attribute cannot be None Type and must'
                             ' be set by hand using set_MWebv before this'
                             'stage, or by using setcoords followed by'
                             'mwEBVfromMaps\n')

        SEDfromSNcosmo.addCCMDust(a_x=ax, b_x=bx, ebv=self.ebvofMW)
        return SEDfromSNcosmo

    def catsimBandFluxes(self, time, bandpassobject):
        """
        return the flux in the bandpass in units of the flux
        the AB magnitude reference spectrum would have in the
        same band.

        Parameters
        ----------
        time: mandatory, float
            MJD at which band fluxes are evaluated
        bandpassobject: mandatory, `lsst.sims.photUtils.BandPass` object
            A particular bandpass which is an instantiation of one of
            (u, g, r, i, z, y)
        Returns
        -------
        float

        Examples
        --------
        .. note: If there is an unphysical value of sed in
        the wavelength range, it produces a flux of  `np.nan`
        """
        SEDfromSNcosmo = self.SNObjectSED(time=time,
                                          bandpass=bandpassobject)
        return SEDfromSNcosmo.calcFlux(bandpass=bandpassobject) / 3631.0

    def catsimBandMagError(self, time, bandpassobject, m5, photParams=None,
                           magnitude=None):
        """
        return the mag uncertainty in the bandpass

        Parameters
        ----------
        time: mandatory, float
            MJD at which band fluxes are evaluated
        bandpassobject: mandatory, `lsst.sims.photUtils.BandPass` object
            A particular bandpass which is an instantiation of one of
            (u, g, r, i, z, y)
        m5 :
        photParams :
        magnitude :
        Returns
        -------
        float

        Examples
        --------
        .. note: If there is an unphysical value of sed in
        the wavelength range, it produces a flux of  `np.nan`
        """

        if magnitude is None:
            mag = self.catsimBandMags(time=time, bandpassobject=bandpassobject)
        else:
            mag = magnitude

        # mag = np.asarray([[mag]])
        bandpass = bandpassobject

        if photParams is None:
            photParams = PhotometricParameters()

        # m5 = np.asarray([[m5]])

        magerr = calcMagError_m5(magnitude=mag, bandpass=bandpass, m5=m5,
                                 photParams=photParams)
        return magerr[0]

    def catsimBandFluxError(self, time, bandpassobject, m5, photParams=None,
                            magnitude=None):
        """
        return the flux uncertainty in the bandpass in units 'maggies'
        (the flux the AB magnitude reference spectrum would have in the
        same band.)

        Parameters
        ----------
        time: mandatory, float
            MJD at which band fluxes are evaluated
        bandpassobject: mandatory, `lsst.sims.photUtils.BandPass` object
            A particular bandpass which is an instantiation of one of
            (u, g, r, i, z, y)
        m5 :
        photParams :
        magnitude :
        Returns
        -------
        float

        Examples
        --------
        .. note: If there is an unphysical value of sed in
        the wavelength range, it produces a flux of  `np.nan`
        """
        # SEDfromSNcosmo = self.SNObjectSED(time=time,
        #                                 bandpass=bandpassobject)

        if magnitude is None:
            mag = self.catsimBandMags(time=time, bandpassobject=bandpassobject)
        else:
            mag = magnitude
        flux = self.catsimBandFluxes(time=time, bandpassobject=bandpassobject)

        # mag = np.asarray([[mag]])
        bandpass = bandpassobject

        if photParams is None:
            photParams = PhotometricParameters()

        # m5 = np.asarray([[m5]])
        SNR, gamma = calcSNR_m5(magnitude=mag, bandpass=bandpass, m5=m5,
                                photParams=photParams)
        return flux / SNR

    def catsimManyBandFluxes(self, time, bandpassDict,
                             observedBandPassInd=None):
        """
        return the flux in the bandpass in units of the flux
        the AB magnitude reference spectrum would have in the
        same band.

        Parameters
        ----------
        time: mandatory, float
            MJD at which band fluxes are evaluated
        bandpassDict: mandatory, `lsst.sims.photUtils.BandpassDict` instance
        observedBandPassInd : optional, list of integers, defaults to None
            integer correspdonding to index of the bandpasses used in the
            observation in the ordered dict bandpassDict
        Returns
        -------
        `~numpy.ndarray` of length =len(observedBandPassInd)

        Examples
        --------
        .. note: If there is an unphysical value of sed in
        the wavelength range, it produces a flux of  `np.nan`
        """
        SEDfromSNcosmo = self.SNObjectSED(time=time,
                                          bandpass=bandpassDict['u'])
        wavelen_step = np.diff(SEDfromSNcosmo.wavelen)[0]
        SEDfromSNcosmo.flambdaTofnu()
        f = SEDfromSNcosmo.manyFluxCalc(bandpassDict.phiArray,
                                        wavelen_step=wavelen_step,
                                        observedBandpassInd=observedBandPassInd)
        return f / 3631.

    def catsimManyBandMags(self, time, bandpassDict,
                           observedBandPassInd=None):
        """
        return the flux in the bandpass in units of the flux
        the AB magnitude reference spectrum would have in the
        same band.

        Parameters
        ----------
        time: mandatory, float
            MJD at which band fluxes are evaluated
        bandpassDict: mandatory, `lsst.sims.photUtils.BandpassDict` instance
        observedBandPassInd : optional, list of integers, defaults to None
            integer correspdonding to index of the bandpasses used in the
            observation in the ordered dict bandpassDict
        Returns
        -------
        `~numpy.ndarray` of length =len(observedBandPassInd)

        Examples
        --------
        .. note: If there is an unphysical value of sed in
        the wavelength range, it produces a flux of  `np.nan`
        """
        f = self.catsimManyBandFluxes(time,
                                      bandpassDict,
                                      observedBandPassInd)

        return -2.5 * np.log10(f)

    def catsimBandMags(self, bandpassobject, time):
        """
        return the magnitude in the bandpass in the AB magnitude system

        Parameters
        ----------
        bandpassobject : mandatory,`sims.photUtils.BandPass` instances
            LSST Catsim bandpass instance for the bandpass
        time : mandatory, float
              MJD at which this is evaluated
        Returns
        -------
            float
        Examples
        --------
        """
        fluxRatio = self.catsimBandFluxes(bandpassobject=bandpassobject,
                                          time=time)
        return -2.5 * np.log10(fluxRatio)

    def catsimADU(self, time, bandpassDict,
                  photParams=None,
                  observedBandPassInds=None):
        """
        time: float, mandatory
            MJD of the observation

        bandpassDict: mandatory,
            Dictionary of instances of `sims.photUtils.Bandpass` for
            filters

        photParams: Instance of `sims.photUtils.PhotometricParameters`, optional,
                    defaults to None
                    Describes the observational parameters used in specifying the
                    photometry of the ovservation
        observedBandPassInd: None
            Not used now
        """
        SEDfromSNcosmo = self.SNObjectSED(time=time,
                                          bandpass=bandpassDict)

        bandpassNames = bandpassDict.keys()
        adus = np.zeros(len(bandpassNames))

        for i, filt in enumerate(bandpassNames):
            bandpass = bandpassDict[filt]
            adus[i] = SEDfromSNcosmo.calcADU(bandpass, photParams=photParams)

        return adus
