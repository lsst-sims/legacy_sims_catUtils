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
from builtins import str
import numpy as np

from lsst.sims.photUtils.Sed import Sed
from lsst.sims.catUtils.dust import EBVbase
from lsst.sims.photUtils.BandpassDict import BandpassDict
from lsst.sims.photUtils.SignalToNoise import calcSNR_m5, calcMagError_m5
from lsst.sims.photUtils.PhotometricParameters import PhotometricParameters

import sncosmo

__all__ = ['SNObject']


_sn_ax_cache = None
_sn_bx_cache = None
_sn_ax_bx_wavelen = None

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
    rectifySED : Bool, True by Default
        if the SED is negative at the requested time and wavelength, return 0.
        instead of the negative value.


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
        Instantiate object

        Parameters
        ----------
        ra : float
            ra of the SN in degrees
        dec : float
            dec of the SN in degrees

        source : instance of `sncosmo.SALT2Source`, optional, defaults to using salt2-extended 
            source class to define the model
        """

        dust = sncosmo.CCM89Dust()
        sncosmo.Model.__init__(self, source=source, effects=[dust, dust],
                               effect_names=['host', 'mw'],
                               effect_frames=['rest', 'obs'])

        # Current implementation of Model has a default value of mwebv = 0.
        # ie. no extinction, but this is not part of the API, so should not
        # depend on it, set explicitly in order to unextincted SED from
        # SNCosmo. We will use catsim extinction from `lsst.sims.photUtils`.

        self.ModelSource = source
        self.set(mwebv=0.)

        # self._ra, self._dec is initialized as None for cases where ra, dec
        # is not provided
        self._ra = None
        self._dec = None

        # ra, dec is input in degrees
        # If provided, set _ra, _dec in radians
        self._hascoords = True
        if dec is None:
            self._hascoords = False
        if ra is None:
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


        # Behavior of model outside temporal range :
        # if 'zero' then all fluxes outside the temporal range of the model
        # are set to 0.
        self._modelOutSideTemporalRange = 'zero'
        
        # SED will be rectified to 0. for negative values of SED if this
        # attribute is set to True
        self.rectifySED = True
        return

    @property
    def SNstate(self):
        """
        Dictionary summarizing the state of SNObject. Can be used to
        serialize to disk, and create SNObject from SNstate

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
        dust = sncosmo.CCM89Dust()
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

    @property
    def modelOutSideTemporalRange(self):
        """
        Defines the behavior of the model when sampled at times beyond the model
        definition.
        """
        return self._modelOutSideTemporalRange

    @modelOutSideTemporalRange.setter
    def modelOutSideTemporalRange(self, value):
        if value != 'zero':
            raise ValueError('Model not implemented, defaulting to zero method\n')
        return self._modelOutSideTemporalRange


    def equivalentSNCosmoModel(self):
        """
        returns an SNCosmo Model which is equivalent to SNObject
        """
        snState = self.SNstate
        dust = sncosmo.CCM89Dust()
        sncosmoModel = sncosmo.Model(source=snState['ModelSource'],
                                     effects=[dust, dust],
                                     effect_names=['host', 'mw'],
                                     effect_frames=['rest', 'obs'])

        sncosmoParams = self.sncosmoParamDict(snState, sncosmoModel)
        sncosmoParams['mwebv'] = snState['MWE(B-V)']
        sncosmoModel.set(**sncosmoParams)
        return sncosmoModel


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
        for param in SNstate:
            if param in SNCosmoModel.param_names:
                sncosmoParams[param] = SNstate[param]
        sncosmoParams['mwebv'] = SNstate['MWE(B-V)']
        return sncosmoParams

    @staticmethod
    def sncosmoParamDict(SNstate, SNCosmoModel):
        """
        return a dictionary that contains the parameters of SNCosmoModel
        that are contained in SNstate. Note that this does not return the
        equivalent SNCosmo  model.

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
        for param in SNstate:
            if param in SNCosmoModel.param_names:
                sncosmoParams[param] = SNstate[param]
        return sncosmoParams



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


    def redshift(self, z, cosmo):
        """
        Redshift the instance holding the intrinsic brightness of the object
        fixed. By intrinsic brightness here, we mean the BessellB band asbolute
        magnitude in rest frame. This requires knowing the cosmology


        Parameters
        ----------
        z : float, mandatory
            redshift at which the object must be placed.
        cosmo : instance of `astropy.cosmology` objects, mandatory
            specifies the cosmology.
        Returns
        -------
        None, but it changes the instance

        """
        import numbers

        # Check that the input redshift is a scalar
        try:
            assert isinstance(z, numbers.Number)
        except:
            raise TypeError('The argument z in method redshift should be'
                            'a scalar Numeric')

        # Ensure that the input redshift is greater than 0.
        try:
            assert z > 0.
        except:
            raise ValueError('The argument z in the method SNObject.redshift'
                             'should be greater than 0.')

        # Find the current value of the rest frame BessellB AB magnitude
        peakAbsMag = self.source_peakabsmag('BessellB', 'AB', cosmo=cosmo)
        self.set(z=z)
        self.set_source_peakabsmag(peakAbsMag, 'BessellB', 'AB', cosmo=cosmo)
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
        # the ranges of SNCosmo.Model in time. Behavior beyond this is 
        # determined by self.modelOutSideTemporalRange
        if (time >= self.mintime()) and (time <= self.maxtime()):
            # If SNCosmo is requested a SED value beyond the wavelength range
            # of model it will crash. Try to prevent that by returning np.nan for
            # such wavelengths. This will still not help band flux calculations
            # but helps us get past this stage.

            flambda = flambda * np.nan

            # Convert to Ang
            wave = wavelen * 10.0
            mask1 = wave >= self.minwave()
            mask2 = wave <= self.maxwave()
            mask = mask1 & mask2
            wave = wave[mask]

            # flux density dE/dlambda returned from SNCosmo in
            # ergs/cm^2/sec/Ang, convert to ergs/cm^2/sec/nm

            flambda[mask] = self.flux(time=time, wave=wave)
            flambda[mask] = flambda[mask] * 10.0
        else:
            # use prescription for modelOutSideTemporalRange
            if self.modelOutSideTemporalRange != 'zero':
                raise NotImplementedError('Model not implemented, change to zero\n')
                # Else Do nothing as flambda is already 0.
                # This takes precedence over being outside wavelength range
                
        if self.rectifySED:
            # Note that this converts nans into 0.
            flambda = np.where(flambda > 0., flambda, 0.)
        SEDfromSNcosmo = Sed(wavelen=wavelen, flambda=flambda)

        if not applyExtinction:
            return SEDfromSNcosmo

        # Apply LSST extinction
        global _sn_ax_cache
        global _sn_bx_cache
        global _sn_ax_bx_wavelen
        if _sn_ax_bx_wavelen is None \
        or len(wavelen)!=len(_sn_ax_bx_wavelen) \
        or (wavelen!=_sn_ax_bx_wavelen).any():

            ax, bx = SEDfromSNcosmo.setupCCM_ab()
            _sn_ax_cache = ax
            _sn_bx_cache = bx
            _sn_ax_bx_wavelen = np.copy(wavelen)
        else:
            ax = _sn_ax_cache
            bx = _sn_bx_cache

        if self.ebvofMW is None:
            raise ValueError('ebvofMW attribute cannot be None Type and must'
                             ' be set by hand using set_MWebv before this'
                             'stage, or by using setcoords followed by'
                             'mwEBVfromMaps\n')

        SEDfromSNcosmo.addDust(a_x=ax, b_x=bx, ebv=self.ebvofMW)
        return SEDfromSNcosmo

    def SNObjectSourceSED(self, time, wavelen=None):
        """
        Return the rest Frame SED of SNObject at the phase corresponding to
        time, at rest frame wavelengths wavelen. If wavelen is None,
        then the SED is sampled at the rest frame wavelengths native to the
        SALT model being used.

        Parameters
        ----------
        time : float, mandatory,
            observer frame time at which the SED has been requested in units
            of days.
        wavelen : `np.ndarray`, optional, defaults to native SALT wavelengths
            array of wavelengths in the rest frame of the supernova in units
            of nm. If None, this defaults to the wavelengths at which the
            SALT model is sampled natively.
        Returns
        -------
        `numpy.ndarray` of dtype float.

        .. note: The result should usually match the SALT source spectrum.
        However, it may be different for the following reasons:
        1. If the time of observation is outside the model range, the values
            have to be inserted using additional models. Here only one model
            is currently implemented, where outside the model range the value
            is set to 0.
        2. If the wavelengths are beyond the range of the SALT model, the SED
            flambda values are set to `np.nan` and these are actually set to 0.
            if `self.rectifySED = True`
        3. If the `flambda` values of the SALT model are negative which happens
            in the less sampled phases of the model, these values are set to 0,
            if `self.rectifySED` = True.  (Note: if `self.rectifySED` = True, then
            care will be taken to make sure that the flux at 500nm is not exactly
            zero, since that will cause PhoSim normalization of the SED to be
            NaN).
        """
        phase = (time - self.get('t0')) / (1. + self.get('z'))
        source = self.source

        # Set the default value of wavelength  input
        if wavelen is None:
            # use native SALT grid in Ang
            wavelen = source._wave
        else:
            #input wavelen in nm, convert to Ang
            wavelen = wavelen.copy()
            wavelen *= 10.0

        flambda = np.zeros(len(wavelen))
        # self.mintime() and self.maxtime() are properties describing
        # the ranges of SNCosmo.Model in time. Behavior beyond this is 
        # determined by self.modelOutSideTemporalRange
        insidephaseRange = (phase <= source.maxphase())and(phase >= source.minphase())
        if insidephaseRange:
            # If SNCosmo is requested a SED value beyond the wavelength range
            # of model it will crash. Try to prevent that by returning np.nan for
            # such wavelengths. This will still not help band flux calculations
            # but helps us get past this stage.

            flambda = flambda * np.nan

            mask1 = wavelen >= source.minwave()
            mask2 = wavelen <= source.maxwave()
            mask = mask1 & mask2
            # Where we have to calculate fluxes because it is not `np.nan`
            wave = wavelen[mask]
            flambda[mask] = source.flux(phase, wave)
        else:
            if self.modelOutSideTemporalRange == 'zero':
                # flambda was initialized as np.zeros before start of
                # conditional
                pass
            else:
                raise NotImplementedError('Only modelOutSideTemporalRange=="zero" implemented')


        # rectify the flux
        if self.rectifySED:
            flux = np.where(flambda>0., flambda, 0.)
        else:
            flux = flambda


        # convert per Ang to per nm
        flux *= 10.0
        # convert ang to nm
        wavelen = wavelen / 10.

        # If there is zero flux at 500nm, set
        # the flux in the slot closest to 500nm
        # equal to 0.01*minimum_non_zero_flux
        # (this is so SEDs used in PhoSim can have
        # finite normalization)
        if self.rectifySED:
            closest_to_500nm = np.argmin(np.abs(wavelen-500.0))
            if flux[closest_to_500nm] == 0.0:
                non_zero_flux = np.where(flux>0.0)
                if len(non_zero_flux[0])>0:
                    min_non_zero = np.min(flux[non_zero_flux])
                    flux[closest_to_500nm] = 0.01*min_non_zero

        sed = Sed(wavelen=wavelen, flambda=flux)
        # This has the cosmology built in.
        return sed

    def catsimBandFlux(self, time, bandpassobject):
        """
        return the flux in the bandpass in units of maggies which is the flux
        the AB magnitude reference spectrum would have in the same band.

        Parameters
        ----------
        time: mandatory, float
            MJD at which band fluxes are evaluated
        bandpassobject: mandatory, `lsst.sims.photUtils.BandPass` object
            A particular bandpass which is an instantiation of one of
            (u, g, r, i, z, y)
        Returns
        -------
        float value for flux in band in units of maggies

        Examples
        --------
        >>> bandpassnames = ['u', 'g', 'r', 'i', 'z', 'y']
        >>> LSST_BandPass = BandpassDict.loadTotalBandpassesFromFiles()
        >>> SN = SNObject(ra=30., dec=-60.)
        >>> SN.set(z=0.96, t0=571181, x1=2.66, c=0.353, x0=1.796112e-06)
        >>> SN.catsimBandFlux(bandpassobject=LSST_BandPass['r'], time=571190.)
        >>> 1.9856857972304903e-11

        .. note: If there is an unphysical value of sed in
        the wavelength range, it produces a flux of  `np.nan`
        """
        # Speedup for cases outside temporal range of model
        if time <= self.mintime() or time >= self.maxtime() :
            return 0.
        SEDfromSNcosmo = self.SNObjectSED(time=time,
                                          bandpass=bandpassobject)
        return SEDfromSNcosmo.calcFlux(bandpass=bandpassobject) / 3631.0
 
    def catsimBandMag(self, bandpassobject, time, fluxinMaggies=None,
                      noNan=False):
        """
        return the magnitude in the bandpass in the AB magnitude system

        Parameters
        ----------
        bandpassobject : mandatory,`sims.photUtils.BandPass` instances
            LSST Catsim bandpass instance for a particular bandpass
        time : mandatory, float
            MJD at which this is evaluated
        fluxinMaggies: float, defaults to None
            provide the flux in maggies, if not provided, this will be evaluated
        noNan : Bool, defaults to False
            If True, an AB magnitude of 200.0 rather than nan values is
            associated with a flux of 0.
        Returns
        -------
        float value of band magnitude in AB system

        Examples
        --------
        """
        if fluxinMaggies is None:
            fluxinMaggies = self.catsimBandFlux(bandpassobject=bandpassobject,
                                                time=time)
        if noNan:
            if fluxinMaggies <= 0.:
                return 200.0
        with np.errstate(divide='ignore', invalid='ignore'):
            return -2.5 * np.log10(fluxinMaggies)

    def catsimBandFluxError(self, time, bandpassobject, m5,
                            fluxinMaggies=None,
                            magnitude=None,
                            photParams=None):
        """
        return the flux uncertainty in the bandpass in units 'maggies'
        (the flux the AB magnitude reference spectrum would have in the
        same band.) for a source of given brightness. The source brightness
        may be calculated, but the need for calculation is overridden by a
        provided flux in bandpass (in units of maggies) which itself may be
        overridden by a provided magnitude. If the provided/calculated flux
        is 0. or negative the magnitude calculated is taken to be 200.0 rather
        than a np.nan.


        Parameters
        ----------
        time: mandatory, float
            MJD at which band fluxes are evaluated
        bandpassobject: mandatory, `lsst.sims.photUtils.BandPass` object
            A particular bandpass which is an instantiation of one of
            (u, g, r, i, z, y)
        m5 : float, mandatory
            fiveSigma Depth for the sky observation
        photParams : instance of `sims.photUtils.PhotometricParameters`, defaults to `None` 
            describes the hardware parameters of the Observing system
        magnitude : float, defaults to None
            AB magnitude of source in bandpass.
        fluxinMaggies : float, defaults to None
            flux in Maggies for source in bandpass
        Returns
        -------
        float

        Examples
        --------
        .. note: If there is an unphysical value of sed the fluxinMaggies might
        be `np.nan`. The magnitude calculated from this is calculated using `noNan`
         and is therefore 200.0 rather than `np.nan`. 
        """
        if fluxinMaggies is None:
            fluxinMaggies = self.catsimBandFlux(time=time,
                                                bandpassobject=bandpassobject)
        if magnitude is None:
            mag = self.catsimBandMag(time=time, fluxinMaggies=fluxinMaggies,
                                     bandpassobject=bandpassobject, noNan=True)
        else:
            mag = magnitude

        # recalculate fluxinMaggies as the previous one might have been `np.nan`
        # the noise is contaminated if this is `np.nan`
        fluxinMaggies = 10.0**(-0.4 * mag)

        if photParams is None:
            photParams = PhotometricParameters()

        SNR, gamma = calcSNR_m5(magnitude=mag, bandpass=bandpassobject,
                                m5=m5, photParams=photParams)
        return fluxinMaggies / SNR

    def catsimBandMagError(self, time, bandpassobject, m5, photParams=None,
                           magnitude=None):
        """
        return the 68 percent uncertainty on the magnitude in the bandpass

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
            mag = self.catsimBandMag(time=time,
                                     bandpassobject=bandpassobject,
                                     noNan=True)
        else:
            mag = magnitude

        bandpass = bandpassobject

        if photParams is None:
            photParams = PhotometricParameters()

        magerr = calcMagError_m5(magnitude=mag,
                                 bandpass=bandpassobject,
                                 m5=m5,
                                 photParams=photParams)
        return magerr[0]


    def catsimManyBandFluxes(self, time, bandpassDict,
                             observedBandPassInd=None):
        """
        return the flux in the multiple bandpasses of a bandpassDict
        indicated by observedBandPassInd in units of maggies

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

        with np.errstate(invalid='ignore', divide='ignore'):
            return -2.5 * np.log10(f)


    def catsimManyBandADUs(self, time, bandpassDict,
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

        bandpassNames = list(bandpassDict.keys())
        adus = np.zeros(len(bandpassNames))

        for i, filt in enumerate(bandpassNames):
            bandpass = bandpassDict[filt]
            adus[i] = SEDfromSNcosmo.calcADU(bandpass, photParams=photParams)

        return adus
