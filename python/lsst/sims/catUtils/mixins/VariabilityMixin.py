"""
This module defines methods which implement the astrophysical
variability models used by CatSim.  InstanceCatalogs apply
variability by calling the applyVariability() method in
the Variability class.  To add a new variability model to
this framework, users should define a method which accepts
as arguments:

valid_dexes -- the output of numpy.where() indicating
which astrophysical objects actually depend on the
variability model.

params -- a dict whose keys are the names of parameters
required by the variability model and whose values
are lists of the parameters required for the
variability model for all astrophysical objects
in the CatSim database (even those objects that do
not depend on the model; these objects can have
None in the parameter lists).

expmjd -- the MJD of the observation.

The method implementing the variability model should be
marked with the decorator @register_method(key) where key
is a string uniquely identifying the variability model.
applyVariability() will call the variability method
by reading in the json-ized dict varParamStr from the
CatSim database.  varParamStr should look like

{'m':method_name, 'p':{'p1': val1, 'p2': val2,...}}

method_name is the register_method() key referring
to the variabilty model. p1, p2, etc. are the parameters
expected by the variability model.
"""

from builtins import range
from builtins import object
import numpy
import linecache
import math
import os
import copy
import numbers
import json as json
from lsst.utils import getPackageDir
from lsst.sims.catalogs.decorators import register_method, compound
from lsst.sims.photUtils import Sed, BandpassDict
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d

__all__ = ["Variability", "VariabilityStars", "VariabilityGalaxies",
           "VariabilityAGN",
           "reset_agn_lc_cache", "StellarVariabilityModels",
           "ExtraGalacticVariabilityModels", "MLTflaringMixin"]

_AGN_LC_CACHE = {} # a global cache of agn light curve calculations
_MLT_LC_CACHE = None

def reset_agn_lc_cache():
    """
    Resets the _AGN_LC_CACHE (a global dict for cacheing time steps in AGN
    light curves) to an empty dict.
    """
    global _AGN_LC_CACHE
    _AGN_LC_CACHE = {}
    return None


class Variability(object):
    """
    Variability class for adding temporal variation to the magnitudes of
    objects in the base catalog.

    This class provides methods that all variability models rely on.
    Actual implementations of variability models will be provided by
    the *VariabilityModels classes.
    """

    _survey_start = 59580.0 # start time of the LSST survey being simulated (MJD)

    # the file wherein light curves for MLT dwarf flares are stored
    _mlt_lc_file = os.path.join(getPackageDir('sims_data'),
                                'catUtilsData', 'mdwarf_flare_light_curves_170412.npz')

    variabilityInitialized = False

    def num_variable_obj(self, params):
        """
        Return the total number of objects in the catalog

        Parameters
        ----------
        params is the dict of parameter arrays passed to a variability method

        Returns
        -------
        The number of objects in the catalog
        """
        params_keys = list(params.keys())
        if len(params_keys) == 0:
            return 0

        return len(params[params_keys[0]])

    def initializeVariability(self, doCache=False):
        """
        It will only be called from applyVariability, and only
        if self.variabilityInitiailized == True (which this method sets)

        @param [in] doCache controls whether or not the code caches calculated
        light curves for future use
        """
        # Docstring is a best approximation of what this method does.
        # This is older code.

        self.variabilityInitialized=True
        #below are variables to cache the light curves of variability models
        self.variabilityLcCache = {}
        self.variabilityCache = doCache
        try:
            self.variabilityDataDir = os.environ.get("SIMS_SED_LIBRARY_DIR")
        except:
            raise RuntimeError("sims_sed_library must be setup to compute variability because it contains"+
                               " the lightcurves")



    def applyVariability(self, varParams_arr):
        """
        Read in an array/list of varParamStr objects taken from the CatSim
        database.  For each varParamStr, call the appropriate variability
        model to calculate magnitude offsets that need to be applied to
        the corresponding astrophysical offsets.  Return a 2-D numpy
        array of magnitude offsets in which each row is an LSST band
        in ugrizy order and each column is an astrophysical object from
        the CatSim database.
        """

        # construct a registry of all of the variability models
        # available to the InstanceCatalog
        if not hasattr(self, '_methodRegistry'):
            self._methodRegistry = {}
            for methodname in dir(self):
                method=getattr(self, methodname)
                if hasattr(method, '_registryKey'):
                    if method._registryKey not in self._methodRegistry:
                        self._methodRegistry[method._registryKey] = method

        if self.variabilityInitialized == False:
            self.initializeVariability(doCache=True)

        # A numpy array of magnitude offsets.  Each row is
        # an LSST band in ugrizy order.  Each columns is an
        # astrophysical object from the CatSim database.
        deltaMag = numpy.zeros((6, len(varParams_arr)))

        # When the InstanceCatalog calls all of its getters
        # with an empty chunk to check column dependencies,
        # call all of the variability models in the
        # _methodRegistry to make sure that all of the column
        # dependencies of the variability models are detected.
        if len(varParams_arr) == 0:
            for method_name in self._methodRegistry:
                self._methodRegistry[method_name]([],{},0)

        # Keep a list of all of the specific variability models
        # that need to be called.  There is one entry for each
        # astrophysical object in the CatSim database.
        method_name_arr = []

        # Keep a dict keyed on all of the method names in
        # method_name_arr.  params[method_name] will be another
        # dict keyed on the names of the parameters required by
        # the method method_name.  The values of this dict will
        # be lists of parameter values for all astrophysical
        # objects in the CatSim database.  Even objects that
        # do no callon method_name will have entries in these
        # lists (they will be set to None).
        params = {}

        for ix, varCmd in enumerate(varParams_arr):
            if str(varCmd) != 'None':
                varCmd = json.loads(varCmd)

                # find the key associated with the name of
                # the specific variability model to be applied
                if 'varMethodName' in varCmd:
                    meth_key = 'varMethodName'
                else:
                    meth_key = 'm'

                # find the key associated with the list of
                # parameters to be supplied to the variability
                # model
                if 'pars' in varCmd:
                    par_key = 'pars'
                else:
                    par_key = 'p'
            else:
                # if there is no varParamStr, setup a null model
                varCmd = {'varMethodName': 'None', 'pars':{}}
                meth_key = 'varMethodName'
                par_key = 'pars'

            # if we have discovered a new variability model
            # that needs to be called, initialize its entries
            # in the params dict
            if varCmd[meth_key] not in method_name_arr:
                params[varCmd[meth_key]] = {}
                for p_name in varCmd[par_key]:
                    params[varCmd[meth_key]][p_name] = [None]*len(varParams_arr)

            method_name_arr.append(varCmd[meth_key])
            for p_name in varCmd[par_key]:
                params[varCmd[meth_key]][p_name][ix] = varCmd[par_key][p_name]

        method_name_arr = numpy.array(method_name_arr)
        for method_name in params:
            for p_name in params[method_name]:
                params[method_name][p_name] = numpy.array(params[method_name][p_name])

        expmjd=None

        # Loop over all of the variability models that need to be called.
        # Call each variability model on the astrophysical objects that
        # require the model.  Add the result to deltaMag.
        for method_name in numpy.unique(method_name_arr):
            if method_name != 'None':

                if expmjd is None:
                    expmjd = self.obs_metadata.mjd.TAI

                if method_name not in self._methodRegistry:
                    raise RuntimeError("Your InstanceCatalog does not contain " \
                                       + "a variability method corresponding to '%s'"
                                       % method_name)

                deltaMag += self._methodRegistry[method_name](numpy.where(numpy.char.equal(method_name, method_name_arr)),
                                                              params[method_name],
                                                              expmjd)

        return deltaMag


    def applyStdPeriodic(self, valid_dexes, params, keymap, expmjd,
                         inDays=True, interpFactory=None):

        """
        Applies a specified variability method.

        The params for the method are provided in the dict params{}

        The keys for those parameters are in the dict keymap{}

        This is because the syntax used here is not necessarily the syntax
        used in the data bases.

        The method will return a dict of magnitude offsets.  The dict will
        be keyed to the filter names.

        @param [in] valid_dexes is the result of numpy.where() indicating
        which astrophysical objects from the CatSim database actually use
        this variability model.

        @param [in] params is a dict of parameters for the variability model.
        The dict is keyed to the names of parameters.  The values are arrays
        of parameter values.

        @param [in] keymap is a dict mapping from the parameter naming convention
        used by the database to the parameter naming convention used by the
        variability methods below.

        @param [in] expmjd is the mjd of the observation

        @param [in] inDays controls whether or not the time grid
        of the light curve is renormalized by the period

        @param [in] interpFactory is the method used for interpolating
        the light curve

        @param [out] magoff is a 2D numpy array of magnitude offsets.  Each
        row is an LSST band in ugrizy order.  Each column is a different
        astrophysical object from the CatSim database.
        """
        magoff = numpy.zeros((6, self.num_variable_obj(params)))
        expmjd = numpy.asarray(expmjd)
        for ix in valid_dexes[0]:
            filename = params[keymap['filename']][ix]
            toff = params[keymap['t0']][ix]

            inPeriod = None
            if 'period' in params:
                inPeriod = params['period'][ix]

            epoch = expmjd - toff
            if filename in self.variabilityLcCache:
                splines = self.variabilityLcCache[filename]['splines']
                period = self.variabilityLcCache[filename]['period']
            else:
                lc = numpy.loadtxt(os.path.join(self.variabilityDataDir,filename), unpack=True, comments='#')
                if inPeriod is None:
                    dt = lc[0][1] - lc[0][0]
                    period = lc[0][-1] + dt
                else:
                    period = inPeriod

                if inDays:
                    lc[0] /= period

                splines  = {}

                if interpFactory is not None:
                    splines['u'] = interpFactory(lc[0], lc[1])
                    splines['g'] = interpFactory(lc[0], lc[2])
                    splines['r'] = interpFactory(lc[0], lc[3])
                    splines['i'] = interpFactory(lc[0], lc[4])
                    splines['z'] = interpFactory(lc[0], lc[5])
                    splines['y'] = interpFactory(lc[0], lc[6])
                    if self.variabilityCache:
                        self.variabilityLcCache[filename] = {'splines':splines, 'period':period}
                else:
                    splines['u'] = interp1d(lc[0], lc[1])
                    splines['g'] = interp1d(lc[0], lc[2])
                    splines['r'] = interp1d(lc[0], lc[3])
                    splines['i'] = interp1d(lc[0], lc[4])
                    splines['z'] = interp1d(lc[0], lc[5])
                    splines['y'] = interp1d(lc[0], lc[6])
                    if self.variabilityCache:
                        self.variabilityLcCache[filename] = {'splines':splines, 'period':period}

            phase = epoch/period - epoch//period
            magoff[0][ix] = splines['u'](phase)
            magoff[1][ix] = splines['g'](phase)
            magoff[2][ix] = splines['r'](phase)
            magoff[3][ix] = splines['i'](phase)
            magoff[4][ix] = splines['z'](phase)
            magoff[5][ix] = splines['y'](phase)

        return magoff


class StellarVariabilityModels(Variability):
    """
    A mixin providing standard stellar variability models.
    """

    @register_method('applyRRly')
    def applyRRly(self, valid_dexes, params, expmjd):

        if len(params) == 0:
            return numpy.array([[],[],[],[],[],[]])

        keymap = {'filename':'filename', 't0':'tStartMjd'}
        return self.applyStdPeriodic(valid_dexes, params, keymap, expmjd,
                interpFactory=InterpolatedUnivariateSpline)

    @register_method('applyCepheid')
    def applyCepheid(self, valid_dexes, params, expmjd):

        if len(params) == 0:
            return numpy.array([[],[],[],[],[],[]])

        keymap = {'filename':'lcfile', 't0':'t0'}
        return self.applyStdPeriodic(valid_dexes, params, keymap, expmjd, inDays=False,
                interpFactory=InterpolatedUnivariateSpline)

    @register_method('applyEb')
    def applyEb(self, valid_dexes, params, expmjd):

        if len(params) == 0:
            return numpy.array([[],[],[],[],[],[]])

        keymap = {'filename':'lcfile', 't0':'t0'}
        dMags = self.applyStdPeriodic(valid_dexes, params, keymap, expmjd,
                                      inDays=False,
                                      interpFactory=InterpolatedUnivariateSpline)
        if len(dMags)>0:
            if dMags.min()<0.0:
                raise RuntimeError("Negative delta flux in applyEb")
        dMags = -2.5*numpy.log10(dMags)
        return dMags

    @register_method('applyMicrolensing')
    def applyMicrolensing(self, valid_dexes, params, expmjd_in):
        return self.applyMicrolens(valid_dexes, params,expmjd_in)

    @register_method('applyMicrolens')
    def applyMicrolens(self, valid_dexes, params, expmjd_in):
        #I believe this is the correct method based on
        #http://www.physics.fsu.edu/Courses/spring98/AST3033/Micro/lensing.htm
        #
        #21 October 2014
        #This method assumes that the parameters for microlensing variability
        #are stored in a varParamStr column in the database.  Actually, the
        #current microlensing event tables in the database store each
        #variability parameter as its own database column.
        #At some point, either this method or the microlensing tables in the
        #database will need to be changed.

        if len(params) == 0:
            return numpy.array([[],[],[],[],[],[]])

        dMags = numpy.zeros((6, self.num_variable_obj(params)))
        expmjd = numpy.asarray(expmjd_in,dtype=float)
        epochs = expmjd - params['t0'][valid_dexes].astype(float)
        umin = params['umin'].astype(float)[valid_dexes]
        that = params['that'].astype(float)[valid_dexes]
        u = numpy.sqrt(umin**2 + ((2.0*epochs/that)**2))
        magnification = (u**2+2.0)/(u*numpy.sqrt(u**2+4.0))
        dmag = -2.5*numpy.log10(magnification)
        for ix in range(6):
            dMags[ix][valid_dexes] += dmag
        return dMags


    @register_method('applyAmcvn')
    def applyAmcvn(self, valid_dexes, params, expmjd_in):
        #21 October 2014
        #This method assumes that the parameters for Amcvn variability
        #are stored in a varParamStr column in the database.  Actually, the
        #current Amcvn event tables in the database store each
        #variability parameter as its own database column.
        #At some point, either this method or the Amcvn tables in the
        #database will need to be changed.

        if not isinstance(expmjd_in, numbers.Number):
            raise RuntimeError("Cannot pass multiple expMJD into applyAmcvn")

        if len(params) == 0:
            return numpy.array([[],[],[],[],[],[]])

        maxyears = 10.
        dMag = numpy.zeros((6, self.num_variable_obj(params)))
        epoch = expmjd_in

        amplitude = params['amplitude'].astype(float)[valid_dexes]
        t0 = params['t0'].astype(float)[valid_dexes]
        period = params['period'].astype(float)[valid_dexes]
        burst_freq = params['burst_freq'].astype(float)[valid_dexes]
        burst_scale = params['burst_scale'].astype(float)[valid_dexes]
        amp_burst = params['amp_burst'].astype(float)[valid_dexes]
        color_excess = params['color_excess_during_burst'].astype(float)[valid_dexes]
        does_burst = params['does_burst'][valid_dexes]

        # get the light curve of the typical variability
        uLc   = amplitude*numpy.cos((epoch - t0)/period)
        gLc   = copy.deepcopy(uLc)
        rLc   = copy.deepcopy(uLc)
        iLc   = copy.deepcopy(uLc)
        zLc   = copy.deepcopy(uLc)
        yLc   = copy.deepcopy(uLc)

        # add in the flux from any bursting
        local_bursting_dexes = numpy.where(does_burst==1)
        for i_burst in local_bursting_dexes[0]:
            adds = 0.0
            for o in numpy.linspace(t0[i_burst] + burst_freq[i_burst],\
                                 t0[i_burst] + maxyears*365.25, \
                                 numpy.ceil(maxyears*365.25/burst_freq[i_burst])):
                tmp = numpy.exp( -1*(epoch - o)/burst_scale[i_burst])/numpy.exp(-1.)
                adds -= amp_burst[i_burst]*tmp*(tmp < 1.0)  ## kill the contribution
            ## add some blue excess during the outburst
            uLc[i_burst] += adds +  2.0*color_excess[i_burst]
            gLc[i_burst] += adds + color_excess[i_burst]
            rLc[i_burst] += adds + 0.5*color_excess[i_burst]
            iLc[i_burst] += adds
            zLc[i_burst] += adds
            yLc[i_burst] += adds

        dMag[0][valid_dexes] += uLc
        dMag[1][valid_dexes] += gLc
        dMag[2][valid_dexes] += rLc
        dMag[3][valid_dexes] += iLc
        dMag[4][valid_dexes] += zLc
        dMag[5][valid_dexes] += yLc
        return dMag

    @register_method('applyBHMicrolens')
    def applyBHMicrolens(self, valid_dexes, params, expmjd_in):
        #21 October 2014
        #This method assumes that the parameters for BHMicrolensing variability
        #are stored in a varParamStr column in the database.  Actually, the
        #current BHMicrolensing event tables in the database store each
        #variability parameter as its own database column.
        #At some point, either this method or the BHMicrolensing tables in the
        #database will need to be changed.

        if len(params) == 0:
            return numpy.array([[],[],[],[],[],[]])

        magoff = numpy.zeros((6, self.num_variable_obj(params)))
        expmjd = numpy.asarray(expmjd_in,dtype=float)
        filename_arr = params['filename']
        toff_arr = params['t0'].astype(float)
        for ix in valid_dexes[0]:
            toff = toff_arr[ix]
            filename = filename_arr[ix]
            epoch = expmjd - toff
            lc = numpy.loadtxt(os.path.join(self.variabilityDataDir, filename), unpack=True, comments='#')
            dt = lc[0][1] - lc[0][0]
            period = lc[0][-1]
            #BH lightcurves are in years
            lc[0] *= 365.
            minage = lc[0][0]
            maxage = lc[0][-1]
            #I'm assuming that these are all single point sources lensed by a
            #black hole.  These also can be used to simulate binary systems.
            #Should be 8kpc away at least.
            magnification = InterpolatedUnivariateSpline(lc[0], lc[1])
            moff = -2.5*numpy.log(magnification(epoch))
            for ii in range(6):
                magoff[ii][ix] = moff

        return magoff


class MLTflaringMixin(Variability):
    """
    A mixin providing the model for cool dwarf stellar flares.
    """

    @register_method('MLT')
    def applyMlTflaring(self, valid_dexes, params, expmjd):

        parallax = self.column_by_name('parallax')
        ebv = self.column_by_name('ebv')

        global _MLT_LC_CACHE

        # this needs to occur before loading the MLT light curve cache,
        # just in case the user wants to override the light curve cache
        # file by hand before generating the catalog
        if len(params) == 0:
            return numpy.array([[],[],[],[],[],[]])

        if not hasattr(self, 'photParams'):
            raise RuntimeError("To apply MLT dwarf flaring, your "
                               "InstanceCatalog must have a member variable "
                               "photParams which is an instantiation of the "
                               "class PhotometricParameters, which can be "
                               "imported from lsst.sims.photUtils. "
                               "This is so that your InstanceCatalog has "
                               "knowledge of the effective area of the LSST "
                               "mirror.")

        if _MLT_LC_CACHE is None:

            if not os.path.exists(self._mlt_lc_file):
                catutils_scripts = os.path.join(getPackageDir('sims_catUtils'), 'support_scripts')
                raise RuntimeError("The MLT flaring light curve file:\n"
                                    + "\n%s\n" % self._mlt_lc_file
                                    + "\ndoes not exist."
                                    +"\n\n"
                                    + "Go into %s " % catutils_scripts
                                    + "and run get_mdwarf_flares.sh "
                                    + "to get the data")

            _MLT_LC_CACHE = numpy.load(self._mlt_lc_file)

        if not hasattr(self, '_mlt_dust_lookup'):
            # Construct a look-up table to determine the factor
            # by which to multiply the flares' flux to account for
            # dust as a function of E(B-V).  Recall that we are
            # modeling all MLT flares as 9000K blackbodies.

            if not hasattr(self, 'lsstBandpassDict'):
                raise RuntimeError('You are asking for MLT dwarf flaring '
                                   'magnitudes in a catalog that has not '
                                   'defined lsstBandpassDict.  The MLT '
                                   'flaring magnitudes model does not know '
                                   'how to apply dust extinction to the '
                                   'flares without the member variable '
                                   'lsstBandpassDict being defined.')

            ebv_grid = numpy.arange(0.0, 7.01, 0.01)
            bb_wavelen = numpy.arange(200.0, 1500.0, 0.1)
            hc_over_k = 1.4387e7  # nm*K
            temp = 9000.0  # black body temperature in Kelvin
            exp_arg = hc_over_k/(temp*bb_wavelen)
            exp_term = 1.0/(numpy.exp(exp_arg) - 1.0)
            ln_exp_term = numpy.log(exp_term)

            log_bb_flambda = -5.0*numpy.log(bb_wavelen) + ln_exp_term
            bb_flambda = numpy.exp(log_bb_flambda)
            bb_sed = Sed(wavelen=bb_wavelen, flambda=bb_flambda)

            base_fluxes = self.lsstBandpassDict.fluxListForSed(bb_sed)

            a_x, b_x = bb_sed.setupCCMab()
            self._mlt_dust_lookup = {}
            self._mlt_dust_lookup['ebv'] = ebv_grid
            list_of_bp = self.lsstBandpassDict.keys()
            for bp in list_of_bp:
                self._mlt_dust_lookup[bp] = numpy.zeros(len(ebv_grid))
            for iebv, ebv_val in enumerate(ebv_grid):
                wv, fl = bb_sed.addCCMDust(a_x, b_x,
                                           ebv=ebv_val,
                                           wavelen=bb_wavelen,
                                           flambda=bb_flambda)

                dusty_bb = Sed(wavelen=wv, flambda=fl)
                dusty_fluxes = self.lsstBandpassDict.fluxListForSed(dusty_bb)
                for ibp, bp in enumerate(list_of_bp):
                    self._mlt_dust_lookup[bp][iebv] = dusty_fluxes[ibp]/base_fluxes[ibp]

        # get the distance to each star in parsecs
        _au_to_parsec = 1.0/206265.0
        dd = _au_to_parsec/parallax

        # get the area of the sphere through which the star's energy
        # is radiating to get to us (in cm^2)
        _cm_per_parsec = 3.08576e16
        sphere_area = 4.0*numpy.pi*numpy.power(dd*_cm_per_parsec, 2)

        flux_factor = self.photParams.effarea/sphere_area

        dMags = numpy.zeros((6, self.num_variable_obj(params)))
        mag_name_tuple = ('u', 'g', 'r', 'i', 'z', 'y')
        base_fluxes = {}
        base_mags = {}
        ss = Sed()
        for mag_name in mag_name_tuple:
            if ('lsst_%s' % mag_name in self._actually_calculated_columns or
                'delta_lsst_%s' % mag_name in self._actually_calculated_columns):

                mm = self.column_by_name('quiescent_lsst_%s' % mag_name)
                base_mags[mag_name] = mm
                base_fluxes[mag_name] = ss.fluxFromMag(mm)

        lc_names_unique = numpy.unique(params['lc'].astype(str))
        for lc_name in lc_names_unique:
            use_this_lc = numpy.where(numpy.char.find(params['lc'], lc_name)==0)

            lc_name = lc_name.replace('.txt', '')
            if 'late' in lc_name:
                lc_name = lc_name.replace('in', '')
            time_arr = self._survey_start + _MLT_LC_CACHE['%s_time' % lc_name]
            dt = time_arr.max() - time_arr.min()

            t_interp = self.obs_metadata.mjd.TAI + params['t0'][use_this_lc]
            while t_interp.max() > time_arr.max():
                bad_dexes = numpy.where(t_interp>time_arr.max())
                t_interp[bad_dexes] -= dt

            for i_mag, mag_name in enumerate(mag_name_tuple):
                if ('lsst_%s' % mag_name in self._actually_calculated_columns or
                    'delta_lsst_%s' % mag_name in self._actually_calculated_columns):
                    flux_arr = _MLT_LC_CACHE['%s_%s' % (lc_name, mag_name)]
                    dflux = numpy.interp(t_interp, time_arr, flux_arr)
                    dflux *= flux_factor[use_this_lc]

                    dust_factor = numpy.interp(ebv[use_this_lc],
                                               self._mlt_dust_lookup['ebv'],
                                               self._mlt_dust_lookup[mag_name])
                    dflux *= dust_factor

                    dMags[i_mag][use_this_lc] = (ss.magFromFlux(base_fluxes[mag_name][use_this_lc] + dflux)
                                                 - base_mags[mag_name][use_this_lc])

        return dMags


class ExtraGalacticVariabilityModels(Variability):
    """
    A mixin providing the model for AGN variability.
    """

    @register_method('applyAgn')
    def applyAgn(self, valid_dexes, params, expmjd):

        global _AGN_LC_CACHE

        if len(params) == 0:
            return numpy.array([[],[],[],[],[],[]])

        dMags = numpy.zeros((6, self.num_variable_obj(params)))
        toff_arr = params['t0_mjd'].astype(float)
        seed_arr = params['seed']
        tau_arr = params['agn_tau'].astype(float)
        sfu_arr = params['agn_sfu'].astype(float)
        sfg_arr = params['agn_sfg'].astype(float)
        sfr_arr = params['agn_sfr'].astype(float)
        sfi_arr = params['agn_sfi'].astype(float)
        sfz_arr = params['agn_sfz'].astype(float)
        sfy_arr = params['agn_sfy'].astype(float)

        for ix in valid_dexes[0]:
            toff = toff_arr[ix]
            seed = seed_arr[ix]
            tau = tau_arr[ix]

            sfint = {}
            sfint['u'] = sfu_arr[ix]
            sfint['g'] = sfg_arr[ix]
            sfint['r'] = sfr_arr[ix]
            sfint['i'] = sfi_arr[ix]
            sfint['z'] = sfz_arr[ix]
            sfint['y'] = sfy_arr[ix]

            # A string made up of this AGNs variability parameters that ought
            # to uniquely identify it.
            #
            agn_ID = '%d_%.12f_%.12f_%.12f_%.12f_%.12f_%.12f_%.12f_%.12f' \
            %(seed, sfint['u'], sfint['g'], sfint['r'], sfint['i'], sfint['z'],
              sfint['y'], tau, toff)

            resumption = False

            # Check to see if this AGN has already been simulated.
            # If it has, see if the previously simulated MJD is
            # earlier than the first requested MJD.  If so,
            # use that previous simulation as the starting point.
            #
            if agn_ID in _AGN_LC_CACHE:
                if _AGN_LC_CACHE[agn_ID]['mjd'] <expmjd:
                    resumption = True

            if resumption:
                rng = copy.deepcopy(_AGN_LC_CACHE[agn_ID]['rng'])
                start_date = _AGN_LC_CACHE[agn_ID]['mjd']
                dx_0 = _AGN_LC_CACHE[agn_ID]['dx']
            else:
                start_date = toff
                rng = numpy.random.RandomState(seed)
                dx_0 = {}
                for k in sfint:
                    dx_0[k]=0.0

            endepoch = expmjd - start_date

            if endepoch < 0:
                raise RuntimeError("WARNING: Time offset greater than minimum epoch.  " +
                                   "Not applying variability. "+
                                   "expmjd: %e should be > toff: %e  " % (expmjd, toff) +
                                   "in applyAgn variability method")

            dt = tau/100.
            nbins = int(math.ceil(endepoch/dt))

            x1 = (nbins-1)*dt
            x2 = (nbins)*dt

            dt = dt/tau
            es = rng.normal(0., 1., nbins)*math.sqrt(dt)
            dx_cached = {}

            for k, ik in zip(('u', 'g', 'r', 'i', 'z', 'y'), range(6)):
                dx2 = dx_0[k]
                for i in range(nbins):
                    #The second term differs from Zeljko's equation by sqrt(2.)
                    #because he assumes stdev = sfint/sqrt(2)
                    dx1 = dx2
                    dx2 = -dx1*dt + sfint[k]*es[i] + dx1

                dx_cached[k] = dx2
                dMags[ik][ix] = (endepoch*(dx1-dx2)+dx2*x1-dx1*x2)/(x1-x2)

            # Reset that AGN light curve cache once it contains
            # one million objects (to prevent it from taking up
            # too much memory).
            if len(_AGN_LC_CACHE)>1000000:
                reset_agn_lc_cache()

            if agn_ID not in _AGN_LC_CACHE:
                _AGN_LC_CACHE[agn_ID] = {}

            _AGN_LC_CACHE[agn_ID]['mjd'] = start_date+x2
            _AGN_LC_CACHE[agn_ID]['rng'] = copy.deepcopy(rng)
            _AGN_LC_CACHE[agn_ID]['dx'] = dx_cached

        return dMags


class _VariabilityPointSources(object):

    @compound('delta_lsst_u', 'delta_lsst_g', 'delta_lsst_r',
             'delta_lsst_i', 'delta_lsst_z', 'delta_lsst_y')
    def get_stellar_variability(self):
        """
        Getter for the change in magnitudes due to stellar
        variability.  The PhotometryStars mixin is clever enough
        to automatically add this to the baseline magnitude.
        """

        varParams = self.column_by_name('varParamStr')
        return self.applyVariability(varParams)


class VariabilityStars(_VariabilityPointSources, StellarVariabilityModels,
                       MLTflaringMixin):
    """
    This is a mixin which wraps the methods from the class
    StellarVariabilityModels into getters for InstanceCatalogs
    (specifically, InstanceCatalogs of stars).  Getters in
    this method should define columns named like

    delta_columnName

    where columnName is the name of the baseline (non-varying) magnitude
    column to which delta_columnName will be added.  The getters in the
    photometry mixins will know to find these columns and add them to
    columnName, provided that the columns here follow this naming convention.

    Thus: merely including VariabilityStars in the inheritance tree of
    an InstanceCatalog daughter class will activate variability for any column
    for which delta_columnName is defined.
    """
    pass


class VariabilityAGN(_VariabilityPointSources, ExtraGalacticVariabilityModels):
    """
    This is a mixin which wraps the methods from the class
    ExtraGalacticVariabilityModels into getters for InstanceCatalogs
    of AGN.  Getters in this method should define columns named like

    delta_columnName

    where columnName is the name of the baseline (non-varying) magnitude
    column to which delta_columnName will be added.  The getters in the
    photometry mixins will know to find these columns and add them to
    columnName, provided that the columns here follow this naming convention.

    Thus: merely including VariabilityStars in the inheritance tree of
    an InstanceCatalog daughter class will activate variability for any column
    for which delta_columnName is defined.
    """
    pass


class VariabilityGalaxies(ExtraGalacticVariabilityModels):
    """
    This is a mixin which wraps the methods from the class
    ExtraGalacticVariabilityModels into getters for InstanceCatalogs
    (specifically, InstanceCatalogs of galaxies).  Getters in this
    method should define columns named like

    delta_columnName

    where columnName is the name of the baseline (non-varying) magnitude
    column to which delta_columnName will be added.  The getters in the
    photometry mixins will know to find these columns and add them to
    columnName, provided that the columns here follow this naming convention.

    Thus: merely including VariabilityStars in the inheritance tree of
    an InstanceCatalog daughter class will activate variability for any column
    for which delta_columnName is defined.
    """

    @compound('delta_uAgn', 'delta_gAgn', 'delta_rAgn',
              'delta_iAgn', 'delta_zAgn', 'delta_yAgn')
    def get_galaxy_variability_total(self):

        """
        Getter for the change in magnitude due to AGN
        variability.  The PhotometryGalaxies mixin is
        clever enough to automatically add this to
        the baseline magnitude.
        """
        varParams = self.column_by_name("varParamStr")
        return self.applyVariability(varParams)
