"""
This module defines methods which implement the astrophysical
variability models used by CatSim.  InstanceCatalogs apply
variability by calling the applyVariability() method in
the Variability class.  To add a new variability model to
this framework, users should define a method which
returns the delta magnitude of the variable source
in the LSST bands and accepts as arguments:

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

    expmjd -- the MJD of the observation.  This must be
    able to accept a float or a numpy array.

If expmjd is a float, the variability model should
return a 2-D numpy array in which the first index
varies over the band and the second index varies
over the object, i.e.

    out_dmag[0][2] is the delta magnitude of the 2nd object
    in the u band

    out_dmag[3][15] is the delta magnitude of the 15th object
    in the i band.

If expmjd is a numpy array, the variability should
return a 3-D numpy array in which the first index
varies over the band, the second index varies over
the object, and the third index varies over the
time step, i.e.

    out_dmag[0][2][15] is the delta magnitude of the 2nd
    object in the u band at the 15th value of expmjd

    out_dmag[3][11][2] is the delta magnitude of the
    11th object in the i band at the 2nd value of
    expmjd

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
import numpy as np
import linecache
import math
import os
import gzip
import copy
import numbers
import json as json
from lsst.utils import getPackageDir
from lsst.sims.catalogs.decorators import register_method, compound
from lsst.sims.photUtils import Sed, BandpassDict
from lsst.sims.utils.CodeUtilities import sims_clean_up
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d

import time

__all__ = ["Variability", "VariabilityStars", "VariabilityGalaxies",
           "VariabilityAGN",
           "reset_agn_lc_cache", "StellarVariabilityModels",
           "ExtraGalacticVariabilityModels", "MLTflaringMixin",
           "ParametrizedLightCurveMixin"]

_AGN_LC_CACHE = {}  # a global cache of agn light curve calculations

_MLT_LC_NPZ = None  # this will be loaded from a .npz file
                    # (.npz files are the result of numpy.savez())

_MLT_LC_NPZ_NAME = None  # the name of the .npz file to beloaded

_MLT_LC_TIME_CACHE = {}  # a dict for storing loaded time grids

_MLT_LC_FLUX_CACHE = {}  # a dict for storing loaded flux grids

_PARAMETRIZED_LC_MODELS = {}  # a dict for storing the parametrized light curve models
_PARAMETRIZED_MODELS_LOADED = []  # a list of all of the files from which models were loaded

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
        if self.variabilityInitiailized == False (which this method then
        sets to True)

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



    def applyVariability(self, varParams_arr, expmjd=None):
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


        if isinstance(expmjd, numbers.Number) or expmjd is None:
            # A numpy array of magnitude offsets.  Each row is
            # an LSST band in ugrizy order.  Each column is an
            # astrophysical object from the CatSim database.
            deltaMag = np.zeros((6, len(varParams_arr)))
        else:
            # the last dimension varies over time
            deltaMag = np.zeros((6, len(varParams_arr), len(expmjd)))

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
        # astrophysical object in the CatSim database.  We will
        # ultimately run np.where on method_name_arr to determine
        # which objects need to be passed through which
        # variability methods.
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

        method_name_arr = np.array(method_name_arr)
        for method_name in params:
            for p_name in params[method_name]:
                params[method_name][p_name] = np.array(params[method_name][p_name])

        # Loop over all of the variability models that need to be called.
        # Call each variability model on the astrophysical objects that
        # require the model.  Add the result to deltaMag.
        for method_name in np.unique(method_name_arr):
            if method_name != 'None':

                if expmjd is None:
                    expmjd = self.obs_metadata.mjd.TAI

                if method_name not in self._methodRegistry:
                    raise RuntimeError("Your InstanceCatalog does not contain " \
                                       + "a variability method corresponding to '%s'"
                                       % method_name)

                deltaMag += self._methodRegistry[method_name](np.where(np.char.equal(method_name, method_name_arr)),
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
        if isinstance(expmjd, numbers.Number):
            magoff = np.zeros((6, self.num_variable_obj(params)))
        else:
            magoff = np.zeros((6, self.num_variable_obj(params), len(expmjd)))
        expmjd = np.asarray(expmjd)
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
                lc = np.loadtxt(os.path.join(self.variabilityDataDir,filename), unpack=True, comments='#')
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
            return np.array([[],[],[],[],[],[]])

        keymap = {'filename':'filename', 't0':'tStartMjd'}
        return self.applyStdPeriodic(valid_dexes, params, keymap, expmjd,
                interpFactory=InterpolatedUnivariateSpline)

    @register_method('applyCepheid')
    def applyCepheid(self, valid_dexes, params, expmjd):

        if len(params) == 0:
            return np.array([[],[],[],[],[],[]])

        keymap = {'filename':'lcfile', 't0':'t0'}
        return self.applyStdPeriodic(valid_dexes, params, keymap, expmjd, inDays=False,
                interpFactory=InterpolatedUnivariateSpline)

    @register_method('applyEb')
    def applyEb(self, valid_dexes, params, expmjd):

        if len(params) == 0:
            return np.array([[],[],[],[],[],[]])

        keymap = {'filename':'lcfile', 't0':'t0'}
        d_fluxes = self.applyStdPeriodic(valid_dexes, params, keymap, expmjd,
                                         inDays=False,
                                          interpFactory=InterpolatedUnivariateSpline)
        if len(d_fluxes)>0:
            if d_fluxes.min()<0.0:
                raise RuntimeError("Negative delta flux in applyEb")
        if isinstance(expmjd, numbers.Number):
            dMags = np.zeros((6, self.num_variable_obj(params)))
        else:
            dMags = np.zeros((6, self.num_variable_obj(params), len(expmjd)))
        dmag_vals = -2.5*np.log10(d_fluxes)
        dMags += np.where(np.logical_not(np.logical_or(np.isnan(dmag_vals), np.isinf(dmag_vals))),
                          dmag_vals, 0.0)
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
            return np.array([[],[],[],[],[],[]])

        expmjd = np.asarray(expmjd_in,dtype=float)
        if isinstance(expmjd_in, numbers.Number):
            dMags = np.zeros((6, self.num_variable_obj(params)))
            epochs = expmjd - params['t0'][valid_dexes].astype(float)
            umin = params['umin'].astype(float)[valid_dexes]
            that = params['that'].astype(float)[valid_dexes]
        else:
            dMags = np.zeros((6, self.num_variable_obj(params), len(expmjd)))
            # cast epochs, umin, that into 2-D numpy arrays; the first index will iterate
            # over objects; the second index will iterate over times in expmjd
            epochs = np.array([expmjd - t0 for t0 in params['t0'][valid_dexes].astype(float)])
            umin = np.array([[uu]*len(expmjd) for uu in params['umin'].astype(float)[valid_dexes]])
            that = np.array([[tt]*len(expmjd) for tt in params['that'].astype(float)[valid_dexes]])

        u = np.sqrt(umin**2 + ((2.0*epochs/that)**2))
        magnification = (u**2+2.0)/(u*np.sqrt(u**2+4.0))
        dmag = -2.5*np.log10(magnification)
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

        if len(params) == 0:
            return np.array([[],[],[],[],[],[]])

        maxyears = 10.
        if isinstance(expmjd_in, numbers.Number):
            dMag = np.zeros((6, self.num_variable_obj(params)))
            amplitude = params['amplitude'].astype(float)[valid_dexes]
            t0_arr = params['t0'].astype(float)[valid_dexes]
            period = params['period'].astype(float)[valid_dexes]
            epoch_arr = expmjd_in
        else:
            dMag = np.zeros((6, self.num_variable_obj(params), len(expmjd_in)))
            n_time = len(expmjd_in)
            t0_arr = np.array([[tt]*n_time for tt in params['t0'].astype(float)[valid_dexes]])
            amplitude = np.array([[aa]*n_time for aa in params['amplitude'].astype(float)[valid_dexes]])
            period = np.array([[pp]*n_time for pp in params['period'].astype(float)[valid_dexes]])
            epoch_arr = np.array([expmjd_in]*len(valid_dexes[0]))

        epoch = expmjd_in

        t0 = params['t0'].astype(float)[valid_dexes]
        burst_freq = params['burst_freq'].astype(float)[valid_dexes]
        burst_scale = params['burst_scale'].astype(float)[valid_dexes]
        amp_burst = params['amp_burst'].astype(float)[valid_dexes]
        color_excess = params['color_excess_during_burst'].astype(float)[valid_dexes]
        does_burst = params['does_burst'][valid_dexes]

        # get the light curve of the typical variability
        uLc   = amplitude*np.cos((epoch_arr - t0_arr)/period)
        gLc   = copy.deepcopy(uLc)
        rLc   = copy.deepcopy(uLc)
        iLc   = copy.deepcopy(uLc)
        zLc   = copy.deepcopy(uLc)
        yLc   = copy.deepcopy(uLc)

        # add in the flux from any bursting
        local_bursting_dexes = np.where(does_burst==1)
        for i_burst in local_bursting_dexes[0]:
            adds = 0.0
            for o in np.linspace(t0[i_burst] + burst_freq[i_burst],\
                                 t0[i_burst] + maxyears*365.25, \
                                 np.ceil(maxyears*365.25/burst_freq[i_burst])):
                tmp = np.exp( -1*(epoch - o)/burst_scale[i_burst])/np.exp(-1.)
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
            return np.array([[],[],[],[],[],[]])

        if isinstance(expmjd_in, numbers.Number):
            magoff = np.zeros((6, self.num_variable_obj(params)))
        else:
            magoff = np.zeros((6, self.num_variable_obj(params), len(expmjd_in)))
        expmjd = np.asarray(expmjd_in,dtype=float)
        filename_arr = params['filename']
        toff_arr = params['t0'].astype(float)
        for ix in valid_dexes[0]:
            toff = toff_arr[ix]
            filename = filename_arr[ix]
            epoch = expmjd - toff
            lc = np.loadtxt(os.path.join(self.variabilityDataDir, filename), unpack=True, comments='#')
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
            mag_val = magnification(epoch)
            # If we are interpolating out of the light curve's domain, set
            # the magnification equal to 1
            mag_val = np.where(np.isnan(mag_val), 1.0, mag_val)
            moff = -2.5*np.log(mag_val)
            for ii in range(6):
                magoff[ii][ix] = moff

        return magoff


class MLTflaringMixin(Variability):
    """
    A mixin providing the model for cool dwarf stellar flares.
    """

    # the file wherein light curves for MLT dwarf flares are stored
    _mlt_lc_file = os.path.join(getPackageDir('sims_data'),
                                'catUtilsData', 'mdwarf_flare_light_curves_170412.npz')

    @register_method('MLT')
    def applyMLTflaring(self, valid_dexes, params, expmjd,
                        parallax=None, ebv=None, quiescent_mags=None):
        """
        parallax, ebv, and quiescent_mags are optional kwargs for use if you are
        calling this method outside the context of an InstanceCatalog (presumably
        with a numpy array of expmjd)

        parallax is the parallax of your objects in radians

        ebv is the E(B-V) value for your objects

        quiescent_mags is a dict keyed on ('u', 'g', 'r', 'i', 'z', 'y')
        with the quiescent magnitudes of the objects
        """

        t_start = time.time()
        t_before = t_start

        if parallax is None:
            parallax = self.column_by_name('parallax')
        if ebv is None:
            ebv = self.column_by_name('ebv')

        global _MLT_LC_NPZ
        global _MLT_LC_NPZ_NAME
        global _MLT_LC_TIME_CACHE
        global _MLT_LC_FLUX_CACHE

        # this needs to occur before loading the MLT light curve cache,
        # just in case the user wants to override the light curve cache
        # file by hand before generating the catalog
        if len(params) == 0:
            return np.array([[],[],[],[],[],[]])

        if quiescent_mags is None:
            quiescent_mags = {}
            for mag_name in ('u', 'g', 'r', 'i', 'z', 'y'):
                if ('lsst_%s' % mag_name in self._actually_calculated_columns or
                    'delta_lsst_%s' % mag_name in self._actually_calculated_columns):

                    quiescent_mags[mag_name] = self.column_by_name('quiescent_lsst_%s' % mag_name)

        if not hasattr(self, 'photParams'):
            raise RuntimeError("To apply MLT dwarf flaring, your "
                               "InstanceCatalog must have a member variable "
                               "photParams which is an instantiation of the "
                               "class PhotometricParameters, which can be "
                               "imported from lsst.sims.photUtils. "
                               "This is so that your InstanceCatalog has "
                               "knowledge of the effective area of the LSST "
                               "mirror.")

        if _MLT_LC_NPZ is None or _MLT_LC_NPZ_NAME != self._mlt_lc_file or _MLT_LC_NPZ.fid is None:
            if not os.path.exists(self._mlt_lc_file):
                catutils_scripts = os.path.join(getPackageDir('sims_catUtils'), 'support_scripts')
                raise RuntimeError("The MLT flaring light curve file:\n"
                                    + "\n%s\n" % self._mlt_lc_file
                                    + "\ndoes not exist."
                                    +"\n\n"
                                    + "Go into %s " % catutils_scripts
                                    + "and run get_mdwarf_flares.sh "
                                    + "to get the data")

            _MLT_LC_NPZ = np.load(self._mlt_lc_file)
            sims_clean_up.targets.append(_MLT_LC_NPZ)
            _MLT_LC_NPZ_NAME = self._mlt_lc_file
            _MLT_LC_TIME_CACHE = {}
            _MLT_LC_FLUX_CACHE = {}

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

            ebv_grid = np.arange(0.0, 7.01, 0.01)
            bb_wavelen = np.arange(200.0, 1500.0, 0.1)
            hc_over_k = 1.4387e7  # nm*K
            temp = 9000.0  # black body temperature in Kelvin
            exp_arg = hc_over_k/(temp*bb_wavelen)
            exp_term = 1.0/(np.exp(exp_arg) - 1.0)
            ln_exp_term = np.log(exp_term)

            # Blackbody f_lambda function;
            # discard normalizing factors; we only care about finding the
            # ratio of fluxes between the case with dust extinction and
            # the case without dust extinction
            log_bb_flambda = -5.0*np.log(bb_wavelen) + ln_exp_term
            bb_flambda = np.exp(log_bb_flambda)
            bb_sed = Sed(wavelen=bb_wavelen, flambda=bb_flambda)

            base_fluxes = self.lsstBandpassDict.fluxListForSed(bb_sed)

            a_x, b_x = bb_sed.setupCCMab()
            self._mlt_dust_lookup = {}
            self._mlt_dust_lookup['ebv'] = ebv_grid
            list_of_bp = self.lsstBandpassDict.keys()
            for bp in list_of_bp:
                self._mlt_dust_lookup[bp] = np.zeros(len(ebv_grid))
            for iebv, ebv_val in enumerate(ebv_grid):
                wv, fl = bb_sed.addCCMDust(a_x, b_x,
                                           ebv=ebv_val,
                                           wavelen=bb_wavelen,
                                           flambda=bb_flambda)

                dusty_bb = Sed(wavelen=wv, flambda=fl)
                dusty_fluxes = self.lsstBandpassDict.fluxListForSed(dusty_bb)
                for ibp, bp in enumerate(list_of_bp):
                    self._mlt_dust_lookup[bp][iebv] = dusty_fluxes[ibp]/base_fluxes[ibp]

        t_init = time.time()-t_before
        t_before = time.time()

        # get the distance to each star in parsecs
        _au_to_parsec = 1.0/206265.0
        dd = _au_to_parsec/parallax

        # get the area of the sphere through which the star's energy
        # is radiating to get to us (in cm^2)
        _cm_per_parsec = 3.08576e18
        sphere_area = 4.0*np.pi*np.power(dd*_cm_per_parsec, 2)

        flux_factor = 1.0/sphere_area

        if isinstance(expmjd, numbers.Number):
            dMags = np.zeros((6, self.num_variable_obj(params)))
        else:
            dMags = np.zeros((6, self.num_variable_obj(params), len(expmjd)))

        mag_name_tuple = ('u', 'g', 'r', 'i', 'z', 'y')
        base_fluxes = {}
        base_mags = {}
        ss = Sed()
        for mag_name in mag_name_tuple:
            if ('lsst_%s' % mag_name in self._actually_calculated_columns or
                'delta_lsst_%s' % mag_name in self._actually_calculated_columns):

                mm = quiescent_mags[mag_name]
                base_mags[mag_name] = mm
                base_fluxes[mag_name] = ss.fluxFromMag(mm)

        lc_name_arr = params['lc'].astype(str)
        lc_names_unique = np.unique(lc_name_arr)

        t_mag_init = time.time()-t_before
        t_before = time.time()

        print('applying MLT to %d -- %d unique' % (len(lc_name_arr), len(lc_names_unique)))
        not_none = 0

        t_use_this = 0.0
        t_flux = 0.0
        t_format_time = 0.0
        t_where = 0.0
        t_load = 0.0

        for lc_name in lc_names_unique:
            t_before = time.time()
            if 'None' in lc_name:
                t_use_this += time.time()-t_before
                continue

            t_before_where = time.time()
            use_this_lc = np.where(np.char.find(lc_name_arr, lc_name)==0)
            t_where += time.time()-t_before_where
            not_none += len(use_this_lc[0])

            lc_name = lc_name.replace('.txt', '')

            # 2017 May 1
            # There isn't supposed to be a 'late_inactive' light curve.
            # Unfortunately, I (Scott Daniel) assigned 'late_inactive'
            # light curves to some of the stars on our database.  Rather
            # than fix the database table (which will take about a week of
            # compute time), I am going to fix the problem here by mapping
            # 'late_inactive' into 'late_active'.
            if 'late' in lc_name:
                lc_name = lc_name.replace('in', '')

            t_before_load = time.time()
            if lc_name in _MLT_LC_TIME_CACHE:
                raw_time_arr = _MLT_LC_TIME_CACHE[lc_name]
            else:
                raw_time_arr = _MLT_LC_NPZ['%s_time' % lc_name]
                _MLT_LC_TIME_CACHE[lc_name] = raw_time_arr

            time_arr = self._survey_start + raw_time_arr
            dt = time_arr.max() - time_arr.min()
            t_load += time.time()-t_before_load

            t_before_time = time.time()
            if isinstance(expmjd, numbers.Number):
                t_interp = (expmjd + params['t0'][use_this_lc]).astype(float)
            else:
                n_obj = len(use_this_lc[0])
                n_time = len(expmjd)
                t_interp = np.ones(shape=(n_obj, n_time))*expmjd
                t_interp += np.array([[tt]*n_time for tt in params['t0'][use_this_lc].astype(float)])

            while t_interp.max() > time_arr.max():
                bad_dexes = np.where(t_interp>time_arr.max())
                t_interp[bad_dexes] -= dt

            t_format_time += time.time()-t_before_time
            t_use_this += time.time()-t_before
            t_before = time.time()

            for i_mag, mag_name in enumerate(mag_name_tuple):
                if ('lsst_%s' % mag_name in self._actually_calculated_columns or
                    'delta_lsst_%s' % mag_name in self._actually_calculated_columns):

                    flux_name = '%s_%s' % (lc_name, mag_name)

                    if flux_name in _MLT_LC_FLUX_CACHE:
                        flux_arr = _MLT_LC_FLUX_CACHE[flux_name]
                    else:
                        flux_arr = _MLT_LC_NPZ[flux_name]
                        _MLT_LC_FLUX_CACHE[flux_name] = flux_arr
                    n_zero = len(np.where(flux_arr==0.0)[0])
                    #print('calculating %s; %d -- %e %e %e --- %d of %d' % (flux_name,len(t_interp),flux_arr.min(),
                    #                                          np.median(flux_arr),flux_arr.max(),n_zero,len(flux_arr)))

                    dflux = np.interp(t_interp, time_arr, flux_arr)

                    if isinstance(expmjd, numbers.Number):
                        dflux *= flux_factor[use_this_lc]
                    else:
                        dflux *= np.array([flux_factor[use_this_lc]]*n_time).transpose()

                    dust_factor = np.interp(ebv[use_this_lc],
                                            self._mlt_dust_lookup['ebv'],
                                            self._mlt_dust_lookup[mag_name])

                    if not isinstance(expmjd, numbers.Number):
                        dust_factor = np.array([dust_factor]*n_time).transpose()

                    dflux *= dust_factor

                    if isinstance(expmjd, numbers.Number):
                        local_base_fluxes = base_fluxes[mag_name][use_this_lc]
                        local_base_mags = base_mags[mag_name][use_this_lc]
                    else:
                        local_base_fluxes = np.array([base_fluxes[mag_name][use_this_lc]]*n_time).transpose()
                        local_base_mags = np.array([base_mags[mag_name][use_this_lc]]*n_time).transpose()

                    dMags[i_mag][use_this_lc] = (ss.magFromFlux(local_base_fluxes + dflux)
                                                 - local_base_mags)

                    #print("local dmags %e %e %e" % (dMags[i_mag][use_this_lc].min(),
                    #                                np.median(dMags[i_mag][use_this_lc]),
                    #                                dMags[i_mag][use_this_lc].max()))

            t_flux += time.time()-t_before

        print('took %.2e\nt_init %.2e\nt_mag_init %.2e\nt_use %.2e\nt_flux %.2e\nformat_time %.2e\nnot_none %d' %
        (time.time()-t_start,t_init,t_mag_init,t_use_this,t_flux,t_format_time,not_none))
        print('t_where %.2e\nt_load %.2e\n' % (t_where, t_load))

        return dMags


class ParametrizedLightCurveMixin(Variability):
    """
    This mixin models variability using parametrized functions fit
    to light curves.

    The parametrized light curves should be stored in an ASCII file
    (or a gzipped ASCII file) whose columns are:

    lc_name -- a string; the original name of the light curve

    n_t_steps -- an int; the number of time steps in the original light curve

    t_span -- a float; t_max - t_min from the original light curve

    n_components -- an int; how many Fourier components were used
    in the parametrization

    chisquared -- this is a series of n_components columns; the nth
    chisquared column is the chisquared of the parametrization after
    n components (i.e. the 5th chisquared value is the chisquared of
    the parametrized light curve with respect to the original light
    curve if you only use the first 5 Fourier components).  This is
    not actually used by this class, but it is expected when parsing
    the parameter file.  It mainly exists if one wishes to perform
    a cut in the parametrization (e.g. only keep as many components
    as are needed to reach some threshold in chisquared/n_t_steps).

    median -- a float; the median flux of the original light cuve

    aa -- a float (see below)
    bb -- a float (see below)
    cc -- a float (see below)
    omega -- a float (see below)
    tau -- a float (see below)

    There will actually be n_components aa, bb, cc, omega, tau
    columns ordered as

    aa_0, bb_0, cc_0, omega_0, tau_0, aa_1, bb_1, cc_1, omega_1, tau_1, ...

    The light curve is parametrized as

    flux = median + \sum_i { aa_i*cos(omega_i*(t-tau_i)) +
                             bb_i*sin(omega_i*(t-tau_i)) +
                             cc_i }
    """

    def load_parametrized_light_curves(self, file_name=None):
        """
        This method will load the parametrized light curve models
        used by the ParametrizedLightCurveMixin and store them in
        a global cache.  It is enough to just run this method from
        any instantiation of ParametrizedLightCurveMixin.

        Parameters
        ----------
        file_name is the absolute path to the file being loaded.
        If None, it will load the default Kepler-based light curve model.
        """
        global _PARAMETRIZED_LC_MODELS
        global _PARAMETRIZED_MODELS_LOADED

        if file_name in _PARAMETRIZED_MODELS_LOADED:
            return

        if len(_PARAMETRIZED_LC_MODELS) == 0:
            sims_clean_up.targets.append(_PARAMETRIZED_LC_MODELS)
            sims_clean_up.targets.append(_PARAMETRIZED_MODELS_LOADED)

        if file_name is None:
            sims_data_dir = getPackageDir('sims_data')
            lc_dir = os.path.join(sims_data_dir, 'catUtilsData')
            file_name = os.path.join(lc_dir, 'kplr_lc_params.txt.gz')

        if file_name.endswith('.gz'):
            open_fn = gzip.open
        else:
            open_fn = open

        if not os.path.exists(file_name):
            if file_name.endswith('kplr_lc_params.txt.gz'):
                download_script = os.path.join(getPackageDir('sims_catUtils'), 'support_scripts',
                                               'get_kepler_light_curves.sh')
                raise RuntimeError('You have not yet downloaded\n%s\n\n' % file_name
                                   + 'Try running the script\n%s' % download_script)
            else:
                raise RuntimeError('The file %s does not exist' % file_name)

        with open_fn(file_name, 'r') as input_file:
            for line in input_file:
                if type(line) == bytes:
                    line = line.decode("utf-8")
                if line[0] == '#':
                    continue
                params = line.strip().split()
                name = params[0]
                tag = int(name.split('_')[0][4:])
                if tag in _PARAMETRIZED_LC_MODELS:
                    # In case multiple sets of models have been loaded that
                    # duplicate identifying integers.
                    raise RuntimeError("You are trying to load light curve with the "
                                       "identifying tag %d.  That has already been " % tag
                                       + "loaded.  I am unsure how to proceed")
                n_c = int(params[3])
                median = float(params[4+n_c])
                local_aa = []
                local_bb = []
                local_cc = []
                local_omega = []
                local_tau = []

                for i_c in range(n_c):
                    base_dex = 5+n_c+i_c*5
                    local_aa.append(float(params[base_dex]))
                    local_bb.append(float(params[base_dex+1]))
                    local_cc.append(float(params[base_dex+2]))
                    local_omega.append(float(params[base_dex+3]))
                    local_tau.append(float(params[base_dex+4]))
                local_aa = np.array(local_aa)
                local_bb = np.array(local_bb)
                local_cc = np.array(local_cc)
                local_omega = np.array(local_omega)
                local_tau = np.array(local_tau)
                _PARAMETRIZED_LC_MODELS[tag] = {}
                _PARAMETRIZED_LC_MODELS[tag]['median'] = median
                _PARAMETRIZED_LC_MODELS[tag]['a'] = local_aa
                _PARAMETRIZED_LC_MODELS[tag]['b'] = local_bb
                _PARAMETRIZED_LC_MODELS[tag]['c'] = local_cc
                _PARAMETRIZED_LC_MODELS[tag]['omega'] = local_omega
                _PARAMETRIZED_LC_MODELS[tag]['tau'] = local_tau

        _PARAMETRIZED_MODELS_LOADED.append(file_name)

    def _calc_dflux(self, lc_id, expmjd):
        """
        Parameters
        ----------
        lc_id is an integer referring to the ID of the light curve in
        the parametrized light curve model (these need to be unique
        across all parametrized light curve catalogs loaded)

        expmjd is either a number or an array referring to the MJD of the
        observations

        Returns
        -------
        baseline_flux is a number indicating the quiescent flux
        of the light curve

        delta_flux is a number or an array of the flux above or below
        the quiescent flux at each of expmjd
        """

        global _PARAMETRIZED_LC_MODELS

        try:
            model = _PARAMETRIZED_LC_MODELS[lc_id]
        except KeyError:
            raise KeyError('A KeyError was raised on the light curve id %d.  ' % lc_id
                           + 'You may not have loaded your parametrized light '
                           + 'curve models, yet.  '
                           + 'See the load_parametrized_light_curves() method in the '
                           + 'ParametrizedLightCurveMixin class')

        tau = model['tau']
        omega = model['omega']
        aa = model['a']
        bb = model['b']
        cc = model['c']

        quiescent_flux = model['median'] + cc.sum()

        omega_t = np.outer(expmjd, omega)
        omega_tau = omega*tau

        # use trig identities to calculate
        # \sum_i a_i*cos(omega_i*(expmjd-tau_i)) + b_i*sin(omega_i*(expmjd-tau_i))
        a_cos_omega_tau = aa*np.cos(omega_tau)
        a_sin_omega_tau = aa*np.sin(omega_tau)
        b_cos_omega_tau = bb*np.cos(omega_tau)
        b_sin_omega_tau = bb*np.sin(omega_tau)

        cos_omega_t = np.cos(omega_t)
        sin_omega_t = np.sin(omega_t)

        delta_flux = np.dot(cos_omega_t, a_cos_omega_tau)
        delta_flux += np.dot(sin_omega_t, a_sin_omega_tau)
        delta_flux += np.dot(sin_omega_t, b_cos_omega_tau)
        delta_flux -= np.dot(cos_omega_t, b_sin_omega_tau)

        if len(delta_flux)==1:
            delta_flux = np.float(delta_flux)
        return quiescent_flux, delta_flux

    @register_method('kplr')  # this 'kplr' tag derives from the fact that default light curves come from Kepler
    def applyParametrizedLightCurve(self, valid_dexes, params, expmjd):

        t_start = time.time()

        if len(params) == 0:
            return np.array([[], [], [], [], [], []])

        n_obj = self.num_variable_obj(params)
        good = np.where(np.logical_not(np.isnan(params['lc'].astype(float))))
        unq_lc_int = np.unique(params['lc'][good])

        # print('applyParamLC %d obj; %d unique' % (n_obj, len(unq_lc_int)))

        if isinstance(expmjd, numbers.Number):
            n_t = 1
            d_mag_out = np.zeros((6, n_obj))
            lc_time = expmjd - params['t0'].astype(float)
        else:
            n_t = len(expmjd)
            d_mag_out = np.zeros((6, n_obj, n_t))
            t0_float = params['t0'].astype(float)
            lc_time = np.zeros(n_t*n_obj)
            i_start = 0
            for i_obj in range(n_obj):
                lc_time[i_start:i_start+n_t] = expmjd - t0_float[i_obj]
                i_start += n_t

        # print('initialized arrays in %e' % (time.time()-t_start))
        # t_assign = 0.0
        # t_flux = 0.0
        # t_use_this = 0.0

        not_none = 0

        lc_int_arr = np.array(params['lc'])
        for lc_int in unq_lc_int:
            if lc_int is None:
                continue
            # t_before = time.time()
            if n_t == 1:
                use_this_lc = np.where(lc_int_arr == lc_int)[0]
                not_none += len(use_this_lc)
            else:
                use_this_lc_unq = np.where(lc_int_arr == lc_int)[0]
                not_none += len(use_this_lc_unq)
                use_this_lc = np.zeros(len(use_this_lc_unq)*n_t, dtype=int)
                template_arange = np.arange(0, n_t, dtype=int)
                for ii, i_lc in enumerate(use_this_lc_unq):
                    use_this_lc[ii*n_t:(ii+1)*n_t] = template_arange + i_lc*n_t

            # t_use_this += time.time()-t_before
            try:
                assert len(use_this_lc) % n_t == 0
            except AssertionError:
                raise RuntimeError("Something went wrong in applyParametrizedLightCurve\n"
                                   "len(use_this_lc) %d ; n_t %d" % (len(use_this_lc), n_t))

            # t_before = time.time()
            q_flux, d_flux = self._calc_dflux(lc_int, lc_time[use_this_lc])
            d_mag = 2.5*np.log10(1.0+d_flux/q_flux)
            # t_flux += time.time()-t_before

            # t_before = time.time()
            if n_t == 1:
                for i_filter in range(6):
                    d_mag_out[i_filter][use_this_lc] = d_mag
            else:
                for i_obj in range(len(use_this_lc)//n_t):
                    i_start = i_obj*n_t
                    obj_dex = use_this_lc_unq[i_obj]
                    for i_filter in range(6):
                        d_mag_out[i_filter][obj_dex] = d_mag[i_start:i_start+n_t]

            # t_assign += time.time()-t_before

        # print('applyParametrized took %.2e\nassignment %.2e\nflux %.2e\nuse %.2e\n' %
        # (time.time()-t_start,t_assign,t_flux,t_use_this))

        print('applying Parametrized LC to %d' % not_none)
        print('took %.2e\n' % (time.time()-t_start))

        return d_mag_out


class ExtraGalacticVariabilityModels(Variability):
    """
    A mixin providing the model for AGN variability.
    """

    @register_method('applyAgn')
    def applyAgn(self, valid_dexes, params, expmjd):

        global _AGN_LC_CACHE

        if len(params) == 0:
            return np.array([[],[],[],[],[],[]])

        if isinstance(expmjd, numbers.Number):
            dMags = np.zeros((6, self.num_variable_obj(params)))
            expmjd_arr = [expmjd]
        else:
            dMags = np.zeros((6, self.num_variable_obj(params), len(expmjd)))
            expmjd_arr = expmjd

        toff_arr = params['t0_mjd'].astype(float)
        seed_arr = params['seed']
        tau_arr = params['agn_tau'].astype(float)
        sfu_arr = params['agn_sfu'].astype(float)
        sfg_arr = params['agn_sfg'].astype(float)
        sfr_arr = params['agn_sfr'].astype(float)
        sfi_arr = params['agn_sfi'].astype(float)
        sfz_arr = params['agn_sfz'].astype(float)
        sfy_arr = params['agn_sfy'].astype(float)

        for i_time, expmjd_val in enumerate(expmjd_arr):
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
                    if _AGN_LC_CACHE[agn_ID]['mjd'] <expmjd_val:
                        resumption = True

                if resumption:
                    rng = copy.deepcopy(_AGN_LC_CACHE[agn_ID]['rng'])
                    start_date = _AGN_LC_CACHE[agn_ID]['mjd']
                    dx_0 = _AGN_LC_CACHE[agn_ID]['dx']
                else:
                    start_date = toff
                    rng = np.random.RandomState(seed)
                    dx_0 = {}
                    for k in sfint:
                        dx_0[k]=0.0

                endepoch = expmjd_val - start_date

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
                    dm_val = (endepoch*(dx1-dx2)+dx2*x1-dx1*x2)/(x1-x2)
                    if isinstance(expmjd, numbers.Number):
                        dMags[ik][ix] = dm_val
                    else:
                        dMags[ik][ix][i_time] = dm_val

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
        dmag = self.applyVariability(varParams)
        if dmag.shape != (6, len(varParams)):
            raise RuntimeError("applyVariability is returning "
                               "an array of shape %s\n" % dmag.shape
                               + "should be (6, %d)" % len(varParams))
        return dmag


class VariabilityStars(_VariabilityPointSources, StellarVariabilityModels,
                       MLTflaringMixin, ParametrizedLightCurveMixin):
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
        dmag = self.applyVariability(varParams)
        if dmag.shape != (6, len(varParams)):
            raise RuntimeError("applyVariability is returning "
                               "an array of shape %s\n" % dmag.shape
                               + "should be (6, %d)" % len(varParams))
        return dmag
