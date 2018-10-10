"""
photUtils -


ljones@astro.washington.edu  (and ajc@astro.washington.edu)

and now (2014 March 28): scott.f.daniel@gmail.com

Collection of utilities to aid usage of Sed and Bandpass with dictionaries.

"""

from builtins import zip
from builtins import range
from builtins import object
import os
import numpy as np
from collections import OrderedDict
from lsst.utils import getPackageDir
from lsst.sims.photUtils import Sed, Bandpass, LSSTdefaults, calcGamma, \
                                calcMagError_m5, calcSNR_m5, PhotometricParameters, magErrorFromSNR, \
                                BandpassDict
from lsst.sims.utils import defaultSpecMap
from lsst.sims.catalogs.decorators import compound
from lsst.sims.photUtils import SedList
from lsst.sims.utils import defaultSpecMap
from lsst.utils import getPackageDir

__all__ = ["PhotometryBase", "PhotometryGalaxies", "PhotometryStars", "PhotometrySSM"]


class PhotometryBase(object):
    """
    This class provides an InstanceCatalog with a member variable photParams,
    which is an instance of the class PhotometricParameters.

    It also provides the method calculateMagnitudeUncertainty, which takes magnitudes,
    a BandpassDict, and an ObservationMetaData as inputs and returns the uncertainties
    on those magnitudes.
    """

    #an object carrying around photometric parameters like readnoise, effective area, plate scale, etc.
    #defaults to LSST values
    photParams = PhotometricParameters()


    def _cacheGamma(self, m5_names, bandpassDict):
        """
        Generate or populate the cache of gamma values used by this InstanceCatalog
        to calculate photometric uncertainties (gamma is defined in equation 5 of
        the LSST overview paper arXiv:0805.2366)

        @param [in] m5_names is a list of the names of keys by which m5 values are
        referred to in the dict self.obs_metadata.m5

        @param [in] bandpassDict is the bandpassDict containing the bandpasses
        corresponding to those m5 values.
        """

        if not hasattr(self, '_gamma_cache'):
            self._gamma_cache = {}

        for mm, bp in zip(m5_names, bandpassDict.values()):
            if mm not in self._gamma_cache and mm in self.obs_metadata.m5:
                self._gamma_cache[mm] = calcGamma(bp, self.obs_metadata.m5[mm], photParams=self.photParams)


    def _magnitudeUncertaintyGetter(self, column_name_list, m5_name_list, bandpassDict_name):
        """
        Generic getter for magnitude uncertainty columns.

        Columns must be named 'sigma_xx' where 'xx' is the column name of
        the associated magnitude

        @param [in] column_name_list is the list of magnitude column names
        associated with the uncertainties calculated by this getter
        (the 'xx' in the 'sigma_xx' above)

        @param [in] m5_name_list are the keys to the self.obs_metadata.m5 dict
        corresponding to the bandpasses in column_names.  For example: in
        the case of galaxies, the magnitude columns

        column_names = ['uBulge', 'gBulge', 'rBulge', 'iBulge', 'zBulge', 'yBulge']

        may correspond to m5 values keyed to

        m5_name_list = ['u', 'g', 'r', 'i', 'z', 'y']

        @param [in] bandpassDict_name is a string indicating the name of
        the InstanceCatalog member variable containing the BandpassDict
        to be used when calculating these magnitude uncertainties.
        The BandpassDict itself will be accessed using
        getattr(self, bandpassDict_name)

        @param [out] returns a 2-D numpy array in which the first index
        is the uncertainty column and the second index is the catalog row
        (i.e. output suitable for a getter in an InstanceCatalog)
        """

        # make sure that the magnitudes associated with any requested
        # uncertainties actually are calculated
        num_elements = None
        mag_dict = {}
        for name in column_name_list:
            if 'sigma_%s' % name in self._actually_calculated_columns:
                mag_dict[name] = self.column_by_name(name)
                if num_elements is None:
                    num_elements = len(mag_dict[name])

        # These lines must come after the preceding lines;
        # the bandpassDict will not be loaded until the magnitude
        # getters are called
        bandpassDict = getattr(self, bandpassDict_name)
        self._cacheGamma(m5_name_list, bandpassDict)


        output = []

        for name, m5_name, bp in zip(column_name_list, m5_name_list, bandpassDict.values()):
            if 'sigma_%s' % name not in self._actually_calculated_columns:
                output.append(np.ones(num_elements)*np.NaN)
            else:
                try:
                    m5 = self.obs_metadata.m5[m5_name]
                    gamma = self._gamma_cache[m5_name]

                    sigma_list, gamma = calcMagError_m5(mag_dict[name], bp, m5, self.photParams, gamma=gamma)

                    output.append(sigma_list)

                except KeyError as kk:
                    msg = 'You got the KeyError: %s' % kk.args[0]
                    raise KeyError('%s \n' % msg \
                                   + 'Is it possible your ObservationMetaData does not have the proper\n'
                                   'm5 values defined?')


        return np.array(output)


    @compound('sigma_lsst_u','sigma_lsst_g','sigma_lsst_r','sigma_lsst_i',
              'sigma_lsst_z','sigma_lsst_y')
    def get_lsst_photometric_uncertainties(self):
        """
        Getter for photometric uncertainties associated with lsst bandpasses
        """

        return self._magnitudeUncertaintyGetter(['lsst_u', 'lsst_g', 'lsst_r',
                                                 'lsst_i', 'lsst_z', 'lsst_y'],
                                                ['u', 'g', 'r', 'i', 'z', 'y'],
                                                 'lsstBandpassDict')


    def calculateVisibility(self, magFilter, sigma=0.1, randomSeed=None, pre_generate_randoms=False):
        """
        Determine (probabilistically) whether a source was detected or not.

        The 'completeness' of source detection at any magnitude is calculated by
          completeness = (1 + e^^(magFilter-m5)/sigma)^(-1)
        For each source with a magnitude magFilter, if a uniform random number [0-1)
          is less than or equal to 'completeness', then it is counted as "detected".
        See equation 24, figure 8 and table 5 from SDSS completeness analysis in
        http://iopscience.iop.org/0004-637X/794/2/120/pdf/apj_794_2_120.pdf
        "THE SLOAN DIGITAL SKY SURVEY COADD: 275 deg2 OF DEEP SLOAN DIGITAL SKY SURVEY IMAGING ON STRIPE 82"

        @ param [in] magFilter is the magnitude of the object in the observed filter.

        @ param [in] sigma is the FWHM of the distribution (default = 0.1)

        @ param [in] randomSeed is an option to set a random seed (default None)
        @ param [in] pre_generate_randoms is an option (default False) to pre-generate a series of 12,000,000 random numbers
           for use throughout the visibility calculation [the random numbers used are randoms[objId]].

        @ param [out] visibility (None/1).
        """
        if len(magFilter) == 0:
            return np.array([])
        # Calculate the completeness at the magnitude of each object.
        completeness = 1.0 / (1 + np.exp((magFilter - self.obs_metadata.m5[self.obs_metadata.bandpass])/sigma))
        # Seed numpy if desired and not previously done.
        if (randomSeed is not None) and (not hasattr(self, 'ssm_random_seeded')):
            np.random.seed(randomSeed)
            self.ssm_random_seeded = True
        # Pre-generate random numbers, if desired and not previously done.
        if pre_generate_randoms and not hasattr(self, 'ssm_randoms'):
            self.ssm_randoms = np.random.rand(14000000)
        # Calculate probability values to compare to completeness.
        if hasattr(self, 'ssm_randoms'):
            # Grab the random numbers from self.randoms.
            probability = self.ssm_randoms[self.column_by_name('objId')]
        else:
            probability = np.random.random_sample(len(magFilter))
        # Compare the random number to the completeness.
        visibility = np.where(probability <= completeness, 1, None)
        return visibility

    def _variabilityGetter(self, columnNames):
        """
        Find columns named 'delta_*' and return them to be added
        to '*' magnitude columns (i.e. look for delta_lsst_u so that
        it can be added to lsst_u)

        Parameters
        ----------
        columnNames is a list of the quiescent columns (lsst_u in the
        example above) whose deltas we are looking for

        Returns
        -------
        A numpy array in which each row is a delta magnitude and each
        column is an astrophysical object/database row
        """

        num_obj = len(self.column_by_name(self.db_obj.idColKey))
        delta = []

        # figure out which of these columns we are actually calculating
        indices = [ii for ii, name in enumerate(columnNames)
                   if name in self._actually_calculated_columns]

        for ix, columnName in enumerate(columnNames):
            if indices is None or ix in indices:
                delta_name = 'delta_' + columnName
                if delta_name in self._all_available_columns:
                    delta.append(self.column_by_name(delta_name))
                else:
                    delta.append(np.zeros(num_obj))
            else:
                delta.append(np.zeros(num_obj))

        return np.array(delta)

class PhotometryGalaxies(PhotometryBase):
    """
    This mixin provides the code necessary for calculating the component magnitudes associated with
    galaxies.  It assumes that we want LSST filters.
    """

    def _hasCosmoDistMod(self):
        """
        Determine whether or not this InstanceCatalog has a column
        specifically devoted to the cosmological distance modulus.
        """
        if 'cosmologicalDistanceModulus' in self._all_available_columns:
            return True
        return False


    def _loadBulgeSedList(self, wavelen_match):
        """
        Load a SedList of galaxy bulge Seds.
        The list will be stored in the variable self._bulgeSedList.

        @param [in] wavelen_match is the wavelength grid (in nm)
        on which the Seds are to be sampled.
        """

        sedNameList = self.column_by_name('sedFilenameBulge')
        magNormList = self.column_by_name('magNormBulge')
        redshiftList = self.column_by_name('redshift')
        internalAvList = self.column_by_name('internalAvBulge')
        cosmologicalDimming = not self._hasCosmoDistMod()

        if len(sedNameList)==0:
            return np.ones((0))

        if not hasattr(self, '_bulgeSedList'):
            self._bulgeSedList = SedList(sedNameList, magNormList,
                                         internalAvList=internalAvList,
                                         redshiftList=redshiftList,
                                         cosmologicalDimming=cosmologicalDimming,
                                         wavelenMatch=wavelen_match,
                                         fileDir=getPackageDir('sims_sed_library'),
                                         specMap=defaultSpecMap)
        else:
            self._bulgeSedList.flush()
            self._bulgeSedList.loadSedsFromList(sedNameList, magNormList,
                                               internalAvList=internalAvList,
                                               redshiftList=redshiftList)


    def _loadDiskSedList(self, wavelen_match):
        """
        Load a SedList of galaxy disk Seds.
        The list will be stored in the variable self._bulgeSedList.

        @param [in] wavelen_match is the wavelength grid (in nm)
        on which the Seds are to be sampled.
        """

        sedNameList = self.column_by_name('sedFilenameDisk')
        magNormList = self.column_by_name('magNormDisk')
        redshiftList = self.column_by_name('redshift')
        internalAvList = self.column_by_name('internalAvDisk')
        cosmologicalDimming = not self._hasCosmoDistMod()

        if len(sedNameList)==0:
            return np.ones((0))

        if not hasattr(self, '_diskSedList'):
            self._diskSedList = SedList(sedNameList, magNormList,
                                        internalAvList=internalAvList,
                                        redshiftList=redshiftList,
                                        cosmologicalDimming=cosmologicalDimming,
                                        wavelenMatch=wavelen_match,
                                        fileDir=getPackageDir('sims_sed_library'),
                                        specMap=defaultSpecMap)
        else:
            self._diskSedList.flush()
            self._diskSedList.loadSedsFromList(sedNameList, magNormList,
                                               internalAvList=internalAvList,
                                               redshiftList=redshiftList)


    def _loadAgnSedList(self, wavelen_match):
        """
        Load a SedList of galaxy AGN Seds.
        The list will be stored in the variable self._bulgeSedList.

        @param [in] wavelen_match is the wavelength grid (in nm)
        on which the Seds are to be sampled.
        """

        sedNameList = self.column_by_name('sedFilenameAgn')
        magNormList = self.column_by_name('magNormAgn')
        redshiftList = self.column_by_name('redshift')
        cosmologicalDimming = not self._hasCosmoDistMod()

        if len(sedNameList)==0:
            return np.ones((0))

        if not hasattr(self, '_agnSedList'):
            self._agnSedList = SedList(sedNameList, magNormList,
                                       redshiftList=redshiftList,
                                       cosmologicalDimming=cosmologicalDimming,
                                       wavelenMatch=wavelen_match,
                                       fileDir=getPackageDir('sims_sed_library'),
                                       specMap=defaultSpecMap)
        else:
            self._agnSedList.flush()
            self._agnSedList.loadSedsFromList(sedNameList, magNormList,
                                               redshiftList=redshiftList)


    def sum_magnitudes(self, disk = None, bulge = None, agn = None):
        """
        Sum the component magnitudes of a galaxy and return the answer

        @param [in] disk is the disk magnitude must be a numpy array or a float

        @param [in] bulge is the bulge magnitude must be a numpy array or a float

        @param [in] agn is the agn magnitude must be a numpy array or a float

        @param [out] outMag is the total magnitude of the galaxy
        """
        with np.errstate(divide='ignore', invalid='ignore'):
            baselineType = type(None)
            if not isinstance(disk, type(None)):
                baselineType = type(disk)
                if baselineType == np.ndarray:
                    elements=len(disk)

            if not isinstance(bulge, type(None)):
                if baselineType == type(None):
                    baselineType = type(bulge)
                    if baselineType == np.ndarray:
                        elements = len(bulge)
                elif not isinstance(bulge, baselineType):
                    raise RuntimeError("All non-None arguments of sum_magnitudes need to be " +
                                       "of the same type (float or numpy array)")

            elif not isinstance(agn, type(None)):
                if baseLineType == type(None):
                    baselineType = type(agn)
                    if baselineType == np.ndarray:
                        elements = len(agn)
                elif not isinstance(agn, baselineType):
                    raise RuntimeError("All non-None arguments of sum_magnitudes need to be " +
                                       "of the same type (float or numpy array)")

            if baselineType is not float and \
               baselineType is not np.ndarray and \
               baselineType is not np.float and \
               baselineType is not np.float64:

                raise RuntimeError("Arguments of sum_magnitudes need to be " +
                                   "either floats or numpy arrays; you appear to have passed %s " % baselineType)

            mm_0 = 22.
            tol = 1.0e-30

            if baselineType == np.ndarray:
                nn = np.zeros(elements)
            else:
                nn = 0.0

            if disk is not None:
                nn += np.where(np.isnan(disk), 0.0, np.power(10, -0.4*(disk - mm_0)))

            if bulge is not None:
                nn += np.where(np.isnan(bulge), 0.0, np.power(10, -0.4*(bulge - mm_0)))

            if agn is not None:
                nn += np.where(np.isnan(agn), 0.0, np.power(10, -0.4*(agn - mm_0)))

            if baselineType == np.ndarray:
                # according to this link
                # http://stackoverflow.com/questions/25087769/runtimewarning-divide-by-zero-error-how-to-avoid-python-numpy
                # we will still get a divide by zero error from log10, but np.where will be
                # circumventing the offending value, so it is probably okay
                return np.where(nn>tol, -2.5*np.log10(nn) + mm_0, np.NaN)
            else:
                if nn>tol:
                    return -2.5*np.log10(nn) + mm_0
                else:
                    return np.NaN


    def _quiescentMagnitudeGetter(self, componentName, bandpassDict, columnNameList):
        """
        A generic getter for quiescent magnitudes of galaxy components.

        @param [in] componentName is either 'bulge', 'disk', or 'agn'

        @param [in] bandpassDict is a BandpassDict of the bandpasses
        in which to calculate the magnitudes

        @param [in] columnNameList is a list of the columns corresponding to
        these magnitudes (for purposes of applying variability).

        @param [out] magnitudes is a 2-D numpy array of magnitudes in which
        rows correspond to bandpasses and columns correspond to astronomical
        objects.
        """

        # figure out which of these columns we are actually calculating
        indices = [ii for ii, name in enumerate(columnNameList)
                   if name in self._actually_calculated_columns]

        if len(indices) == len(columnNameList):
            indices = None

        if componentName == 'bulge':
            self._loadBulgeSedList(bandpassDict.wavelenMatch)
            if not hasattr(self, '_bulgeSedList'):
                sedList = None
            else:
                sedList = self._bulgeSedList
        elif componentName == 'disk':
            self._loadDiskSedList(bandpassDict.wavelenMatch)
            if not hasattr(self, '_diskSedList'):
                sedList = None
            else:
                sedList = self._diskSedList
        elif componentName == 'agn':
            self._loadAgnSedList(bandpassDict.wavelenMatch)
            if not hasattr(self, '_agnSedList'):
                sedList = None
            else:
                sedList = self._agnSedList
        else:
            raise RuntimeError('_quiescentMagnitudeGetter does not understand component %s ' \
                               % componentName)

        if sedList is None:
            magnitudes = np.ones((len(columnNameList), 0))
        else:
            magnitudes = bandpassDict.magListForSedList(sedList, indices=indices).transpose()

        if self._hasCosmoDistMod():
            cosmoDistMod = self.column_by_name('cosmologicalDistanceModulus')
            if len(cosmoDistMod)>0:
               for ix in range(magnitudes.shape[0]):
                   magnitudes[ix] += cosmoDistMod

        return magnitudes


    @compound('sigma_uBulge', 'sigma_gBulge', 'sigma_rBulge',
              'sigma_iBulge', 'sigma_zBulge', 'sigma_yBulge')
    def get_photometric_uncertainties_bulge(self):
        """
        Getter for photometric uncertainties associated with galaxy bulges
        """

        return self._magnitudeUncertaintyGetter(['uBulge', 'gBulge', 'rBulge',
                                                'iBulge', 'zBulge', 'yBulge'],
                                                ['u', 'g', 'r', 'i', 'z', 'y'],
                                                'lsstBandpassDict')


    @compound('sigma_uDisk', 'sigma_gDisk', 'sigma_rDisk',
              'sigma_iDisk', 'sigma_zDisk', 'sigma_yDisk')
    def get_photometric_uncertainties_disk(self):
        """
        Getter for photometeric uncertainties associated with galaxy disks
        """

        return self._magnitudeUncertaintyGetter(['uDisk', 'gDisk', 'rDisk',
                                                'iDisk', 'zDisk', 'yDisk'],
                                                ['u', 'g', 'r', 'i', 'z', 'y'],
                                                'lsstBandpassDict')


    @compound('sigma_uAgn', 'sigma_gAgn', 'sigma_rAgn',
              'sigma_iAgn', 'sigma_zAgn', 'sigma_yAgn')
    def get_photometric_uncertainties_agn(self):
        """
        Getter for photometric uncertainties associated with Agn
        """

        return self._magnitudeUncertaintyGetter(['uAgn', 'gAgn', 'rAgn',
                                                'iAgn', 'zAgn', 'yAgn'],
                                                ['u', 'g', 'r', 'i', 'z', 'y'],
                                                'lsstBandpassDict')


    @compound('uBulge', 'gBulge', 'rBulge', 'iBulge', 'zBulge', 'yBulge')
    def get_lsst_bulge_mags(self):
        """
        Getter for bulge magnitudes in LSST bandpasses
        """

        # load a BandpassDict of LSST bandpasses, if not done already
        if not hasattr(self, 'lsstBandpassDict'):
            self.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()

        # actually calculate the magnitudes
        mag = self._quiescentMagnitudeGetter('bulge', self.lsstBandpassDict,
                                             self.get_lsst_bulge_mags._colnames)

        mag += self._variabilityGetter(self.get_lsst_bulge_mags._colnames)
        return mag


    @compound('uDisk', 'gDisk', 'rDisk', 'iDisk', 'zDisk', 'yDisk')
    def get_lsst_disk_mags(self):
        """
        Getter for galaxy disk magnitudes in the LSST bandpasses
        """

        # load a BandpassDict of LSST bandpasses, if not done already
        if not hasattr(self, 'lsstBandpassDict'):
            self.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()

        # actually calculate the magnitudes
        mag = self._quiescentMagnitudeGetter('disk', self.lsstBandpassDict,
                                             self.get_lsst_disk_mags._colnames)

        mag += self._variabilityGetter(self.get_lsst_disk_mags._colnames)
        return mag


    @compound('uAgn', 'gAgn', 'rAgn', 'iAgn', 'zAgn', 'yAgn')
    def get_lsst_agn_mags(self):
        """
        Getter for AGN magnitudes in the LSST bandpasses
        """

        # load a BandpassDict of LSST bandpasses, if not done already
        if not hasattr(self, 'lsstBandpassDict'):
            self.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()

        # actually calculate the magnitudes
        mag = self._quiescentMagnitudeGetter('agn', self.lsstBandpassDict,
                                             self.get_lsst_agn_mags._colnames)

        mag += self._variabilityGetter(self.get_lsst_agn_mags._colnames)
        return mag

    @compound('lsst_u', 'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y')
    def get_lsst_total_mags(self):
        """
        Getter for total galaxy magnitudes in the LSST bandpasses
        """

        idList = self.column_by_name('uniqueId')
        numObj = len(idList)
        output = []

        # Loop over the columns calculated by this getter.  For each
        # column, calculate the bluge, disk, and agn magnitude in the
        # corresponding bandpass, then sum them using the
        # sum_magnitudes method.
        for columnName in self.get_lsst_total_mags._colnames:
            if columnName not in self._actually_calculated_columns:
                sub_list = [np.NaN]*numObj
            else:
                bandpass = columnName[-1]
                bulge = self.column_by_name('%sBulge' % bandpass)
                disk = self.column_by_name('%sDisk' % bandpass)
                agn = self.column_by_name('%sAgn' % bandpass)
                sub_list = self.sum_magnitudes(bulge=bulge, disk=disk, agn=agn)

            output.append(sub_list)
        return np.array(output)




class PhotometryStars(PhotometryBase):
    """
    This mixin provides the infrastructure for doing photometry on stars

    It assumes that we want LSST filters.
    """

    def _loadSedList(self, wavelen_match):
        """
        Method to load the member variable self._sedList, which is a SedList.
        If self._sedList does not already exist, this method sets it up.
        If it does already exist, this method flushes its contents and loads a new
        chunk of Seds.
        """

        sedNameList = self.column_by_name('sedFilename')
        magNormList = self.column_by_name('magNorm')
        galacticAvList = self.column_by_name('galacticAv')

        if len(sedNameList)==0:
            return np.ones((0))

        if not hasattr(self, '_sedList'):
            self._sedList = SedList(sedNameList, magNormList,
                                    galacticAvList=galacticAvList,
                                    wavelenMatch=wavelen_match,
                                    fileDir=getPackageDir('sims_sed_library'),
                                    specMap=defaultSpecMap)
        else:
            self._sedList.flush()
            self._sedList.loadSedsFromList(sedNameList, magNormList,
                                          galacticAvList=galacticAvList)


    def _quiescentMagnitudeGetter(self, bandpassDict, columnNameList):
        """
        This method gets the magnitudes for an InstanceCatalog, returning them
        in a 2-D numpy array in which rows correspond to bandpasses and columns
        correspond to astronomical objects.

        @param [in] bandpassDict is a BandpassDict containing the bandpasses
        whose magnitudes are to be calculated

        @param [in] columnNameList is a list of the names of the magnitude columns
        being calculated

        @param [out] magnitudes is a 2-D numpy array of magnitudes in which
        rows correspond to bandpasses in bandpassDict and columns correspond
        to astronomical objects.
        """

        # figure out which of these columns we are actually calculating
        indices = [ii for ii, name in enumerate(columnNameList)
                   if name in self._actually_calculated_columns]

        if len(indices) == len(columnNameList):
            indices = None

        self._loadSedList(bandpassDict.wavelenMatch)

        if not hasattr(self, '_sedList'):
            magnitudes = np.ones((len(columnNameList),0))
        else:
            magnitudes = bandpassDict.magListForSedList(self._sedList, indices=indices).transpose()

        return magnitudes

    @compound('quiescent_lsst_u', 'quiescent_lsst_g', 'quiescent_lsst_r',
              'quiescent_lsst_i', 'quiescent_lsst_z', 'quiescent_lsst_y')
    def get_quiescent_lsst_magnitudes(self):

        if not hasattr(self, 'lsstBandpassDict'):
            self.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()

        return self._quiescentMagnitudeGetter(self.lsstBandpassDict,
                                              self.get_quiescent_lsst_magnitudes._colnames)

    @compound('lsst_u','lsst_g','lsst_r','lsst_i','lsst_z','lsst_y')
    def get_lsst_magnitudes(self):
        """
        getter for LSST stellar magnitudes
        """

        magnitudes = np.array([self.column_by_name('quiescent_lsst_u'),
                                  self.column_by_name('quiescent_lsst_g'),
                                  self.column_by_name('quiescent_lsst_r'),
                                  self.column_by_name('quiescent_lsst_i'),
                                  self.column_by_name('quiescent_lsst_z'),
                                  self.column_by_name('quiescent_lsst_y')])

        delta = self._variabilityGetter(self.get_lsst_magnitudes._colnames)
        magnitudes += delta

        return magnitudes

class PhotometrySSM(PhotometryBase):
    """
    A mixin to calculate photometry for solar system objects.
    """
    # Because solar system objects will not have dust extinctions, we should be able to read in every
    # SED exactly once, calculate the colors and magnitudes, and then get actual magnitudes by adding
    # an offset based on magNorm.

    def _quiescentMagnitudeGetter(self, bandpassDict, columnNameList, bandpassTag='lsst'):
        """
        Method that actually does the work calculating magnitudes for solar system objects.

        Because solar system objects have no dust extinction, this method works by loading
        each unique Sed once, normalizing it, calculating its magnitudes in the desired
        bandpasses, and then storing the normalizing magnitudes and the bandpass magnitudes
        in a dict.  Magnitudes for subsequent objects with identical Seds will be calculated
        by adding an offset to the magnitudes.  The offset is determined by comparing normalizing
        magnitues.

        @param [in] bandpassDict is an instantiation of BandpassDict representing the bandpasses
        to be integrated over

        @param [in] columnNameList is a list of the names of the columns being calculated
        by this getter

        @param [in] bandpassTag (optional) is a string indicating the name of the bandpass system
        (i.e. 'lsst', 'sdss', etc.).  This is in case the user wants to calculate the magnitudes
        in multiple systems simultaneously.  In that case, the dict will store magnitudes for each
        Sed in each magnitude system separately.

        @param [out] a numpy array of magnitudes corresponding to bandpassDict.
        """

        # figure out which of these columns we are actually calculating
        indices = [ii for ii, name in enumerate(columnNameList)
                   if name in self._actually_calculated_columns]

        if len(indices) == len(columnNameList):
            indices = None

        if not hasattr(self, '_ssmMagDict'):
            self._ssmMagDict = {}
            self._ssmMagNormDict = {}
            self._file_dir = getPackageDir('sims_sed_library')
            self._spec_map = defaultSpecMap
            self._normalizing_bandpass = Bandpass()
            self._normalizing_bandpass.imsimBandpass()

        sedNameList = self.column_by_name('sedFilename')
        magNormList = self.column_by_name('magNorm')

        if len(sedNameList)==0:
            # need to return something when InstanceCatalog goes through
            # it's "dry run" to determine what columns are required from
            # the database
            return np.zeros((len(bandpassDict.keys()),0))

        magListOut = []

        for sedName, magNorm in zip(sedNameList, magNormList):
            magTag = bandpassTag+'_'+sedName
            if sedName not in self._ssmMagNormDict or magTag not in self._ssmMagDict:
                dummySed = Sed()
                dummySed.readSED_flambda(os.path.join(self._file_dir, self._spec_map[sedName]))
                fnorm = dummySed.calcFluxNorm(magNorm, self._normalizing_bandpass)
                dummySed.multiplyFluxNorm(fnorm)
                magList = bandpassDict.magListForSed(dummySed, indices=indices)
                self._ssmMagDict[magTag] = magList
                self._ssmMagNormDict[sedName] = magNorm
            else:
                dmag = magNorm - self._ssmMagNormDict[sedName]
                magList = self._ssmMagDict[magTag] + dmag
            magListOut.append(magList)

        return np.array(magListOut).transpose()


    @compound('lsst_u','lsst_g','lsst_r','lsst_i','lsst_z','lsst_y')
    def get_lsst_magnitudes(self):
        """
        getter for LSST magnitudes of solar system objects
        """
        if not hasattr(self, 'lsstBandpassDict'):
            self.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()

        return self._quiescentMagnitudeGetter(self.lsstBandpassDict, self.get_lsst_magnitudes._colnames)


    def get_magFilter(self):
        """
        Generate the magnitude in the filter of the observation.
        """
        magFilter = 'lsst_' + self.obs_metadata.bandpass
        return self.column_by_name(magFilter)

    def get_magSNR(self):
        """
        Calculate the SNR for the observation, given m5 from obs_metadata and the trailing losses.
        """
        magFilter = self.column_by_name('magFilter')
        bandpass = self.lsstBandpassDict[self.obs_metadata.bandpass]
        # Get m5 for the visit
        m5 = self.obs_metadata.m5[self.obs_metadata.bandpass]
        # Adjust the magnitude of the source for the trailing losses.
        dmagSNR = self.column_by_name('dmagTrailing')
        magObj = magFilter - dmagSNR
        if len(magObj) == 0:
            snr = []
        else:
            snr, gamma = calcSNR_m5(magObj, bandpass, m5, self.photParams)
        return snr

    def get_visibility(self):
        """
        Generate a None/1 flag indicating whether the object was detected or not.

        Sets the random seed for 'calculateVisibility' using the obs_metadata.obsHistId
        """
        magFilter = self.column_by_name('magFilter')
        dmagDetect = self.column_by_name('dmagDetection')
        magObj = magFilter - dmagDetect
        # Adjusted m5 value, accounting for the fact these are moving objects.
        mjdSeed = np.int(self.obs_metadata.mjd.TAI * 1000000) % 4294967295
        visibility = self.calculateVisibility(magObj, randomSeed=mjdSeed, pre_generate_randoms=True)
        return visibility


    @compound('dmagTrailing', 'dmagDetection')
    def get_ssm_dmag(self):
        """
        This getter will calculate:

        dmagTrailing: the offset in m5 used to represent the loss in signal to noise
        resulting from the fact that the object's motion smears out its PSF

        dmagDetection: the offset in m5 used to represent the shift in detection
        threshold resulting from the fact that the object's motion smears out
        its PSF
        """

        if self.obs_metadata.seeing is None:
            raise RuntimeError("Cannot calculate dmagTraling/dmagDetection. "
                               "Your catalog's ObservationMetaData does not "
                               "specify seeing.")

        if len(self.obs_metadata.seeing)>1:
            valueList = list(self.obs_metadata.seeing.values())
            for ix in range(1, len(valueList)):
                if np.abs(valueList[ix]-valueList[0])>0.0001:

                    raise RuntimeError("dmagTrailing/dmagDetection calculation is confused. "
                                       "Your catalog's ObservationMetaData contains multiple "
                                       "seeing values.  Re-create your catalog with only one seeing value.")

        if not hasattr(self, 'photParams') or self.photParams is None:
            raise RuntimeError("You cannot calculate dmagTrailing/dmagDetection. "
                               "Your catalog does not have an associated PhotometricParameters "
                               "member variable.  It is impossible to know what the exposure time is.")

        dradt = self.column_by_name('velRa') # in radians per day (actual sky velocity;
                                             # i.e., no need to divide by cos(dec))

        ddecdt = self.column_by_name('velDec') # in radians per day

        if len(dradt)==0:
            return np.zeros((2,0))

        a_trail = 0.76
        b_trail = 1.16
        a_det = 0.42
        b_det = 0.00
        seeing = self.obs_metadata.seeing[self.obs_metadata.bandpass] # this will be in arcsec
        texp = self.photParams.nexp*self.photParams.exptime  # in seconds
        velocity = np.sqrt(np.power(np.degrees(dradt),2) + np.power(np.degrees(ddecdt),2)) # in degrees/day
        x = velocity*texp/(24.0*seeing)
        xsq = np.power(x,2)
        dmagTrail = 1.25*np.log10(1.0 + a_trail * xsq/(1.0+b_trail*x))
        dmagDetect = 1.25*np.log10(1.0 + a_det * xsq/(1.0 + b_det*x))

        return np.array([dmagTrail, dmagDetect])
