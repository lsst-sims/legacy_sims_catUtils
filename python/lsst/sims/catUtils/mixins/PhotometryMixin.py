"""
photUtils -


ljones@astro.washington.edu  (and ajc@astro.washington.edu)

and now (2014 March 28): scott.f.daniel@gmail.com

Collection of utilities to aid usage of Sed and Bandpass with dictionaries.

"""

import os
import numpy
from collections import OrderedDict
from lsst.utils import getPackageDir
from lsst.sims.photUtils import Sed, Bandpass, LSSTdefaults, calcGamma, \
                                calcMagError_m5, PhotometricParameters, magErrorFromSNR, \
                                BandpassDict
from lsst.sims.utils import defaultSpecMap
from lsst.sims.catalogs.measures.instance import compound
from lsst.sims.photUtils import SedList

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

    def calculateMagnitudeUncertainty(self, magnitudes, bandpassDict, obs_metadata=None):
        """
        Calculate the uncertainty in magnitudes using the model from equation (5) of arXiv:0805.2366

        @param [in] magnitudes is a numpy array containing the object magnitudes.  Every row corresponds to
        a bandpass, which is to say that magnitudes[i][j] is the magnitude of the jth object in the
        bandpass characterized by self.bandpassDict.values()[i]

        @param [in] bandpassDict is a BandpassDict characterizing the bandpasses being used

        @param [in] obs_metadata is the metadata of this observation (mostly desired because
        it will contain information about m5, the magnitude at which objects are detected
        at the 5-sigma limit in each bandpass)

        @param [out] a 1-dimensional numpy array containing the magntiude uncertainty
        corresponding to the magnitudes passed into this method.
        """

        if obs_metadata is None:
            raise RuntimeError("Need to pass an ObservationMetaData into calculatePhotometricUncertainty")

        if magnitudes.shape[0] != len(bandpassDict):
            raise RuntimeError("Passed %d magnitudes to " % magnitudes.shape[0] + \
                                " PhotometryBase.calculatePhotometricUncertainty; " + \
                                "needed %d " % len(bandpassDict))

        #if we have not run this method before, calculate and cache the m5 and gamma parameter
        #values (see eqn 5 of arXiv:0805.2366) for future use
        m5Defaults = None

        if not hasattr(self, '_gammaList') or \
        len(self._gammaList) != len(bandpassDict):

            mm = []
            gg = []
            for b in bandpassDict.keys():
                if b in obs_metadata.m5:
                    mm.append(obs_metadata.m5[b])
                    gg.append(calcGamma(bandpassDict[b], obs_metadata.m5[b],
                              photParams=self.photParams))
                else:
                    if m5Defaults is None:
                        m5Defaults = LSSTdefaults()

                    try:
                        mm.append(m5Defaults.m5(b))
                        gg.append(m5Defaults.gamma(b))
                    except:
                        raise RuntimeError("No way to calculate gamma or m5 for filter %s " % b)

            self._m5List = numpy.array(mm)
            self._gammaList = numpy.array(gg)

        error = calcMagError_m5(magnitudes, bandpassDict.values(), self._m5List,
                                gamma=self._gammaList, photParams=self.photParams)

        return error


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
            return numpy.ones((0))

        if not hasattr(self, '_bulgeSedList'):
            self._bulgeSedList = SedList(sedNameList, magNormList,
                                               internalAvList=internalAvList,
                                               redshiftList=redshiftList,
                                               cosmologicalDimming=cosmologicalDimming,
                                               wavelenMatch=wavelen_match)
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
            return numpy.ones((0))

        if not hasattr(self, '_diskSedList'):
            self._diskSedList = SedList(sedNameList, magNormList,
                                               internalAvList=internalAvList,
                                               redshiftList=redshiftList,
                                               cosmologicalDimming=cosmologicalDimming,
                                               wavelenMatch=wavelen_match)
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
            return numpy.ones((0))

        if not hasattr(self, '_agnSedList'):
            self._agnSedList = SedList(sedNameList, magNormList,
                                               redshiftList=redshiftList,
                                               cosmologicalDimming=cosmologicalDimming,
                                               wavelenMatch=wavelen_match)
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

        baselineType = type(None)
        if not isinstance(disk, type(None)):
            baselineType = type(disk)
            if baselineType == numpy.ndarray:
                elements=len(disk)

        if not isinstance(bulge, type(None)):
            if baselineType == type(None):
                baselineType = type(bulge)
                if baselineType == numpy.ndarray:
                    elements = len(bulge)
            elif not isinstance(bulge, baselineType):
                raise RuntimeError("All non-None arguments of sum_magnitudes need to be " +
                                   "of the same type (float or numpy array)")

        elif not isinstance(agn, type(None)):
            if baseLineType == type(None):
                baselineType = type(agn)
                if baselineType == numpy.ndarray:
                    elements = len(agn)
            elif not isinstance(agn, baselineType):
                raise RuntimeError("All non-None arguments of sum_magnitudes need to be " +
                                   "of the same type (float or numpy array)")

        if baselineType is not float and \
           baselineType is not numpy.ndarray and \
           baselineType is not numpy.float and \
           baselineType is not numpy.float64:

            raise RuntimeError("Arguments of sum_magnitudes need to be " +
                               "either floats or numpy arrays; you appear to have passed %s " % baselineType)

        mm_0 = 22.
        tol = 1.0e-30

        if baselineType == numpy.ndarray:
            nn = numpy.zeros(elements)
        else:
            nn = 0.0

        if disk is not None:
            nn += numpy.where(numpy.isnan(disk), 0.0, numpy.power(10, -0.4*(disk - mm_0)))

        if bulge is not None:
            nn += numpy.where(numpy.isnan(bulge), 0.0, numpy.power(10, -0.4*(bulge - mm_0)))

        if agn is not None:
            nn += numpy.where(numpy.isnan(agn), 0.0, numpy.power(10, -0.4*(agn - mm_0)))

        if baselineType == numpy.ndarray:
            # according to this link
            # http://stackoverflow.com/questions/25087769/runtimewarning-divide-by-zero-error-how-to-avoid-python-numpy
            # we will still get a divide by zero error from log10, but numpy.where will be
            # circumventing the offending value, so it is probably okay
            return numpy.where(nn>tol, -2.5*numpy.log10(nn) + mm_0, numpy.NaN)
        else:
            if nn>tol:
                return -2.5*numpy.log10(nn) + mm_0
            else:
                return numpy.NaN


    def _magnitudeGetter(self, componentName, bandpassDict,
                         columnNameList, indices=None):
        """
        A generic getter for magnitudes of galaxy components.

        @param [in] componentName is either 'bulge', 'disk', or 'agn'

        @param [in] bandpassDict is a BandpassDict of the bandpasses
        in which to calculate the magnitudes

        @param [in] columnNameList is a list of the columns corresponding to
        these magnitudes (for purposes of applying variability).

        @param [in] indices is an optional list of ints corresponding to
        the index in bandpassDict of the magnitudes that will actually
        be output by the InstanceCatalog

        @param [out] magnitudes is a 2-D numpy array of magnitudes in which
        rows correspond to bandpasses and columns correspond to astronomical
        objects.
        """

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
            raise RuntimeError('_magnitudeGetter does not understand component %s ' \
                               % componentName)

        if sedList is None:
            magnitudes = numpy.ones((len(columnNameList), 0))
        else:
            magnitudes = bandpassDict.magListForSedList(sedList, indices=indices).transpose()

        if self._hasCosmoDistMod():
            cosmoDistMod = self.column_by_name('cosmologicalDistanceModulus')
            if len(cosmoDistMod)>0:
               for ix in range(magnitudes.shape[0]):
                   magnitudes[ix] += cosmoDistMod

        for ix, columnName in enumerate(columnNameList):
            if indices is None or ix in indices:
                delta_name = 'delta_' + columnName
                if delta_name in self._all_available_columns:
                    delta = self.column_by_name(delta_name)
                    magnitudes[ix] += delta

        return magnitudes


    @compound('sigma_lsst_u','sigma_lsst_g','sigma_lsst_r',
              'sigma_lsst_i','sigma_lsst_z','sigma_lsst_y')
    def get_photometric_uncertainties_total(self):
        """
        Getter for total photometric uncertainties associated with galaxies
        """
        magnitudes = numpy.array([self.column_by_name('lsst_u'),
                                  self.column_by_name('lsst_g'),
                                  self.column_by_name('lsst_r'),
                                  self.column_by_name('lsst_i'),
                                  self.column_by_name('lsst_z'),
                                  self.column_by_name('lsst_y')])

        return self.calculateMagnitudeUncertainty(magnitudes, self.lsstBandpassDict,
                                                  obs_metadata=self.obs_metadata)

    @compound('sigma_uBulge', 'sigma_gBulge', 'sigma_rBulge',
              'sigma_iBulge', 'sigma_zBulge', 'sigma_yBulge')
    def get_photometric_uncertainties_bulge(self):
        """
        Getter for photometric uncertainties associated with galaxy bulges
        """
        magnitudes = numpy.array([self.column_by_name('uBulge'),
                                  self.column_by_name('gBulge'),
                                  self.column_by_name('rBulge'),
                                  self.column_by_name('iBulge'),
                                  self.column_by_name('zBulge'),
                                  self.column_by_name('yBulge')])

        return self.calculateMagnitudeUncertainty(magnitudes, self.lsstBandpassDict,
                                                  obs_metadata=self.obs_metadata)

    @compound('sigma_uDisk', 'sigma_gDisk', 'sigma_rDisk',
              'sigma_iDisk', 'sigma_zDisk', 'sigma_yDisk')
    def get_photometric_uncertainties_disk(self):
        """
        Getter for photometeric uncertainties associated with galaxy disks
        """
        magnitudes = numpy.array([self.column_by_name('uDisk'),
                                  self.column_by_name('gDisk'),
                                  self.column_by_name('rDisk'),
                                  self.column_by_name('iDisk'),
                                  self.column_by_name('zDisk'),
                                  self.column_by_name('yDisk')])

        return self.calculateMagnitudeUncertainty(magnitudes, self.lsstBandpassDict,
                                                  obs_metadata=self.obs_metadata)


    @compound('sigma_uAgn', 'sigma_gAgn', 'sigma_rAgn',
              'sigma_iAgn', 'sigma_zAgn', 'sigma_yAgn')
    def get_photometric_uncertainties_agn(self):
        """
        Getter for photometric uncertainties associated with Agn
        """
        magnitudes = numpy.array([self.column_by_name('uAgn'),
                                  self.column_by_name('gAgn'),
                                  self.column_by_name('rAgn'),
                                  self.column_by_name('iAgn'),
                                  self.column_by_name('zAgn'),
                                  self.column_by_name('yAgn')])

        return self.calculateMagnitudeUncertainty(magnitudes, self.lsstBandpassDict,
                                                  obs_metadata=self.obs_metadata)



    @compound('uBulge', 'gBulge', 'rBulge', 'iBulge', 'zBulge', 'yBulge')
    def get_lsst_bulge_mags(self):
        """
        Getter for bulge magnitudes in LSST bandpasses
        """

        # load a BandpassDict of LSST bandpasses, if not done already
        if not hasattr(self, 'lsstBandpassDict'):
            self.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()

        # totally optional; list the indices in self.lsstBandpassDict
        # of the columns which will actually be output (so that we don't
        # spend time calculating things we do not need)
        indices = [ii for ii, name in enumerate(self.get_lsst_bulge_mags._colnames) \
                   if name in self._actually_calculated_columns]

        if len(indices)==6:
            indices = None

        # actually calculate the magnitudes
        return self._magnitudeGetter('bulge', self.lsstBandpassDict,
                                     self.get_lsst_bulge_mags._colnames,
                                     indices=indices)



    @compound('uDisk', 'gDisk', 'rDisk', 'iDisk', 'zDisk', 'yDisk')
    def get_lsst_disk_mags(self):
        """
        Getter for galaxy disk magnitudes in the LSST bandpasses
        """

        # load a BandpassDict of LSST bandpasses, if not done already
        if not hasattr(self, 'lsstBandpassDict'):
            self.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()

        # totally optional; list the indices in self.lsstBandpassDict
        # of the columns which will actually be output (so that we don't
        # spend time calculating things we do not need)
        indices = [ii for ii, name in enumerate(self.get_lsst_disk_mags._colnames) \
                   if name in self._actually_calculated_columns]

        if len(indices)==6:
            indices = None

        # actually calculate the magnitudes
        return self._magnitudeGetter('disk', self.lsstBandpassDict,
                                     self.get_lsst_disk_mags._colnames,
                                     indices=indices)


    @compound('uAgn', 'gAgn', 'rAgn', 'iAgn', 'zAgn', 'yAgn')
    def get_lsst_agn_mags(self):
        """
        Getter for AGN magnitudes in the LSST bandpasses
        """

        # load a BandpassDict of LSST bandpasses, if not done already
        if not hasattr(self, 'lsstBandpassDict'):
            self.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()

        # totally optional; list the indices in self.lsstBandpassDict
        # of the columns which will actually be output (so that we don't
        # spend time calculating things we do not need)
        indices = [ii for ii, name in enumerate(self.get_lsst_agn_mags._colnames) \
                   if name in self._actually_calculated_columns]

        if len(indices)==6:
            indices = None

        # actually calculate the magnitudes
        return self._magnitudeGetter('agn', self.lsstBandpassDict,
                                     self.get_lsst_agn_mags._colnames,
                                     indices=indices)


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
                sub_list = [numpy.NaN]*numObj
            else:
                bandpass = columnName[-1]
                bulge = self.column_by_name('%sBulge' % bandpass)
                disk = self.column_by_name('%sDisk' % bandpass)
                agn = self.column_by_name('%sAgn' % bandpass)
                sub_list = self.sum_magnitudes(bulge=bulge, disk=disk, agn=agn)

            output.append(sub_list)
        return numpy.array(output)




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
            return numpy.ones((0))

        if not hasattr(self, '_sedList'):
            self._sedList = SedList(sedNameList, magNormList,
                                         galacticAvList=galacticAvList,
                                         wavelenMatch=wavelen_match)
        else:
            self._sedList.flush()
            self._sedList.loadSedsFromList(sedNameList, magNormList,
                                          galacticAvList=galacticAvList)


    def _magnitudeGetter(self, bandpassDict, columnNameList, indices=None):
        """
        This method gets the magnitudes for an InstanceCatalog, returning them
        in a 2-D numpy array in which rows correspond to bandpasses and columns
        correspond to astronomical objects.

        @param [in] bandpassDict is a BandpassDict containing the bandpasses
        whose magnitudes are to be calculated

        @param [in] columnNameList is a list of the names of the magnitude columns
        being calculated (so that any variability model can determine if a
        'delta_column_name' exists

        @param [in] indices is an optional list of indices indicating which bandpasses to actually
        calculate magnitudes for

        @param [out] magnitudes is a 2-D numpy array of magnitudes in which
        rows correspond to bandpasses in bandpassDict and columns correspond
        to astronomical objects.
        """

        self._loadSedList(bandpassDict.wavelenMatch)

        if not hasattr(self, '_sedList'):
            magnitudes = numpy.ones((len(columnNameList),0))
        else:
            magnitudes = bandpassDict.magListForSedList(self._sedList, indices=indices).transpose()

        for ix, columnName in enumerate(columnNameList):
            if indices is None or ix in indices:
                delta_name = 'delta_' + columnName
                if delta_name in self._all_available_columns:
                    delta = self.column_by_name(delta_name)
                    magnitudes[ix] += delta

        return magnitudes


    @compound('sigma_lsst_u','sigma_lsst_g','sigma_lsst_r','sigma_lsst_i',
              'sigma_lsst_z','sigma_lsst_y')
    def get_lsst_photometric_uncertainties(self):
        """
        Getter for photometric uncertainties associated with stellar
        magnitudes
        """

        magnitudes = numpy.array([self.column_by_name('lsst_u'),
                                  self.column_by_name('lsst_g'),
                                  self.column_by_name('lsst_r'),
                                  self.column_by_name('lsst_i'),
                                  self.column_by_name('lsst_z'),
                                  self.column_by_name('lsst_y')])

        return self.calculateMagnitudeUncertainty(magnitudes, self.lsstBandpassDict,
                                                  obs_metadata=self.obs_metadata)


    @compound('lsst_u','lsst_g','lsst_r','lsst_i','lsst_z','lsst_y')
    def get_lsst_magnitudes(self):
        """
        getter for LSST stellar magnitudes
        """
        if not hasattr(self, 'lsstBandpassDict'):
            self.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()

        indices = [ii for ii, name in enumerate(self.get_lsst_magnitudes._colnames) \
                   if name in self._actually_calculated_columns]

        if len(indices) == 6:
            indices = None

        return self._magnitudeGetter(self.lsstBandpassDict, self.get_lsst_magnitudes._colnames, indices=indices)


class PhotometrySSM(PhotometryBase):
    """
    A mixin to calculate photometry for solar system objects.
    """
    # Because solar system objects will not have dust extinctions, we should be able to read in every
    # SED exactly once, calculate the colors and magnitudes, and then get actual magnitudes by adding
    # an offset based on magNorm.

    def _magnitudeGetter(self, bandpassDict, bandpassTag='lsst'):
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

        @param [in] bandpassTag (optional) is a string indicating the name of the bandpass system
        (i.e. 'lsst', 'sdss', etc.).  This is in case the user wants to calculate the magnitudes
        in multiple systems simultaneously.  In that case, the dict will store magnitudes for each
        Sed in each magnitude system separately.

        @param [out] a numpy array of magnitudes corresponding to bandpassDict.
        """

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
            return numpy.zeros((len(bandpassDict.keys()),0))

        magListOut = []

        for sedName, magNorm in zip(sedNameList, magNormList):
            magTag = bandpassTag+'_'+sedName
            if sedName not in self._ssmMagNormDict or magTag not in self._ssmMagDict:
                dummySed = Sed()
                dummySed.readSED_flambda(os.path.join(self._file_dir, self._spec_map[sedName]))
                fnorm = dummySed.calcFluxNorm(magNorm, self._normalizing_bandpass)
                dummySed.multiplyFluxNorm(fnorm)
                magList = bandpassDict.magListForSed(dummySed)
                self._ssmMagDict[magTag] = magList
                self._ssmMagNormDict[sedName] = magNorm
            else:
                dmag = magNorm - self._ssmMagNormDict[sedName]
                magList = self._ssmMagDict[magTag] + dmag
            magListOut.append(magList)

        return numpy.array(magListOut).transpose()


    @compound('lsst_u','lsst_g','lsst_r','lsst_i','lsst_z','lsst_y')
    def get_lsst_magnitudes(self):
        """
        getter for LSST magnitudes of solar system objects
        """

        if not hasattr(self, 'lsstBandpassDict'):
            self.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()

        return self._magnitudeGetter(self.lsstBandpassDict)


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
            valueList = self.obs_metadata.seeing.values()
            for ix in range(1, len(valueList)):
                if numpy.abs(valueList[ix]-valueList[0])>0.0001:

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
            return numpy.zeros((2,0))

        a_trail = 0.76
        b_trail = 1.16
        a_det = 0.42
        b_det = 0.00
        seeing = self.obs_metadata.seeing[self.obs_metadata.bandpass] # this will be in arcsec
        texp = self.photParams.nexp*self.photParams.exptime
        velocity = numpy.sqrt(numpy.power(numpy.degrees(dradt),2) + numpy.power(numpy.degrees(ddecdt),2))
        x = velocity*texp/(24.0*seeing)
        xsq = numpy.power(x,2)
        dmagTrail = 1.25*numpy.log10(1.0 + a_trail * xsq/(1.0+b_trail*x))
        dmagDetect = 1.25*numpy.log10(1.0 + a_det * xsq/(1.0 + b_det*x))

        return numpy.array([dmagTrail, dmagDetect])
