from __future__ import print_function
import numpy as np
import copy

from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from  lsst.sims.catUtils.mixins import PhotometryStars, VariabilityStars
from lsst.sims.catalogs.measures.instance import InstanceCatalog, compound
from lsst.sims.utils import haversine

import time

__all__ = ["StellarLightCurveGenerator"]

class _stellarLightCurveCatalog(InstanceCatalog, VariabilityStars, PhotometryStars):
    """
    This class wraps a basic stellar variability InstanceCatalog.  It provides its
    own photometry getter that

    1) only returns the magnitude and uncertainty in the bandpass specified by
    self.obs_metadata

    2) caches all of the SEDs read in so that they can be reused when sampling the
    objects in this catalog at a different MJD.

    It provides its own iter_catalog() method that allows users to pass the cached
    results of database queries in so that the catalog does not have to query the
    database multiple times to get the same result.

    It should only be used in the context of the LightCurveGenerator class.
    """

    column_outputs = ["uniqueId", "raJ2000", "decJ2000",
                      "lightCurveMag", "sigma_lightCurveMag"]

    _sedList_cache = None # this will become a list of the SedList objects read in by the getter

    _sedList_to_use = None # this will be filled in by an SedList instantiation to be used by
                           # the getter

    def _reset_sed_cache(self):
        """
        Reset the caches used to store SEDs across iterations on the same
        patch of sky.
        """
        self._sedList_cache = None
        self._sedList_to_use = None


    def _loadSedList(self, wavelen_match):
        """
        Wraps the PhotometryStars._loadSedList method.

        If self._sedList_to_use is None, this will call the base method.
        defined in PhotometryStars.

        Otherwise, it will set self._sedList=self._sedList_to_use.
        That way, the photometry getters defined in PhotometryStars will
        not read in SEDs that have already been cached.
        """
        if self._sedList_to_use is None:
            PhotometryStars._loadSedList(self, wavelen_match)
        else:
            self._sedList = self._sedList_to_use


    def iter_catalog(self, chunk_size=None, query_cache=None, sed_cache=None):
        """
        chunk_size (optional) is an int specifying the number of rows to return
        from the database at a time

        query_cache (optional) is the result of calling db_obj.query_columns().
        DO NOT use this unless you know what you are doing.  It is an optional
        input for those who want to repeatedly examine the same patch of sky
        without actually querying the database over and over again.  If it is set
        to 'None' (default), this method will handle the database query.

        sed_cache (optional) is a list of SedLists corresponding to the data
        chunks in query_cache.  Passing this in will cause the photometry
        getters in this InstanceCatalog class to use the cached SedLists, rather
        than reading in new ones when calculating magnitudes and uncertainties.

        Returns an iterator over rows of the catalog.
        """

        if query_cache is None:
            yield InstanceCatalog.iter_catalog(self)

        if sed_cache is None:
            sed_cache = [None]*len(query_cache)

        for chunk, sed_list in zip(query_cache, sed_cache):
            self._set_current_chunk(chunk)
            self._sedList_to_use = sed_list
            chunk_cols = [self.transformations[col](self.column_by_name(col))
                          if col in self.transformations.keys() else
                          self.column_by_name(col)
                          for col in self.iter_column_names()]
            for line in zip(*chunk_cols):
                yield line


    @compound("lightCurveMag", "sigma_lightCurveMag")
    def get_lightCurvePhotometry(self):
        """
        A getter which returns the magnitudes and uncertainties in magnitudes
        in the bandpass specified by self.obs_metdata.

        As it runs, this method will cache the SedLists it reads in so that
        they can be used later.
        """

        if len(self.obs_metadata.bandpass) != 1:
            raise RuntimeError("_stellarLightCurveCatalog cannot handle bandpass "
                               "%s" % str(self.obs_metadata.bandpass))

        mag = self.column_by_name("lsst_%s" % self.obs_metadata.bandpass)
        sigma = self.column_by_name("sigma_lsst_%s" % self.obs_metadata.bandpass)

        if self._sedList_cache is None and len(mag)>0:
            self._sedList_cache = []

        if self._sedList_to_use is None and len(mag)>0:
            self._sedList_cache.append(copy.copy(self._sedList))

        return np.array([mag, sigma])


class LightCurveGenerator(object):
    """
    This class will find all of the OpSim pointings in a particular region
    of the sky in a particular filter and then return light curves for all
    of the objects observed in that region of sky.

    Input parameters:
    -----------------
    catalogdb is a CatalogDBObject instantiation connecting to the database
    of objects to be observed.

    opsimdb is the path to the OpSim database of observation.

    opsimdriver (optional; default 'sqlite') indicates the database driver to
    be used when connecting to opsimdb.
    """

    _lightCurveCatalogClass = None

    def __init__(self, catalogdb, opsimdb, opsimdriver="sqlite"):
        self._generator = ObservationMetaDataGenerator(database=opsimdb,
                                                       driver=opsimdriver)

        self._catalogdb = catalogdb


    def generate_light_curves(self, ra, dec, bandpass, chunk_size=10000):
        """
        Generate light curves for all of the objects in a particular region
        of sky in a particular bandpass.

        Input parameters:
        -----------------
        ra is a tuple indicating the (min, max) values of RA in degrees.

        dec is a tuple indicating the (min, max) values of Dec in degrees.

        bandpass is a char (i.e. 'u', 'g', 'r', etc.) indicating which filter
        you want the light curves in.

        chunk_size (optional; default=10000) is an int specifying how many
        objects to pull in from the database at a time.  Note: the larger
        this is, the faster the LightCurveGenerator will run, because it
        will be handling more objects in memory at once.

        Output:
        -------
        A dict keyed on the object's uniqeId (an int).  Each entry in the dict
        is a 2-D numpy array.  The first row is a numpy array of MJDs.  The
        second row is a numpy array of magnitudes.  The third row is a numpy
        array of magnitude uncertainties.
        """

        # First get the list of ObservationMetaData objects corresponding
        # to the OpSim pointings in the region and bandpass of interest
        obs_list = self._generator.getObservationMetaData(
                                     fieldRA=ra,
                                     fieldDec=dec,
                                     telescopeFilter=bandpass,
                                     boundLength=1.75)

        mjd_dict = {}
        mag_dict = {}
        sig_dict = {}

        if len(obs_list) == 0:
            print("No observations found matching your criterion")
            return None

        t_start = time.clock()
        print('starting light curve generation')

        # Group the OpSim pointings so that all of the pointings centered on the same
        # point in the sky are in a list together (this will allow us to generate the
        # light curves one pointing at a time without having to query the database for
        # the same results more than once.
        tol = 1.0e-12

        obs_groups = [] # a list of list of the indices of the ObservationMetaDatas
                        # in obs_list.  All of the ObservationMetaData in
                        # obs_list[i] will point to the same point on the sky.

        for iobs, obs in enumerate(obs_list):
            group_dex = -1

            for ix, obs_g in enumerate(obs_groups):
                dd = haversine(obs._pointingRA, obs._pointingDec,
                               obs_list[obs_g[0]]._pointingRA, obs_list[obs_g[0]]._pointingDec)
                if dd<tol:
                    group_dex = ix
                    break

            if group_dex == -1:
                obs_groups.append([iobs])
            else:
                obs_groups[group_dex].append(iobs)


        cat = None # a placeholder for the _stellarLightCurveCatalog

        # Loop over the list of groups ObservationMetaData objects,
        # querying the database and generating light curves.
        for grp in obs_groups:

            dataCache = None # a cache for the results of database queries
            sedCache = None # a cache for SedList objects loaded by the photometry code

            # loop over the group of like-pointed ObservationMetaDatas, generating
            # the light curves
            for ix in grp:
                obs = obs_list[ix]
                if cat is None:
                    cat = self._lightCurveCatalogClass(self._catalogdb, obs_metadata=obs)
                    db_required_columns = cat.db_required_columns()

                # query the database, caching the results.
                if dataCache is None:
                    query_result = cat.db_obj.query_columns(colnames=cat._active_columns,
                                                            obs_metadata=obs,
                                                            constraint=cat.constraint,
                                                            chunk_size=chunk_size)
                    dataCache = []
                    for chunk in query_result:
                       dataCache.append(chunk)

                # reset global variables to reflect the specific observing conditions
                cat.obs_metadata = obs
                cat._gamma_cache = {}

                ct = 0
                for star_obj in \
                cat.iter_catalog(chunk_size=chunk_size, query_cache=dataCache, sed_cache=sedCache):

                    if star_obj[0] not in mjd_dict:
                        mjd_dict[star_obj[0]] = []
                        mag_dict[star_obj[0]] = []
                        sig_dict[star_obj[0]] = []

                    mjd_dict[star_obj[0]].append(obs.mjd.TAI)
                    mag_dict[star_obj[0]].append(star_obj[3])
                    sig_dict[star_obj[0]].append(star_obj[4])


                # If the _stellarLightCurveCatalog produced a cache of SedLists, copy that
                # into SedCache so that you don't have to load them again.
                if cat._sedList_cache is not None:
                    sedCache = cat._sedList_cache

            # clear the SedList caches in the _stellarLightCurveCatalog so that you don't
            # accidentally use the wrong SEDs on the wrong patch of sky
            cat._reset_sed_cache()

        output_dict = {}
        for unique_id in mjd_dict:
            output_dict[unique_id] = np.array([mjd_dict[unique_id],
                                               mag_dict[unique_id],
                                               sig_dict[unique_id]])

        print('that took %e; grps %d' % (time.clock()-t_start, len(obs_groups)))
        print('len obs_list %d' % len(obs_list))
        return output_dict


class StellarLightCurveGenerator(LightCurveGenerator):

    def __init__(self, *args, **kwargs):
        self._lightCurveCatalogClass = _stellarLightCurveCatalog
        super(StellarLightCurveGenerator, self).__init__(*args, **kwargs)
