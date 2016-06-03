from __future__ import print_function
import numpy as np
import copy

from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from  lsst.sims.catUtils.mixins import PhotometryStars, VariabilityStars
from lsst.sims.catalogs.measures.instance import InstanceCatalog, compound
from lsst.sims.utils import haversine

import time

__all__ = ["StellarLightCurveGenerator"]

_sed_cache = {} # a global cache to store SedLists loaded by the light curve catalogs


class _baseLightCurveCatalog(InstanceCatalog):

    column_outputs = ["uniqueId", "raJ2000", "decJ2000",
                      "lightCurveMag", "sigma_lightCurveMag"]

    def iter_catalog(self, chunk_size=None, query_cache=None):
        """
        chunk_size (optional) is an int specifying the number of rows to return
        from the database at a time

        query_cache (optional) is the result of calling db_obj.query_columns().
        DO NOT use this unless you know what you are doing.  It is an optional
        input for those who want to repeatedly examine the same patch of sky
        without actually querying the database over and over again.  If it is set
        to 'None' (default), this method will handle the database query.

        Returns an iterator over rows of the catalog.
        """

        if query_cache is None:
            for line in InstanceCatalog.iter_catalog(self, chunk_size=chunk_size):
                yield line
        else:
            for chunk in query_cache:
                self._set_current_chunk(chunk)
                chunk_cols = [self.transformations[col](self.column_by_name(col))
                              if col in self.transformations.keys() else
                              self.column_by_name(col)
                              for col in self.iter_column_names()]
                for line in zip(*chunk_cols):
                    yield line




class _stellarLightCurveCatalog(_baseLightCurveCatalog, VariabilityStars, PhotometryStars):
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

    def _loadSedList(self, wavelen_match):
        """
        Wraps the PhotometryStars._loadSedList method.

        If current chunk of objects is not represetned in the global
        _sed_cache, this will call the base method defined in
        PhotometryStars.

        Otherwise, it will read self._sedList from the cache.
        That way, the photometry getters defined in PhotometryStars will
        not read in SEDs that have already been cached.
        """

        global _sed_cache

        object_names = self.column_by_name("uniqueId")

        if len(object_names)>0:
            cache_name = "%s_%s" % (object_names[0], object_names[-1])
        else:
            cache_name = None

        if cache_name not in _sed_cache:
            PhotometryStars._loadSedList(self, wavelen_match)

            if cache_name is not None:
                _sed_cache[cache_name] = copy.copy(self._sedList)
        else:
            self._sedList = _sed_cache[cache_name]


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

        return np.array([self.column_by_name("lsst_%s" % self.obs_metadata.bandpass),
                         self.column_by_name("sigma_lsst_%s" % self.obs_metadata.bandpass)])


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


    def generate_light_curves(self, ra, dec, bandpass, expMJD=None, chunk_size=10000):
        """
        Generate light curves for all of the objects in a particular region
        of sky in a particular bandpass.

        Input parameters:
        -----------------
        ra is a tuple indicating the (min, max) values of RA in degrees.

        dec is a tuple indicating the (min, max) values of Dec in degrees.

        bandpass is a char (i.e. 'u', 'g', 'r', etc.) indicating which filter
        you want the light curves in.

        expMJD is an optional tuple indicating a (min, max) range in MJD.
        Defaults to None, in which case, the light curves over the entire
        10 year survey are returned.

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

        global _sed_cache

        # First get the list of ObservationMetaData objects corresponding
        # to the OpSim pointings in the region and bandpass of interest
        obs_list = self._generator.getObservationMetaData(
                                     fieldRA=ra,
                                     fieldDec=dec,
                                     telescopeFilter=bandpass,
                                     expMJD=expMJD,
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

                for star_obj in \
                cat.iter_catalog(chunk_size=chunk_size, query_cache=dataCache):

                    if star_obj[0] not in mjd_dict:
                        mjd_dict[star_obj[0]] = []
                        mag_dict[star_obj[0]] = []
                        sig_dict[star_obj[0]] = []

                    mjd_dict[star_obj[0]].append(obs.mjd.TAI)
                    mag_dict[star_obj[0]].append(star_obj[3])
                    sig_dict[star_obj[0]].append(star_obj[4])

            _sed_cache = {} # reset sed cache before moving onto the next group of
                            # ObservationMetaData

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
