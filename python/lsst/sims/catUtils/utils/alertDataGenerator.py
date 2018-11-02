import numpy as np
import os
import re
import sqlite3
from collections import OrderedDict
import time
import gc
import multiprocessing.synchronize as mprocSync
from lsst.utils import getPackageDir
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.utils import trixelFromHtmid, getAllTrixels
from lsst.sims.utils import levelFromHtmid, halfSpaceFromRaDec
from lsst.sims.utils import angularSeparation, ObservationMetaData
from lsst.sims.utils import arcsecFromRadians
from lsst.sims.catUtils.utils import _baseLightCurveCatalog
from lsst.sims.utils import _pupilCoordsFromRaDec
from lsst.sims.coordUtils import chipNameFromPupilCoordsLSST
from lsst.sims.coordUtils import pixelCoordsFromPupilCoordsLSST

from lsst.sims.catalogs.decorators import compound, cached
from lsst.sims.photUtils import BandpassDict, Sed, calcSNR_m5
from lsst.sims.photUtils import PhotometricParameters
from lsst.sims.catUtils.mixins import VariabilityStars, AstrometryStars
from lsst.sims.catUtils.mixins import VariabilityGalaxies, AstrometryGalaxies
from lsst.sims.catUtils.mixins import CameraCoordsLSST, PhotometryBase
from lsst.sims.catUtils.mixins import ParametrizedLightCurveMixin
from lsst.sims.catUtils.mixins import create_variability_cache

from lsst.sims.catUtils.baseCatalogModels import StarObj, GalaxyAgnObj
from sqlalchemy.sql import text
from lsst.sims.catalogs.db import ChunkIterator

__all__ = ["AlertDataGenerator",
           "AlertStellarVariabilityCatalog",
           "AlertAgnVariabilityCatalog",
           "_baseAlertCatalog",
           "StellarAlertDBObj",
           "AgnAlertDBObj",
           "StellarAlertDBObjMixin"]


class _lock_context(object):
    """
    A class to handle lock acquisition through a context manager
    which allows for the possiblity that lock is None
    """

    def __init__(self, lock_obj):
        self._lock = lock_obj

    def __enter__(self):
        if isinstance(self._lock, mprocSync.Lock):
            self._lock.acquire()

    def __exit__(self, *args):
        if isinstance(self._lock, mprocSync.Lock):
            self._lock.release()


class StellarAlertDBObjMixin(object):
    """
    Mimics StarObj class, except it allows you to directly query
    all objects in a trixel specified by an htmid.
    """
    def query_columns_htmid(self, colnames=None, chunk_size=None,
                            constraint=None,
                            limit=None, htmid=None):
        """Execute a query from the primary catsim database

        Execute a query, taking advantage of the spherical geometry library and
        htmid indexes on all catalog tables in the UW catsim database

        **Parameters**

            * colnames : list or None
              a list of valid column names, corresponding to entries in the
              `columns` class attribute.  If not specified, all columns are
              queried.
            * chunk_size : int (optional)
              if specified, then return an iterator object to query the database,
              each time returning the next `chunk_size` elements.  If not
              specified, all matching results will be returned.
            * constraint : str (optional)
              a string which is interpreted as SQL and used as a predicate on the query
            * limit : int (optional)
              limits the number of rows returned by the query
            * htmid is the htmid to be queried

        **Returns**

            * result : list or iterator
              If chunk_size is not specified, then result is a list of all
              items which match the specified query.  If chunk_size is specified,
              then result is an iterator over lists of the given size.
        """

        # find the minimum and maximum htmid
        # (level=21 since that is what is implemented
        # on fatboy) that we are asking for
        #
        # Note that sqlalchemy does not like np.int64
        # as a data type
        current_level = levelFromHtmid(htmid)
        n_bits_off = 2*(21-current_level)
        htmid_min = int(htmid << n_bits_off)
        htmid_max = int((htmid+1) << n_bits_off)

        query = self._get_column_query(colnames)

        # add spatial constraints to query.

        # Hint sql engine to seek on htmid
        if not self.tableid.endswith('forceseek'):
            query = query.with_hint(self.table, ' WITH(FORCESEEK)', 'mssql')

        # SQL is not case sensitive but python is:
        if 'htmID' in self.columnMap:
            htmid_name = 'htmID'
        elif 'htmid' in self.columnMap:
            htmid_name = 'htmid'
        else:
            htmid_name = 'htmId'

        # Range join on htmid ranges
        query = query.filter(self.table.c[htmid_name].between(htmid_min, htmid_max))

        if constraint is not None:
            query = query.filter(text(constraint))

        if limit is not None:
            query = query.limit(limit)

        return ChunkIterator(self, query, chunk_size)


class StellarAlertDBObj(StellarAlertDBObjMixin, StarObj):
    pass


class AgnAlertDBObj(GalaxyAgnObj):
    """
    Mimics GalaxyAgnObj class, except it allows you to directly query
    all objects in a trixel specified by an htmid.
    """

    columns = [('htmid', 0, np.int64),
               ('galtileid', None, np.int64),
               ('galid', None, str, 30),
               ('componentra', 'agnra*PI()/180.'),
               ('componentdec', 'agndec*PI()/180.'),
               #: This is actually a problem with the stored procedure.
               #: We need to be able to map columns other than
               #: just ra/dec to raJ2000/decJ2000.  This gets
               #: important when we start perturbing the three galaxy components
               ('raJ2000', 'ra'),
               ('decJ2000', 'dec'),
               ('magNorm', 'magnorm_agn'),
               ('magNormAgn', 'magnorm_agn'),
               ('sedFilename', 'sedname_agn', str, 40),
               ('sedFilenameAgn', 'sedname_agn', str, 40),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('lsst_u', 'u_ab'),
               ('lsst_g', 'g_ab'),
               ('lsst_r', 'r_ab'),
               ('lsst_i', 'i_ab'),
               ('lsst_z', 'z_ab'),
               ('lsst_y', 'y_ab')]

    def query_columns_htmid(self, colnames=None, chunk_size=None,
                            constraint=None,
                            limit=None, htmid=None):
        """Execute a query from the primary catsim database

        Execute a query, taking advantage of the spherical geometry library and
        htmid indexes on all catalog tables in the UW catsim database

        **Parameters**

            * colnames : list or None
              a list of valid column names, corresponding to entries in the
              `columns` class attribute.  If not specified, all columns are
              queried.
            * chunk_size : int (optional)
              if specified, then return an iterator object to query the database,
              each time returning the next `chunk_size` elements.  If not
              specified, all matching results will be returned.
            * constraint : str (optional)
              a string which is interpreted as SQL and used as a predicate on the query
            * limit : int (optional)
              limits the number of rows returned by the query
            * htmid is the htmid to be queried

        **Returns**

            * result : list or iterator
              If chunk_size is not specified, then result is a list of all
              items which match the specified query.  If chunk_size is specified,
              then result is an iterator over lists of the given size.
        """

        trixel = trixelFromHtmid(htmid)
        ra_0, dec_0 = trixel.get_center()
        new_obs = ObservationMetaData(pointingRA=ra_0, pointingDec=dec_0, boundType='circle',
                                      boundLength=trixel.get_radius()+0.1)

        self._queried_trixel = trixel
        self._queried_htmid_level = levelFromHtmid(htmid)

        return self.query_columns(colnames=colnames, chunk_size=chunk_size,
                                  obs_metadata=new_obs, constraint=constraint,
                                  limit=limit)

    def _final_pass(self, results):
        """Modify the results of raJ2000 and decJ2000 to be in radians.

        Also filter the results so that any objects outside of the
        trixel specified in query_columns_htmid are returned with
        htmid=0.

        **Parameters**

            * results : Structured array of results from query

        **Returns**

            * results : Modified structured array

        """

        if hasattr(self, '_queried_trixel'):
            htmid = self._queried_trixel.htmid
            htmid_21 = htmid << 2*(21-self._queried_htmid_level)
            assert levelFromHtmid(htmid_21) == 21
            contains_arr = self._queried_trixel.contains(results['raJ2000'], results['decJ2000'])
            results['htmid'] = np.where(contains_arr, htmid_21, 0)

        results['raJ2000'] = np.radians(results['raJ2000'])
        results['decJ2000'] = np.radians(results['decJ2000'])
        return results


class _baseAlertCatalog(PhotometryBase, CameraCoordsLSST, _baseLightCurveCatalog):

    column_outputs = ['htmid', 'uniqueId', 'raICRS', 'decICRS',
                      'flux', 'SNR', 'dflux',
                      'chipNum', 'xPix', 'yPix']

    default_formats = {'f': '%.4g'}

    default_columns = [('properMotionRa', 0.0, float),
                       ('properMotionDec', 0.0, float),
                       ('parallax', 0.0, float)]

    def iter_catalog_chunks(self, chunk_size=None, query_cache=None, column_cache=None):
        """
        Returns an iterator over chunks of the catalog.

        Parameters
        ----------
        chunk_size : int, optional, defaults to None
            the number of rows to return from the database at a time. If None,
            returns the entire database query in one chunk.

        query_cache : iterator over database rows, optional, defaults to None
            the result of calling db_obj.query_columns().  If query_cache is not
            None, this method will iterate over the rows in query_cache and produce
            an appropriate InstanceCatalog. DO NOT set to non-None values
            unless you know what you are doing.  It is an optional
            input for those who want to repeatedly examine the same patch of sky
            without actually querying the database over and over again.  If it is set
            to None (default), this method will handle the database query.

        column_cache : a dict that will be copied over into the catalogs self._column_cache.
            Should be left as None, unless you know what you are doing.
        """

        if query_cache is None:
            # Call the original version of iter_catalog defined in the
            # InstanceCatalog class.  This version of iter_catalog includes
            # the call to self.db_obj.query_columns, which the user would have
            # used to generate query_cache.
            for line in InstanceCatalog.iter_catalog_chunks(self, chunk_size=chunk_size):
                yield line
        else:
            # Otherwise iterate over the query cache
            transform_keys = list(self.transformations.keys())
            for chunk in query_cache:
                self._set_current_chunk(chunk, column_cache=column_cache)
                chunk_cols = [self.transformations[col](self.column_by_name(col))
                              if col in transform_keys else
                              self.column_by_name(col)
                              for col in self.iter_column_names()]

                if not hasattr(self, '_chunkColMap_output'):

                    self._chunkColMap_output = dict([(col, i)
                                                     for i, col in
                                                     enumerate(self.iter_column_names())])

                yield chunk_cols, self._chunkColMap_output

        self._column_cache = {}
        self._current_chunk = None

    @cached
    def get_chipName(self):
        if len(self.column_by_name('uniqueId')) == 0:
            return np.array([])
        raise RuntimeError("Should not get this far in get_chipName")

    @compound('x_pupil', 'y_pupil')
    def get_pupilFromSky(self):
        if len(self.column_by_name('uniqueId')) == 0:
            return np.array([[], []])
        raise RuntimeError("Should not get this far in get_pupilFromSky")

    @cached
    def get_chipNum(self):
        """
        Concatenate the digits in 'R:i,j S:m,n' to make the chip number ijmn
        """
        chip_name = self.column_by_name('chipName')
        return np.array([int(''.join(re.findall(r'\d+', name))) if name is not None else 0
                        for name in chip_name])

    @compound('xPix', 'yPix')
    def get_pixelCoordinates(self):
        xpup = self.column_by_name('x_pupil')
        ypup = self.column_by_name('y_pupil')
        chip_name = self.column_by_name('chipName')
        xpix, ypix = pixelCoordsFromPupilCoordsLSST(xpup, ypup, chipName=chip_name,
                                                    band=self.obs_metadata.bandpass,
                                                    includeDistortion=True)
        return np.array([xpix, ypix])

    @compound('delta_umag', 'delta_gmag', 'delta_rmag',
              'delta_imag', 'delta_zmag', 'delta_ymag')
    def get_deltaMagAvro(self):
        ra = self.column_by_name('raJ2000')
        if len(ra) == 0:
            return np.array([[], [], [], [], [], []])

        raise RuntimeError("Should not have gotten this far in delta mag getter")

    @compound('lsst_u', 'lsst_g', 'lsst_r', 'lsst_i', 'lsst_z', 'lsst_y')
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

        delta = np.array([self.column_by_name('delta_umag'),
                          self.column_by_name('delta_gmag'),
                          self.column_by_name('delta_rmag'),
                          self.column_by_name('delta_imag'),
                          self.column_by_name('delta_zmag'),
                          self.column_by_name('delta_ymag')])
        magnitudes += delta

        return magnitudes

    @compound('mag', 'dmag', 'quiescent_mag')
    def get_alertPhotometry(self):
        mag = self.column_by_name('lsst_%s' % self.obs_metadata.bandpass)
        quiescent_mag = self.column_by_name('quiescent_lsst_%s' % self.obs_metadata.bandpass)
        dmag = mag - quiescent_mag

        return np.array([mag, dmag, quiescent_mag])

    @cached
    def get_flux(self):
        if len(self.column_by_name('uniqueId')) == 0:
            return np.array([])
        raise RuntimeError("Should not get this far in get_flux")

    @cached
    def get_dflux(self):
        if len(self.column_by_name('uniqueId')) == 0:
            return np.array([])
        raise RuntimeError("Should not get this far in get_dflux")

    @cached
    def get_SNR(self):
        if len(self.column_by_name('uniqueId')) == 0:
            return np.array([])
        raise RuntimeError("Should not get this far in get_SNR")


class AlertStellarVariabilityCatalog(_baseAlertCatalog,
                                     VariabilityStars,
                                     AstrometryStars):

    @compound('quiescent_lsst_u', 'quiescent_lsst_g', 'quiescent_lsst_r',
              'quiescent_lsst_i', 'quiescent_lsst_z', 'quiescent_lsst_y')
    def get_quiescent_lsst_magnitudes(self):
        return np.array([self.column_by_name('umag'), self.column_by_name('gmag'),
                         self.column_by_name('rmag'), self.column_by_name('imag'),
                         self.column_by_name('zmag'), self.column_by_name('ymag')])


class AlertAgnVariabilityCatalog(_baseAlertCatalog,
                                 VariabilityGalaxies,
                                 AstrometryGalaxies):

    @compound('quiescent_lsst_u', 'quiescent_lsst_g', 'quiescent_lsst_r',
              'quiescent_lsst_i', 'quiescent_lsst_z', 'quiescent_lsst_y')
    def get_quiescent_lsst_magnitudes(self):
        return np.array([self.column_by_name('u_ab'), self.column_by_name('g_ab'),
                         self.column_by_name('r_ab'), self.column_by_name('i_ab'),
                         self.column_by_name('z_ab'), self.column_by_name('y_ab')])


class AlertDataGenerator(object):
    """
    This class will read in astrophysical sources and variability
    models from CatSim, observe them with a simulated OpSim
    cadence, and write a series of sqlite files containing all
    of the simulated observations that could trigger an alert.

    In order to make this calculation as efficient as possible,
    the class works by partitioning the sky according to the
    Hierarchical Triangular Mesh (HTM) of

    Kunszt P., Szalay A., Thakar A. (2006) in "Mining The Sky",
    Banday A, Zaroubi S, Bartelmann M. eds.
    ESO Astrophysics Symposia
    https://www.researchgate.net/publication/226072008_The_Hierarchical_Triangular_Mesh

    Szalay A. et al. (2007)
    "Indexing the Sphere with the Hierarchical Triangular Mesh"
    arXiv:cs/0701164

    and simulating all of the observations in a given trixel (the
    elementary unit of the HTM) at once.  Accordintly, the outputs
    produced by this class are files named like

    prefix_NNNN_sqlite.db

    where prefix is specified by theuser and NNNN is the htmid, the
    unique identifying integer, corresponding to each simulated trixel.

    The proper way to run this class is to instantiate it, run
    subdivide_obs on a list of ObservationMetaData corresponding
    to the OpSim pointings to be simulated, and then running
    alert_data_from_htmid on each of the htmid in the class property
    htmid_list.  This last step can easily be parallelized using python's
    multiprocessing module, with each process handling a different htmid.

    The sqlite files produced by alert_data_from_htmid will each contain
    four tables.  They are as follows.  Columns are listed below the
    tables.

    alert_data
    ----------
        uniqueId -- int -- a unique identifier for each astrophysical object

        obshistId -- int -- a unique identifier for each OpSim pointing

        xPix -- float -- the x pixel coordinate of the source on the focal
        plane

        yPix -- float -- the y pixel coordinate of the source on the focal
        plane

        dflux -- float -- the difference in flux between the source's current
        flux and its quiescent flux (the source's quiescent flux can be found
        in the quiescent_flux table).  This is in units of Janskys.

        snr -- float -- the signal to noise of the current detection of the
        source (not the signal to noise of the source's detection in a simulated
        difference image).

        ra -- float -- the current RA of the source in degrees

        dec -- float -- the current Declination of the source in degrees.

        The alert_data table has a mult-column index on uniqueId and obshistId.

    metadata
    --------
        obshistId -- int --- a unique identifier for each OpSim pointing

        TAI -- float -- the International Atomic Time of the observation
        as an MJD (in days)

        band -- int -- denotes the filter used for the observation
        (0=u, 1=g, 2=r, etc.)

        The metadata table is indexed on obshistId

    quiescent_flux
    --------------
        uniqueId -- int -- a unique identifier for each astrophysical source

        band -- int -- an integer denoting each LSST filter (0=u, 1=g, 2=r, etc.)

        flux -- float -- the flux of the source through the filter specified by
        band (in units of Janskys)

        snr -- float -- the signal to noise ratio of the source in the given
        band with m5 taken from Table 2 of the overview paper (arXiv:0805.2366)

        The quiescent_flux table has a multi-column index on uniqueId and band.

    baseline_astrometry
    -------------------
        uniqueId -- int -- a unique identifier for each astrophysical source

        TAI -- float -- the International Atomic Time of the baseline astrometric
        measurements below as a MJD (in days)

        ra -- float -- the RA of the source at TAI in degrees

        dec -- float -- the Declination of the source at TAI in degrees

        pmRA -- float -- the RA proper motion of the source in milliarcseconds/year

        pmDec -- float -- the Declination proper motion of the source in
        milliarcseconds/year

        parallax -- float -- the parallax of the source in milliarcseconds

        The baseline_astrometry table is indexed on uniqueId

    """
    _output_columns = ('uniqueId', 'raICRS', 'decICRS', 'flux', 'dflux', 'SNR',
                       'chipNum', 'xPix', 'yPix')
    def __init__(self,
                 testing=False):
        """
        Parameters
        ----------
        testing as a boolean that should only be True when running the unit tests.
        If True, it prevents the AlertDataGenerator from pre-caching variability
        models, which aids performance, but uses more memory than we want to use
        in a unit test.
        """

        self.mag_name_to_int = {'u': 0, 'g': 1, 'r': 2,
                                'i': 3, 'z': 4, 'y': 5}

        self.phot_params = PhotometricParameters()

        self._variability_cache = create_variability_cache()
        self._stdout_lock = None
        if not testing:
            plm = ParametrizedLightCurveMixin()
            plm.load_parametrized_light_curves(variability_cache = self._variability_cache)
        self.bp_dict = BandpassDict.loadTotalBandpassesFromFiles()

        # This is a file that lists the maximum amplitude of variability
        # for each of the Kepler-derived light curve models.  It will be
        # used by the stellar variability model to figure out which
        # stars can be skipped because they will never vary above
        # the alert-triggering threshold.
        self._dmag_lookup_file_exists = True
        self._dmag_lookup_file = os.path.join(getPackageDir('sims_data'),
                                              'catUtilsData',
                                              'kplr_dmag_171204.txt')

        if not os.path.exists(self._dmag_lookup_file):
            if not testing:
                script_name = os.path.join(getPackageDir('sims_catUtils'), 'support_scripts',
                                           'get_kepler_dmag.sh')
                raise RuntimeError('\n%s does not exist; run the script\n\n%s\n\n' %
                                   script_name)
            else:
                self._dmag_lookup_file_exists = False

    def subdivide_obs(self, obs_list, htmid_level=6):
        """
        Take a list of ObservationMetaData and subdivide
        them according to which trixels (see htmModule.py
        in sims_utils) they intersect.

        Parameters
        ----------
        obs_list is a list of ObservationMetaData

        htmid_level is an int denoting the level of
        the HTM mesh you want to use to tile the sky
        (higher htmid_level corresponds to a finer
        tiling).  Default = 6

        Returns
        -------
        Nothing.

        After running this method, this AlertGenerator
        will contain the following data.

        - a list of the htmid of every trixel intersected
        by the fields of view specified in obs_list.  This
        list is accessible from the property
        AlertGenerator.htmid_list

        - a dict mapping each htmid to the ObservationMetaData
        from obs_list that intersect it.  The method
        AlertGenerator.obs_from_htmid(htmid) will return a
        list of all of the ObservationMetaData that intersect
        the trixel specified by htmid.
        """
        self._trixel_dict = getAllTrixels(htmid_level)
        valid_htmid = []
        for htmid in self._trixel_dict:
            if levelFromHtmid(htmid) == htmid_level:
                valid_htmid.append(htmid)

        obs_list = np.array(obs_list)
        self._obs_list = obs_list
        obs_ra_list = []
        obs_dec_list = []
        halfspace_list = []
        for obs in obs_list:
            obs_ra_list.append(obs.pointingRA)
            obs_dec_list.append(obs.pointingDec)
            hs = halfSpaceFromRaDec(obs.pointingRA,
                                    obs.pointingDec,
                                    obs.boundLength)
            halfspace_list.append(hs)

        obs_ra_list = np.array(obs_ra_list)
        obs_dec_list = np.array(obs_dec_list)
        halfspace_list = np.array(halfspace_list)

        self._htmid_dict = {}
        self._htmid_list = []
        n_obs_list = []
        fov_radius = 1.75
        for i_htmid, htmid in enumerate(valid_htmid):
            trixel = self._trixel_dict[htmid]
            ra_c, dec_c = trixel.get_center()
            radius = trixel.get_radius()
            obs_distance = angularSeparation(ra_c, dec_c, obs_ra_list, obs_dec_list)
            valid_obs = np.where(obs_distance < radius + fov_radius)
            if len(valid_obs[0]) > 0:
                final_obs_list = []
                for obs_dex in valid_obs[0]:
                    hs = halfspace_list[obs_dex]
                    obs = obs_list[obs_dex]
                    if hs.contains_trixel(trixel) != 'outside':
                        final_obs_list.append(obs_dex)

                if len(final_obs_list) == 0:
                    continue

                self._htmid_dict[htmid] = np.array(final_obs_list)
                self._htmid_list.append(htmid)
                n_obs_list.append(len(final_obs_list))

        n_obs_list = np.array(n_obs_list)
        self._htmid_list = np.array(self._htmid_list)
        sorted_dex = np.argsort(-1.0*n_obs_list)
        self._htmid_list = self._htmid_list[sorted_dex]
        print('done subdividing obs list -- %d htmid' %
              len(self._htmid_list))

    @property
    def htmid_list(self):
        """
        A list of the unique htmids corresponding to the trixels
        that need to be queried to generate the alert data
        """
        return self._htmid_list

    def n_obs(self, htmid):
        """
        Return the number of observations that intersect
        the trixel specified by htmid.

        Must run subdivide_obs in order for this method to
        work.
        """
        return len(self._htmid_dict[htmid])

    def obs_from_htmid(self, htmid):
        """
        Return a numpy array containing all of the ObservationMetaData
        that intersect the trixel specified by htmid.

        Must run subdivide_obs in order for this method to
        work.
        """
        return self._obs_list[self._htmid_dict[htmid]]

    def _output_alert_data(self, conn, data_cache):
        """
        Write a cache of alert data to the sqlite file currently open.

        Parameters
        ----------
        conn is the connection to the sqlite file (already open)

        data_cache is a dict containing all of the data to be written.
        It will keyed on a string like 'i_j' where i is the obshistID
        of an OpSim pointing and j is an arbitrary integer.  That key
        will lead to another dict keyed on the columns being output to
        the sqlite file.  The values of this second layer of dict are
        numpy arrays.

        Returns
        -------
        The number of rows written to the sqlite file
        """

        cursor = conn.cursor()
        n_rows_0 = cursor.execute('SELECT COUNT(uniqueId) FROM alert_data').fetchall()

        chunk_lengths = np.zeros(len(data_cache))

        for i_cache_tag, cache_tag in enumerate(data_cache):
            obshistid = int(cache_tag.split('_')[0])
            n_obj = len(data_cache[cache_tag]['uniqueId'])
            chunk_lengths[i_cache_tag] = n_obj

            values = ((int(data_cache[cache_tag]['uniqueId'][i_obj]),
                      obshistid,
                      data_cache[cache_tag]['xPix'][i_obj],
                      data_cache[cache_tag]['yPix'][i_obj],
                      int(data_cache[cache_tag]['chipNum'][i_obj]),
                      data_cache[cache_tag]['dflux'][i_obj],
                      data_cache[cache_tag]['SNR'][i_obj],
                      np.degrees(data_cache[cache_tag]['raICRS'][i_obj]),
                      np.degrees(data_cache[cache_tag]['decICRS'][i_obj]))
                      for i_obj in range(n_obj))

            cursor.executemany('INSERT INTO alert_data VALUES (?,?,?,?,?,?,?,?,?)', values)

        conn.commit()

        n_rows_1 = cursor.execute('SELECT COUNT(uniqueId) FROM alert_data').fetchall()
        conn.commit()
        n_written = (n_rows_1[0][0]-n_rows_0[0][0])

        return n_written

    def _filter_on_photometry_then_chip_name(self, chunk, column_query,
                                             photometry_catalog,
                                             q_m_dict,
                                             q_snr_dict,
                                             dmag_cutoff,
                                             snr_cutoff):
        """
        Determine which simulated observations are actually worth storing
        by first figuring out which observations of which objects are
        photometrically detectable and alert-worthy, then determining
        which of those actually fall on an LSST detector.

        Parameters
        ----------
        chunk is the output yielded from a CatSim ChunkIterator.  It is
        a numpy recarray representing one chunk_size query from the
        underlying simulations database

        column_query is a list of the columns that were queried from
        the database

        photometry_catalog is an instantiation of the InstanceCatalog class
        being used to calculate magnitudes for these variable sources.

        q_m_dict is a dict keyed on i_filter (u=0, g=1, etc.)
        that returns a numpy array of the quiescent magnitudes
        of all of the sources in the current chunk.

        q_snr_dict is a dict keyed on i_filter (u=0, g=1, etc.)
        that returns a numpy array of the quiescent SNR of all
        of the sources in the current chunk.

        dmag_cutoff is the minimum delta magnitude required to warrant
        a simulated alert.

        snr_cutoff is the minimum difference image signal to noise
        required to warrant a simulated alert.

        Outputs
        -------
        chip_name_dict is a dict keyed on i_obs (which is the index of
        an ObservationMetaData's position in self._full_obs_data['obs_dex'],
        NOT its position in self._obs_list).  The values of chip_name_dict
        are tuples containing:
            - a list of the names of the detectors that objects from chunk
              landed on (including Nones for those objects that did not land
              on any detector)

             - a list of the xPupil coords for every object in chunk

             - a list of the yPupil coords for every object in chunk

             - a list of the indexes in chunk of those objects which actually
               landed on a detector

        dmag_arr is a numpy array of the delta_magnitudes of every object
        in chunk.  dmag_arr[11][3][4] is the delta_magnitude of chunk[4]
        in the 3rd band (i.e. the i band) at TAI = expmjd[11].

        dmag_arr_transpose is dmag_arr with the time and object columns
        transposed so that dmag_arr_transpose[4][3][11] == dmag_arr[11][3][4].

        flux_arr is a numpy array with the total flux of the sources (quiescent +
        delta magnitude).  The first index corresponds to the ObservationMetaData
        of the observation.  The second index corresponds to the object.

        dflux_arr is a numpy array with the quiescent_flux-tot_flux of the sources.
        The first index corresponds to the ObservationMetaData
        of the observation.  The second index corresponds to the object.

        snr_arr is a numpy array with the total SNR of the sources (quiescent +
        delta magnitude).  The first index corresponds to the ObservationMetaData
        of the observation.  The second index corresponds to the object.

        time_arr is an array of integers with shape == (len(chunk), len(self._full_obs_data['obs_dex'])).
        A -1 in time_arr means that that combination of object and observation did
        not yield a valid observation.  A +1 means that the object and observation
        combination are valid.
        """

        dummy_sed = Sed()

        ######################################################
        # Calculate the delta_magnitude for all of the sources
        #
        photometry_catalog._set_current_chunk(chunk)
        dmag_arr = photometry_catalog.applyVariability(chunk['varParamStr'],
                                                       variability_cache=self._variability_cache,
                                                       expmjd=self._full_obs_data['mjd']).transpose((2, 0, 1))

        n_time = len(self._full_obs_data['mjd'])
        n_obj = len(chunk)
        tot_snr_array = np.zeros((n_time, n_obj), dtype=float)
        diff_snr_array = np.zeros((n_time, n_obj), dtype=float)
        flux_array = np.zeros((n_time, n_obj), dtype=float)
        dflux_array = np.zeros((n_time, n_obj), dtype=float)

        for i_obs, obs_dex in enumerate(self._full_obs_data['obs_dex']):
            obs = self._obs_list[obs_dex]
            i_filter = self.mag_name_to_int[obs.bandpass]
            gamma = self._full_obs_data['gamma'][i_obs]

            dmag = dmag_arr[i_obs][i_filter]
            q_mag = q_m_dict[i_filter]

            snr, new_gamma = calcSNR_m5(q_mag+dmag,
                                        self.bp_dict[obs.bandpass],
                                        obs.m5[obs.bandpass],
                                        self.phot_params,
                                        gamma=gamma)

            self._full_obs_data['gamma'][i_obs] = new_gamma
            tot_snr_array[i_obs] = snr

            tot_flux = dummy_sed.fluxFromMag(q_mag+dmag)
            tot_noise = tot_flux/snr
            q_flux = dummy_sed.fluxFromMag(q_mag)
            q_noise = q_flux/q_snr_dict[i_filter]

            diff_noise = np.sqrt(tot_noise**2 + q_noise**2)
            dflux = tot_flux-q_flux
            diff_snr = np.abs(dflux)/diff_noise
            diff_snr_array[i_obs] = diff_snr
            flux_array[i_obs] = tot_flux
            dflux_array[i_obs] = dflux

        dmag_arr_transpose = dmag_arr.transpose(2, 1, 0)
        diff_snr_transpose = diff_snr_array.transpose()

        snr_invalid = 0
        n_raw_obj = len(chunk)
        photometrically_valid = -1*np.ones(n_raw_obj, dtype=int)
        for i_obj in range(n_raw_obj):
            keep_it = False
            if snr_cutoff is None or diff_snr_transpose[i_obj].max() >= snr_cutoff:
                for i_filter in range(6):
                    if np.abs(dmag_arr_transpose[i_obj][i_filter]).max() >= dmag_cutoff:
                        keep_it = True
                        break

            if keep_it:
                photometrically_valid[i_obj] = 1

        photometrically_valid = np.where(photometrically_valid >= 0)

        if 'properMotionRa'in column_query:
            pmra = chunk['properMotionRa'][photometrically_valid]
            pmdec = chunk['properMotionDec'][photometrically_valid]
            px = chunk['parallax'][photometrically_valid]
            vrad = chunk['radialVelocity'][photometrically_valid]
        else:
            pmra = None
            pmdec = None
            px = None
            vrad = None

        ###################################################################
        # Figure out which sources actually land on an LSST detector during
        # the observations in question
        #
        chip_name_dict = {}

        # time_arr will keep track of which objects appear in which observations;
        # 1 means the object appears; -1 means it does not
        time_arr_transpose = -1*np.ones((len(self._full_obs_data['obs_dex']), len(chunk['raJ2000'])),
                                        dtype=int)

        for i_obs, obs_dex in enumerate(self._full_obs_data['obs_dex']):
            obs = self._obs_list[obs_dex]
            chip_name_list = np.array([None]*n_raw_obj)
            xpup_list = np.zeros(n_raw_obj, dtype=float)
            ypup_list = np.zeros(n_raw_obj, dtype=float)
            chip_int_arr = -1*np.ones(len(chip_name_list), dtype=int)

            if len(photometrically_valid[0]) > 0:
                xpup_list_val, ypup_list_val = _pupilCoordsFromRaDec(chunk['raJ2000'][photometrically_valid],
                                                                     chunk['decJ2000'][photometrically_valid],
                                                                     pm_ra=pmra, pm_dec=pmdec,
                                                                     parallax=px, v_rad=vrad,
                                                                     obs_metadata=obs)

                xpup_list[photometrically_valid] = xpup_list_val
                ypup_list[photometrically_valid] = ypup_list_val

                chip_name_list[photometrically_valid] = chipNameFromPupilCoordsLSST(xpup_list_val,
                                                                                    ypup_list_val)

                for i_chip, name in enumerate(chip_name_list):
                    if name is not None:
                        chip_int_arr[i_chip] = 1

            valid_obj = np.where(chip_int_arr > 0)
            time_arr_transpose[i_obs][valid_obj] = 1

            chip_name_dict[i_obs] = (chip_name_list,
                                     xpup_list,
                                     ypup_list,
                                     valid_obj)

        time_arr = time_arr_transpose.transpose()
        assert len(chip_name_dict) == len(self._full_obs_data['obs_dex'])

        return (chip_name_dict, dmag_arr, dmag_arr_transpose,
                flux_array, dflux_array, tot_snr_array, time_arr)

    def alert_data_from_htmid(self, htmid, dbobj,
                              dmag_cutoff=0.005,
                              snr_cutoff=None,
                              chunk_size=1000, write_every=10000,
                              output_dir='.', output_prefix='',
                              log_file_name=None,
                              photometry_class=None,
                              chunk_cutoff=-1,
                              lock=None):

        """
        Generate an sqlite file with all of the alert data for a given
        trixel.

        Parameters
        ----------
        htmid is an integer denoting the trixel from self.htmid_list that should
        be simulated

        dbobj is a CatalogDBObject connecting to the data underlying the simulation

        dmag_cutoff indicates the minimum change magnitude needed to trigger a
        simulated alert

        snr_cutoff indicates the minimum difference image SNR needed to trigger
        a simulated alert

        chunk_size denotes the number of objects to query from the database and
        process at one time

        write_every indicates how often to write to the sqlite file (i.e. the
        code will pause the simulation process and write to the sqlite file
        when it has accumulated this many valid observations)

        output_dir is the directory in which to create the sqlite file

        output_prefix is the prefix of the sqlite file's name

        log_file_name is the name of a text file where progress will be written

        photometry_class is a InstanceCatalog class (not an instantiation) that
        contains the methods for calculating the photometry associated with the
        simulated alerts (see AlertStellarVariabilityCatalog and
        AlertAgnVariabilityCatalog in this module)

        chunk_cutoff is an optional int; stop the simulation after this many
        chunks have been processed.  This is for testing purposes.
        If chunk_cutoff == -1, the code will process all of the astrophysical
        objects in the trixel.

        lock is a multiprocessing.Lock() for use if running multiple
        instances of alert_data_from_htmid.  This will prevent multiple processes
        from writing to the log file or stdout simultaneously.
        """

        htmid_level = levelFromHtmid(htmid)
        if log_file_name is None:
            raise RuntimeError('must specify log_file_name')

        if ('_PARAMETRIZED_LC_DMAG_LOOKUP' not in self._variability_cache and
            self._dmag_lookup_file_exists):

            self._variability_cache['_PARAMETRIZED_LC_DMAG_CUTOFF'] = dmag_cutoff
            self._variability_cache['_PARAMETRIZED_LC_DMAG_LOOKUP'] = {}

            with open(self._dmag_lookup_file, 'r') as in_file:
                for line in in_file:
                    if line[0] == '#':
                        continue
                    params = line.split()
                    self._variability_cache['_PARAMETRIZED_LC_DMAG_LOOKUP'][int(params[0])] = float(params[1])

        self._stdout_lock = lock
        this_pid = os.getpid()

        t_start = time.time()  # so that we can get a sense of how long the full
                               # simulation will take

        desired_columns = []
        desired_columns.append('simobjid')
        desired_columns.append('variabilityParameters')
        desired_columns.append('varParamStr')
        desired_columns.append('raJ2000')
        desired_columns.append('decJ2000')
        desired_columns.append('properMotionRa')
        desired_columns.append('properMotionDec')
        desired_columns.append('parallax')
        desired_columns.append('radialVelocity')
        desired_columns.append('ebv')
        desired_columns.append('redshift')
        desired_columns.append('htmid')

        if 'umag' in dbobj.columnMap:
            desired_columns.append('umag')
            desired_columns.append('gmag')
            desired_columns.append('rmag')
            desired_columns.append('imag')
            desired_columns.append('zmag')
            desired_columns.append('ymag')
        elif 'u_ab' in dbobj.columnMap:
            desired_columns.append('u_ab')
            desired_columns.append('g_ab')
            desired_columns.append('r_ab')
            desired_columns.append('i_ab')
            desired_columns.append('z_ab')
            desired_columns.append('y_ab')
        else:
            raise RuntimeError('Not sure what quiescent '
                               'LSST magnitudes are called '
                               'in this CatalogDBObject')

        if photometry_class is None:
            raise RuntimeError('Must specify photometry_class')

        if os.path.exists(output_dir) and not os.path.isdir(output_dir):
            raise RuntimeError('%s is not a dir' % output_dir)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        dummy_sed = Sed()

        # a dummy call to make sure that the initialization
        # is done before we attempt to parallelize calls
        # to chipNameFromRaDecLSST
        chipNameFromPupilCoordsLSST(0.0, 0.0)

        mag_names = ('u', 'g', 'r', 'i', 'z', 'y')

        # from Table 2 of the overview paper
        obs_mag_cutoff = (23.68, 24.89, 24.43, 24.0, 24.45, 22.60)

        gamma_template = {}
        for i_filter in range(6):
            gamma_template[i_filter] = None

        local_obs_valid_dex = self._htmid_dict[htmid]
        print('n valid obs %d' % len(local_obs_valid_dex))

        # store sorted lists of MJD, indices of self._obs_list
        # and photometric catalog classes in self._full_obs_data

        self._full_obs_data = {}

        local_cat_list = []
        local_expmjd_list = []

        for obs_dex in local_obs_valid_dex:
            obs = self._obs_list[obs_dex]
            cat = photometry_class(dbobj, obs_metadata=obs)
            cat.lsstBandpassDict = self.bp_dict
            cat.phot_params = self.phot_params
            local_cat_list.append(cat)
            local_expmjd_list.append(obs.mjd.TAI)

        local_expmjd_list = np.array(local_expmjd_list)
        local_cat_list = np.array(local_cat_list)
        sorted_dex = np.argsort(local_expmjd_list)

        local_expmjd_list = local_expmjd_list[sorted_dex]
        local_cat_list = local_cat_list[sorted_dex]
        local_obs_valid_dex = local_obs_valid_dex[sorted_dex]

        self._full_obs_data['mjd'] = local_expmjd_list
        self._full_obs_data['obs_dex'] = local_obs_valid_dex
        self._full_obs_data['cat'] = local_cat_list
        self._full_obs_data['gamma'] = [None]*len(local_cat_list)

        available_columns = list(dbobj.columnMap.keys())
        column_query = []
        for col in desired_columns:
            if col in available_columns:
                column_query.append(col)

        n_bits_off = 2*(21-htmid_level)

        data_iter = dbobj.query_columns_htmid(colnames=column_query,
                                              htmid=htmid,
                                              chunk_size=chunk_size)

        photometry_catalog = photometry_class(dbobj, self._obs_list[self._full_obs_data['obs_dex'][0]],
                                              column_outputs=['lsst_u',
                                                              'lsst_g',
                                                              'lsst_r',
                                                              'lsst_i',
                                                              'lsst_z',
                                                              'lsst_y'])
        photometry_catalog.lsstBandpassDict = self.bp_dict

        i_chunk = 0

        output_data_cache = {}
        n_rows_cached = 0

        n_obj = 0
        n_actual_obj = 0
        n_time_last = 0
        n_rows = 0

        t_before_obj = time.time()  # so that we can get a sense of how long the
                                    # "iterating over astrophysical objects" part
                                    # of the simulation will take

        db_name = os.path.join(output_dir, '%s_%d_sqlite.db' % (output_prefix, htmid))
        with sqlite3.connect(db_name, isolation_level='EXCLUSIVE') as conn:
            creation_cmd = '''CREATE TABLE alert_data
                           (uniqueId int, obshistId int, xPix float, yPix float,
                            chipNum int, dflux float, snr float, ra float, dec float)'''

            cursor = conn.cursor()
            cursor.execute('PRAGMA journal_mode=WAL;')
            conn.commit()
            cursor.execute(creation_cmd)
            conn.commit()

            creation_cmd = '''CREATE TABLE metadata
                           (obshistId int, TAI float, band int)'''
            cursor.execute(creation_cmd)
            conn.commit()

            for obs_dex in self._full_obs_data['obs_dex']:
                obs = self._obs_list[obs_dex]
                cmd = '''INSERT INTO metadata
                      VALUES(%d, %.5f, %d)''' % (obs.OpsimMetaData['obsHistID'],
                                                 obs.mjd.TAI,
                                                 self.mag_name_to_int[obs.bandpass])

                cursor.execute(cmd)
            conn.commit()

            creation_cmd = '''CREATE TABLE quiescent_flux
                          (uniqueId int, band int, flux float, snr float)'''

            cursor.execute(creation_cmd)
            conn.commit()

            creation_cmd = '''CREATE TABLE baseline_astrometry
                           (uniqueId int, ra real, dec real, pmRA real,
                            pmDec real, parallax real, TAI real)'''

            cursor.execute(creation_cmd)
            conn.commit()

            for chunk in data_iter:
                n_raw_obj = len(chunk)
                i_chunk += 1

                if chunk_cutoff > 0 and i_chunk >= chunk_cutoff:
                    break

                n_time_last = 0

                # filter the chunk so that we are only considering sources that are in
                # the trixel being considered
                reduced_htmid = chunk['htmid'] >> n_bits_off

                valid_htmid = np.where(reduced_htmid == htmid)
                if len(valid_htmid[0]) == 0:
                    continue
                n_htmid_trim = n_raw_obj-len(valid_htmid[0])
                chunk = chunk[valid_htmid]
                n_obj += len(valid_htmid[0])

                q_f_dict = {}
                q_m_dict = {}

                photometry_catalog._set_current_chunk(chunk)

                q_m_dict[0] = photometry_catalog.column_by_name('quiescent_lsst_u')
                q_m_dict[1] = photometry_catalog.column_by_name('quiescent_lsst_g')
                q_m_dict[2] = photometry_catalog.column_by_name('quiescent_lsst_r')
                q_m_dict[3] = photometry_catalog.column_by_name('quiescent_lsst_i')
                q_m_dict[4] = photometry_catalog.column_by_name('quiescent_lsst_z')
                q_m_dict[5] = photometry_catalog.column_by_name('quiescent_lsst_y')

                q_f_dict[0] = dummy_sed.fluxFromMag(q_m_dict[0])
                q_f_dict[1] = dummy_sed.fluxFromMag(q_m_dict[1])
                q_f_dict[2] = dummy_sed.fluxFromMag(q_m_dict[2])
                q_f_dict[3] = dummy_sed.fluxFromMag(q_m_dict[3])
                q_f_dict[4] = dummy_sed.fluxFromMag(q_m_dict[4])
                q_f_dict[5] = dummy_sed.fluxFromMag(q_m_dict[5])

                q_pmra = 1000.0*arcsecFromRadians(photometry_catalog.column_by_name('properMotionRa'))
                q_pmdec = 1000.0*arcsecFromRadians(photometry_catalog.column_by_name('properMotionDec'))
                q_parallax = 1000.0*arcsecFromRadians(photometry_catalog.column_by_name('parallax'))
                q_ra = np.degrees(photometry_catalog.column_by_name('raICRS'))
                q_dec = np.degrees(photometry_catalog.column_by_name('decICRS'))
                q_tai = photometry_catalog.obs_metadata.mjd.TAI

                q_snr_dict = {}
                for i_filter in range(6):

                    snr_template, local_gamma = calcSNR_m5(q_m_dict[i_filter],
                                                           self.bp_dict[mag_names[i_filter]],
                                                           obs_mag_cutoff[i_filter],
                                                           self.phot_params, gamma=gamma_template[i_filter])
                    q_snr_dict[i_filter] = snr_template
                    gamma_template[i_filter] = local_gamma

                unq = photometry_catalog.column_by_name('uniqueId')

                # find chipName, but only for objects with variability amplitude that
                # exceeds dmag_cutoff

                (chip_name_dict,
                 dmag_arr,
                 dmag_arr_transpose,
                 flux_arr,
                 dflux_arr,
                 snr_arr,
                 time_arr) = self._filter_on_photometry_then_chip_name(chunk, column_query,
                                                                       photometry_catalog,
                                                                       q_m_dict,
                                                                       q_snr_dict,
                                                                       dmag_cutoff, snr_cutoff)

                try:
                    assert dmag_arr_transpose.shape == (len(chunk), len(mag_names), len(self._full_obs_data['mjd']))
                except AssertionError:
                    print('dmag_arr_transpose_shape %s' % str(dmag_arr_transpose.shape))
                    print('should be (%d, %d, %d)' % (len(chunk), len(mag_names), len(self._full_obs_data['mjd'])))
                    raise

                # only include those sources for which np.abs(delta_mag) >= dmag_cutoff
                # at some point in their history (note that delta_mag is defined with
                # respect to the quiescent magnitude)
                #
                # also demand that the magnitude at some point is less than obs_mag_cutoff
                #
                # This is different from the dmag>dmag_cutoff check done in
                #
                # self._filter_on_photometry_then_chip_name()
                #
                # Now we have information about which observations actually detected
                # each object (in self._filter_on_photometry_then_chip_name(),
                # we assumed that every object was detected at every time step).

                photometrically_valid_obj = []
                for i_obj in range(len(chunk)):
                    keep_it = False
                    valid_times = np.where(time_arr[i_obj] > 0)
                    if len(valid_times[0]) == 0:
                        continue
                    for i_filter in range(len(mag_names)):
                        if np.abs(dmag_arr_transpose[i_obj][i_filter][valid_times]).max() > dmag_cutoff:
                            dmag_min = dmag_arr_transpose[i_obj][i_filter][valid_times].min()
                            if q_m_dict[i_filter][i_obj] + dmag_min <= obs_mag_cutoff[i_filter]:
                                keep_it = True
                                break
                    if keep_it:
                        photometrically_valid_obj.append(i_obj)
                photometrically_valid_obj = np.array(photometrically_valid_obj)

                del dmag_arr_transpose
                gc.collect()

                if np.abs(dmag_arr).max() < dmag_cutoff:
                    continue

                completely_valid = np.zeros(len(chunk), dtype=int)

                ############################
                # Process and output sources
                #
                for i_obs, obs_dex in enumerate(self._full_obs_data['obs_dex']):
                    obs = self._obs_list[obs_dex]
                    obshistid = obs.OpsimMetaData['obsHistID']

                    obs_mag = obs.bandpass
                    actual_i_mag = self.mag_name_to_int[obs_mag]
                    assert mag_names[actual_i_mag] == obs_mag

                    # only include those sources which fall on a detector for this pointing
                    valid_chip_name, valid_xpup, valid_ypup, chip_valid_obj = chip_name_dict[i_obs]

                    actually_valid_obj = np.intersect1d(photometrically_valid_obj, chip_valid_obj)
                    if len(actually_valid_obj) == 0:
                        continue

                    try:
                        completely_valid[actually_valid_obj] += 1
                    except:
                        print('failed')
                        print(actually_valid_obj)
                        print(completely_valid)
                        raise

                    valid_sources = chunk[actually_valid_obj]
                    local_column_cache = {}
                    local_column_cache['deltaMagAvro'] = OrderedDict([('delta_%smag' % mag_names[i_mag],
                                                                      dmag_arr[i_obs][i_mag][actually_valid_obj])
                                                                      for i_mag in range(len(mag_names))])

                    local_column_cache['flux'] = flux_arr[i_obs][actually_valid_obj]
                    local_column_cache['dflux'] = dflux_arr[i_obs][actually_valid_obj]
                    local_column_cache['SNR'] = snr_arr[i_obs][actually_valid_obj]
                    local_column_cache['chipName'] = valid_chip_name[actually_valid_obj]
                    local_column_cache['pupilFromSky'] = OrderedDict([('x_pupil', valid_xpup[actually_valid_obj]),
                                                                      ('y_pupil', valid_ypup[actually_valid_obj])])

                    cat = self._full_obs_data['cat'][i_obs]
                    i_valid_chunk = 0
                    for valid_chunk, chunk_map in cat.iter_catalog_chunks(query_cache=[valid_sources],
                                                                          column_cache=local_column_cache):
                        i_valid_chunk += 1
                        assert i_valid_chunk == 1
                        n_time_last += len(valid_chunk[0])
                        length_of_chunk = len(valid_chunk[chunk_map['uniqueId']])
                        cache_tag = '%d_%d' % (obshistid, i_chunk)
                        output_data_cache[cache_tag] = {}

                        for col_name in self._output_columns:
                            output_data_cache[cache_tag][col_name] = valid_chunk[chunk_map[col_name]]

                        n_rows_cached += length_of_chunk

                completely_valid = np.where(completely_valid > 0)
                for i_filter in range(6):
                    values = ((int(unq[completely_valid][i_q]),
                               i_filter,
                               q_f_dict[i_filter][completely_valid][i_q],
                               q_snr_dict[i_filter][completely_valid][i_q])
                              for i_q in range(len(completely_valid[0])))
                    cursor.executemany('INSERT INTO quiescent_flux VALUES (?,?,?,?)', values)
                    conn.commit()

                values = ((int(unq[completely_valid][i_q]),
                           q_ra[completely_valid][i_q],
                           q_dec[completely_valid][i_q],
                           q_pmra[completely_valid][i_q],
                           q_pmdec[completely_valid][i_q],
                           q_parallax[completely_valid][i_q],
                           q_tai)
                          for i_q in range(len(completely_valid[0])))

                cursor.executemany('INSERT INTO baseline_astrometry VALUES (?,?,?,?,?,?,?)', values)

                if n_rows_cached >= write_every:
                    with _lock_context(self._stdout_lock):
                        with open(log_file_name, 'a') as out_file:
                            out_file.write('%d is writing \n' % os.getpid())

                            print('%d is writing' % os.getpid())

                    n_rows += self._output_alert_data(conn, output_data_cache)
                    output_data_cache = {}
                    n_rows_cached = 0

                    if n_rows > 0:
                        with _lock_context(self._stdout_lock):
                            with open(log_file_name, 'a') as out_file:
                                elapsed = (time.time()-t_before_obj)/3600.0
                                elapsed_per = elapsed/n_rows
                                rows_per_chunk = float(n_rows)/float(i_chunk)
                                total_projection = 1000.0*rows_per_chunk*elapsed_per
                                out_file.write('\n    %d n_obj %d %d trimmed %d\n' %
                                               (this_pid, n_obj, n_actual_obj, n_htmid_trim))
                                out_file.write('    elapsed %.2e hrs per row %.2e total %2e\n' %
                                               (elapsed, elapsed_per, total_projection))
                                out_file.write('    n_time_last %d; rows %d\n' % (n_time_last, n_rows))

                                out_file.write('%d is done writing\n' % os.getpid())

                                print('\n    %d n_obj %d %d trimmed %d' %
                                      (this_pid, n_obj, n_actual_obj, n_htmid_trim))
                                print('    elapsed %.2e hrs per row %.2e total %2e' %
                                      (elapsed, elapsed_per, total_projection))
                                print('    n_time_last %d; rows %d\n' % (n_time_last, n_rows))
                                print('%d is done writing' % os.getpid())

            if len(output_data_cache) > 0:
                n_rows += self._output_alert_data(conn, output_data_cache)
                output_data_cache = {}

            print('htmid %d that took %.2e hours; n_obj %d n_rows %d' %
                  (htmid, (time.time()-t_start)/3600.0, n_obj, n_rows))

            with _lock_context(self._stdout_lock):
                print("INDEXING %d" % htmid)

            cursor.execute('CREATE INDEX unq_obs ON alert_data (uniqueId, obshistId)')
            cursor.execute('CREATE INDEX unq_flux ON quiescent_flux (uniqueId, band)')
            cursor.execute('CREATE INDEX obs ON metadata (obshistid)')
            cursor.execute('CREATE INDEX unq_ast ON baseline_astrometry (uniqueId)')
            conn.commit()

            with _lock_context(self._stdout_lock):
                with open(log_file_name, 'a') as out_file:
                    out_file.write('done with htmid %d -- %e %d\n' %
                                   (htmid, (time.time()-t_start)/3600.0, n_obj))

        return n_rows
