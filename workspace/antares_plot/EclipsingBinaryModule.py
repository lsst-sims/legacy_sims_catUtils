from lsst.sims.catUtils.baseCatalogModels import StarBase
from lsst.sims.catUtils.utils import StellarAlertDBObjMixin

class OBAFGKObj(StellarAlertDBObjMixin, StarBase):
    objid = 'obafgkstars'
    tableid = 'StarOBAFGKForceseek'
    objectTypeId = 14
    doRunTest = True
    #These types should be matched to the database.
    #: Default map is float.  If the column mapping is the same as the column name, None can be specified
    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('glon', 'gal_l*PI()/180.'),
               ('glat', 'gal_b*PI()/180.'),
               ('magNorm', '(-2.5*log(flux_scale)/log(10.)) - 18.402732642'),
               ('properMotionRa', '(mura/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('parallax', 'parallax*PI()/648000000.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', str, 40)]

from lsst.sims.catUtils.utils import AlertDataGenerator
from lsst.sims.coordUtils import chipNameFromPupilCoordsLSST
from lsst.sims.utils import _pupilCoordsFromRaDec
from lsst.utils import getPackageDir
import json
import time
import os
import numpy as np

class EclipsingBinaryAlertDataGenerator(AlertDataGenerator):

    def _filter_on_photometry_then_chip_name(self, chunk, column_query,
                                             obs_valid_dex, expmjd_list,
                                             photometry_catalog,
                                             dmag_cutoff):
        """
        Determine which simulated observations are actually worth storing
        by first figuring out which observations of which objects are
        photoemtrically detectable and alert-worthy, then determining
        which of those actually fall on an LSST detector.

        Parameters
        ----------
        chunk is the output yielded from a CatSim ChunkIterator.  It is
        a numpy recarray representing one chunk_size query from the
        underlying simulations database

        column_query is a list of the columns that were queried from
        the database

        obs_valid_dex is a list of integers corresponding to indexes in
        self._obs_list of the ObservationMetaData that are actually valid
        for the trixel currently being simulated

        expmjd_list is a numpy array of the TAI dates of ObservtionMetaData
        represented by obs_valid_dex

        photometry_catalog is an instantiation of the InstanceCatalog class
        being used to calculate magnitudes for these variable sources.

        Outputs
        -------
        chip_name_dict is a dict keyed on i_obs (which is the index of
        an ObservationMetaData's position in obs_valid_dex, NOT its
        position in self._obs_list).  The values of chip_name_dict are
        tuples containing:
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

        time_arr is an array of integers with shape == (len(chunk), len(obs_valid_dex)).
        A -1 in time_arr means that that combination of object and observation did
        not yield a valid observation.  A +1 means that the object and observation
        combination are valid.
        """
        if not hasattr(self, '_eb_id_set'):
            self._eb_id_set = set()

            eb_cat = os.path.join(getPackageDir('sims_catUtils'),
                                  'workspace', 'antares_plot', 'data',
                                  'villanova_eb_catalog.txt')

            with open(eb_cat, 'r') as in_file:
                for line in in_file:
                    if line[0] == '#':
                        continue
                    params = line.strip().split()
                    self._eb_id_set.add(int(params[0]))

        for i_star in range(len(chunk)):
            param_dict = json.loads(chunk['varParamStr'][i_star])
            kep_id = param_dict['p']['lc']
            if kep_id not in self._eb_id_set:
                chunk['varParamStr'][i_star] = 'None'

        photometry_catalog._set_current_chunk(chunk)
        dmag_arr = photometry_catalog.applyVariability(chunk['varParamStr'],
                                                       variability_cache=self._variability_cache,
                                                       expmjd=expmjd_list,).transpose((2,0,1))

        dmag_arr_transpose = dmag_arr.transpose(2,1,0)

        n_raw_obj = len(chunk)
        photometrically_valid = -1*np.ones(n_raw_obj, dtype=int)
        for i_obj in range(n_raw_obj):
            keep_it = False
            for i_filter in range(6):
               if np.abs(dmag_arr_transpose[i_obj][i_filter]).max() >= dmag_cutoff:
                   keep_it = True
                   break
            if keep_it:
                photometrically_valid[i_obj] = 1

        photometrically_valid = np.where(photometrically_valid>=0)

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
        t_before_chip_name = time.time()
        chip_name_dict = {}

        # time_arr will keep track of which objects appear in which observations;
        # 1 means the object appears; -1 means it does not
        time_arr_transpose = -1*np.ones((len(obs_valid_dex), len(chunk['raJ2000'])),
                                        dtype=int)

        for i_obs, obs_dex in enumerate(obs_valid_dex):
            obs = self._obs_list[obs_dex]
            chip_name_list = np.array([None]*n_raw_obj)
            xPup_list = np.zeros(n_raw_obj, dtype=float)
            yPup_list = np.zeros(n_raw_obj, dtype=float)
            chip_int_arr = -1*np.ones(len(chip_name_list), dtype=int)

            if len(photometrically_valid[0])>0:
                xPup_list_val, yPup_list_val = _pupilCoordsFromRaDec(chunk['raJ2000'][photometrically_valid],
                                                                     chunk['decJ2000'][photometrically_valid],
                                                                     pm_ra=pmra, pm_dec=pmdec,
                                                                     parallax=px, v_rad=vrad,
                                                                     obs_metadata=obs)

                xPup_list[photometrically_valid] = xPup_list_val
                yPup_list[photometrically_valid] = yPup_list_val

                chip_name_list[photometrically_valid] = chipNameFromPupilCoordsLSST(xPup_list_val,
                                                                                    yPup_list_val)

                for i_chip, name in enumerate(chip_name_list):
                    if name is not None:
                        chip_int_arr[i_chip] = 1

            valid_obj = np.where(chip_int_arr>0)
            time_arr_transpose[i_obs][valid_obj] = 1

            chip_name_dict[i_obs] = (chip_name_list,
                                     xPup_list,
                                     yPup_list,
                                     valid_obj)

        time_arr = time_arr_transpose.transpose()
        assert len(chip_name_dict) == len(obs_valid_dex)

        ######################################################
        # Calculate the delta_magnitude for all of the sources
        #
        t_before_phot = time.time()

        # only calculate photometry for objects that actually land
        # on LSST detectors

        return chip_name_dict, dmag_arr, dmag_arr_transpose, time_arr
