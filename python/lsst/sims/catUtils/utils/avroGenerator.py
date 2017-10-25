import numpy as np
import h5py
import multiprocessing as mproc
from collections import OrderedDict
import time
from lsst.sims.utils import findHtmid, trixelFromHtmid
from lsst.sims.utils import angularSeparation, ObservationMetaData
from lsst.sims.catUtils.utils import _baseLightCurveCatalog
from lsst.sims.coordUtils import _chipNameFromRaDecLSST
from lsst.sims.coordUtils import chipNameFromPupilCoordsLSST

from lsst.sims.catalogs.decorators import compound
from lsst.sims.photUtils import BandpassDict
from lsst.sims.catUtils.mixins import VariabilityStars, AstrometryStars
from lsst.sims.catUtils.mixins import CameraCoordsLSST, PhotometryBase
from lsst.sims.catUtils.mixins import ParametrizedLightCurveMixin
from lsst.sims.catUtils.mixins import create_variability_cache

__all__ = ["AvroGenerator"]


class _baseAvroCatalog(_baseLightCurveCatalog):

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
            # Call the originalversion of iter_catalog defined in the
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
                    self._chunkColMap_output = dict([(col, i) for i, col in enumerate(self.iter_column_names())])
                yield chunk_cols, self._chunkColMap_output


class StellarVariabilityCatalog(VariabilityStars, AstrometryStars, PhotometryBase,
                                CameraCoordsLSST, _baseAvroCatalog):
    column_outputs = ['uniqueId', 'raICRS', 'decICRS',
                      'mag','mag_uncertainty', 'dmag', 'chipName', 'varParamStr']

    default_formats = {'f':'%.4g'}

    @compound('sigma_lsst_u', 'sigma_lsst_g', 'sigma_lsst_r',
              'sigma_lsst_i', 'sigma_lsst_z', 'sigma_lsst_y')
    def get_lsst_photometric_uncertainties(self):
        if not hasattr(self, 'lsstBandpassDict'):
            self.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()
        return self._magnitudeUncertaintyGetter(['lsst_u', 'lsst_g', 'lsst_r',
                                                 'lsst_i', 'lsst_z', 'lsst_y'],
                                                ['u', 'g', 'r', 'i', 'z', 'y'],
                                                'lsstBandpassDict')

    @compound('delta_umag', 'delta_gmag', 'delta_rmag',
              'delta_imag', 'delta_zmag', 'delta_ymag')
    def get_deltaMagAvro(self):
        ra = self.column_by_name('raJ2000')
        if len(ra)==0:
            return np.array([[],[],[],[],[],[]])

        raise RuntimeError("Should not have gotten this far in delta mag getter")


    @compound('quiescent_lsst_u', 'quiescent_lsst_g', 'quiescent_lsst_r',
              'quiescent_lsst_i', 'quiescent_lsst_z', 'quiescent_lsst_y')
    def get_quiescent_lsst_magnitudes(self):
        return np.array([self.column_by_name('umag'), self.column_by_name('gmag'),
                         self.column_by_name('rmag'), self.column_by_name('imag'),
                         self.column_by_name('zmag'), self.column_by_name('ymag')])

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

        delta = np.array([self.column_by_name('delta_umag'),
                          self.column_by_name('delta_gmag'),
                          self.column_by_name('delta_rmag'),
                          self.column_by_name('delta_imag'),
                          self.column_by_name('delta_zmag'),
                          self.column_by_name('delta_ymag')])
        magnitudes += delta

        return magnitudes

    @compound('mag', 'mag_uncertainty', 'dmag')
    def get_avroPhotometry(self):
        mag = self.column_by_name('lsst_%s' % self.obs_metadata.bandpass)
        mag_unc = self.column_by_name('sigma_lsst_%s' % self.obs_metadata.bandpass)
        dmag = self.column_by_name('%smag' % self.obs_metadata.bandpass) - mag

        return np.array([mag, mag_unc, dmag])


def _find_chipNames_parallel(ra, dec, pm_ra=None, pm_dec=None, parallax=None,
                             v_rad=None, obs_metadata_list=None, i_obs_list=None, out_dict=None):

    for i_obs, obs in zip(i_obs_list, obs_metadata_list):
        chip_name_list = _chipNameFromRaDecLSST(ra, dec, pm_ra=pm_ra, pm_dec=pm_dec,
                                                parallax=parallax, v_rad=v_rad,
                                                obs_metadata=obs)

        chip_int_arr = -1*np.ones(len(chip_name_list), dtype=int)
        for i_chip, name in enumerate(chip_name_list):
            if name is not None:
                chip_int_arr[i_chip] = 1
        valid_obj = np.where(chip_int_arr>0)

        out_dict[i_obs] = (chip_name_list[valid_obj], valid_obj)


class AvroGenerator(object):

    def __init__(self, obs_list, n_proc_max=4):
        self._t_chip_name=0.0
        self._t_mlt = 0.0
        self._t_param_lc = 0.0
        self._t_apply_var = 0.0
        self._t_setup = 0.0
        self._t_phot = 0.0
        self._t_out = 0.0
        self._t_filter_phot = 0.0

        self._output_prefix = 'test_hdf5'
        self._dmag_cutoff = 0.005
        self._n_proc_max = n_proc_max
        self._variability_cache = create_variability_cache()
        plm = ParametrizedLightCurveMixin()
        plm.load_parametrized_light_curves(variability_cache = self._variability_cache)
        self.obs_list = np.array(obs_list)
        htmid_level = 7
        self._htmid_list = []
        for obs in self.obs_list:
            htmid = findHtmid(obs.pointingRA, obs.pointingDec, htmid_level)
            self._htmid_list.append(htmid)
        self._htmid_list = np.array(self._htmid_list)
        self._unq_htmid_list = np.unique(self._htmid_list)
        self.bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
        self.chunk_size = 10000
        self._desired_columns = []
        self._desired_columns.append('simobjid')
        self._desired_columns.append('variabilityParameters')
        self._desired_columns.append('varParamStr')
        self._desired_columns.append('raJ2000')
        self._desired_columns.append('decJ2000')
        self._desired_columns.append('properMotionRa')
        self._desired_columns.append('properMotionDec')
        self._desired_columns.append('parallax')
        self._desired_columns.append('radialVelocity')
        self._desired_columns.append('umag')
        self._desired_columns.append('gmag')
        self._desired_columns.append('rmag')
        self._desired_columns.append('imag')
        self._desired_columns.append('zmag')
        self._desired_columns.append('ymag')
        self._desired_columns.append('ebv')
        self._desired_columns.append('redshift')
        print('initialized with %d %d' %
              (len(self.obs_list), len(self._unq_htmid_list)))

    @property
    def htmid_list(self):
        """
        A list of the unique htmid's corresponding to the fields
        of view that need to be queried to generate the alert data
        """
        return self._unq_htmid_list

    def alerts_from_db(self, dbobj):
        n_obs_total = 0
        t_0 = time.time()

        for i_h, htmid in enumerate(self._unq_htmid_list):
            print('processing %d --- %d of %d' % (htmid, i_h, len(self._unq_htmid_list)))
            t_start = time.time()
            n_obs = self.alert_data_from_htmid(htmid, dbobj)
            n_obs_total += n_obs
            print("that took %e hours" % ((time.time()-t_start)/3600.0))
            print("total should take %e hours" %
            (len(self.obs_list)*(time.time()-t_0)/(3600.0*n_obs_total)))
            print('total took %.2e chip name %.2e MLT %.2e paramLC %.2e' %
            (time.time()-t_0, self._t_chip_name,self._t_mlt,self._t_param_lc))
            print('applyVar %.2e setup %.2e' % (self._t_apply_var, self._t_setup))
            print('phot %.2e output %.2e' % (self._t_phot, self._t_out))
            print('filter_phot %.2e' % self._t_filter_phot)
            if i_h>2:
                exit()

    def output_to_hdf5(self, hdf5_file, data_cache):
        """
        Cache will be keyed first on the obsHistID, then all of the columns
        """
        self._output_ct += 1
        for obsHistID in data_cache.keys():
            for col_name in data_cache[obsHistID].keys():
                data_tag = '%d_%d_%s' % (obsHistID, self._output_ct, col_name)
                hdf5_file.create_dataset(data_tag, data=np.array(data_cache[obsHistID][col_name]))

        hdf5_file.flush()

    def alert_data_from_htmid(self, htmid, dbobj, radius=1.75):

        t_start = time.time()

        self._output_ct = -1
        # out_file = h5py.File('%s_%d.hdf5' % (self._output_prefix, htmid), 'w')

        # a dummy call to make sure that the initialization
        # is done before we attempt to parallelize calls
        # to chipNameFromRaDecLSST
        dummy_name = chipNameFromPupilCoordsLSST(0.0, 0.0)

        mag_names = ('u', 'g', 'r', 'i', 'z', 'y')
        valid_dexes = np.where(self._htmid_list == htmid)
        print('valid_dexes %s ' % str(valid_dexes))
        obs_valid = self.obs_list[valid_dexes]
        center_trixel = trixelFromHtmid(htmid)
        center_ra, center_dec = center_trixel.get_center()

        ra_list = []
        dec_list = []
        cat_list = []
        expmjd_list = []
        obshistid_list = []
        for obs in obs_valid:
            ra_list.append(obs.pointingRA)
            dec_list.append(obs.pointingDec)
            cat = StellarVariabilityCatalog(dbobj, obs_metadata=obs)
            cat.lsstBandpassDict =  self.bp_dict
            cat_list.append(cat)
            expmjd_list.append(obs.mjd.TAI)
            obshistid_list.append(obs.OpsimMetaData['obsHistID'])

        expmjd_list = np.array(expmjd_list)
        obshistid_list = np.array(obshistid_list)
        sorted_dex = np.argsort(expmjd_list)

        # out_file.create_dataset('obshistID', data=obshistid_list)
        # out_file.create_dataset('TAI', data=expmjd_list)
        # out_file.flush()

        expmjd_list = expmjd_list[sorted_dex]
        ra_list = np.array(ra_list)[sorted_dex]
        dec_list = np.array(dec_list)[sorted_dex]
        cat_list = np.array(cat_list)[sorted_dex]

        print('built list')

        dist_list = angularSeparation(center_ra, center_dec, ra_list, dec_list)
        radius += dist_list.max()
        center_obs = ObservationMetaData(pointingRA=center_ra,
                                         pointingDec=center_dec,
                                         boundType='circle',
                                         boundLength=radius)

        print('radius %e' % radius)

        available_columns = list(dbobj.columnMap.keys())
        column_query = []
        for col in self._desired_columns:
            if col in available_columns:
                column_query.append(col)

        data_iter = dbobj.query_columns(colnames=column_query,
                                        obs_metadata=center_obs,
                                        chunk_size=self.chunk_size)


        photometry_catalog = StellarVariabilityCatalog(dbobj, obs_valid[0],
                                                       column_outputs=['lsst_u',
                                                                       'lsst_g',
                                                                       'lsst_r',
                                                                       'lsst_i',
                                                                       'lsst_z',
                                                                       'lsst_y'])

        self._t_setup += time.time()-t_start

        print('chunking')
        i_chunk = 0
        t_chipName = 0.0

        n_proc_possible = int(np.ceil(len(obs_valid)/5.0))
        n_proc_chipName = min(n_proc_possible, self._n_proc_max)

        output_data_cache = {}
        ct_to_write = 0

        for chunk in data_iter:
            i_chunk += 1
            if 'properMotionRa'in column_query:
                pmra = chunk['properMotionRa']
                pmdec = chunk['properMotionDec']
                px = chunk['parallax']
                vrad = chunk['radialVelocity']
            else:
                pmra = None
                pmdec = None
                px = None
                vrad = None

            #for ii in range(6):
            #    print('dmag %d: %e %e %e' % (ii,dmag_arr[ii].min(),np.median(dmag_arr[ii]),dmag_arr[ii].max()))
            #exit()

            ###################################################################
            # Figure out which sources actually land on an LSST detector during
            # the observations in question
            #
            t_before_chip_name = time.time()
            if n_proc_chipName == 1:
                chip_name_dict = {}
            else:
                mgr = mproc.Manager()
                chip_name_dict = mgr.dict()
                iobs_sub_list = []
                obs_sub_list = []
                for i_obs in range(n_proc_chipName):
                    iobs_sub_list.append([])
                    obs_sub_list.append([])
                sub_list_ct = 0

            for i_obs, obs in enumerate(obs_valid):
                if n_proc_chipName == 1:
                    chip_name_list = _chipNameFromRaDecLSST(chunk['raJ2000'],
                                                            chunk['decJ2000'],
                                                            pm_ra=pmra,
                                                            pm_dec=pmdec,
                                                            parallax=px,
                                                            v_rad=vrad,
                                                            obs_metadata=obs)

                    chip_int_arr = -1*np.ones(len(chip_name_list), dtype=int)
                    for i_chip, name in enumerate(chip_name_list):
                        if name is not None:
                            chip_int_arr[i_chip] = 1

                    valid_obj = np.where(chip_int_arr>0)
                    chip_name_dict[i_obs] = (chip_name_list[valid_obj], valid_obj)

                else:
                    iobs_sub_list[sub_list_ct].append(i_obs)
                    obs_sub_list[sub_list_ct].append(obs)
                    sub_list_ct += 1
                    if sub_list_ct >= n_proc_chipName:
                        sub_list_ct = 0

            if n_proc_chipName>1:
                process_list = []
                for sub_list_ct in range(len(iobs_sub_list)):
                    p = mproc.Process(target=_find_chipNames_parallel,
                                      args=(chunk['raJ2000'], chunk['decJ2000']),
                                      kwargs={'pm_ra': pmra,
                                              'pm_dec': pmdec,
                                              'parallax': px,
                                              'v_rad': vrad,
                                              'obs_metadata_list': obs_sub_list[sub_list_ct],
                                              'i_obs_list': iobs_sub_list[sub_list_ct],
                                              'out_dict': chip_name_dict})
                    p.start()
                    process_list.append(p)

                for p in process_list:
                    p.join()

                assert len(chip_name_dict) == len(obs_valid)

            self._t_chip_name += time.time()-t_before_chip_name

            ######################################################
            # Calculate the delta_magnitude for all of the sources
            #
            t_before_phot = time.time()

            # only calculate photometry for objects that actually land
            # on LSST detectors

            valid_photometry = -1*np.ones(len(chunk))

            t_before_filter = time.time()
            for i_obs in range(len(obs_valid)):
                name_list, valid_obj = chip_name_dict[i_obs]
                valid_photometry[valid_obj] += 2
            invalid_dex = np.where(valid_photometry<0)
            chunk['varParamStr'][invalid_dex] = 'None'
            self._t_filter_phot += time.time()-t_before_filter

            photometry_catalog._set_current_chunk(chunk)
            dmag_arr = photometry_catalog.applyVariability(chunk['varParamStr'],
                                                           variability_cache=self._variability_cache,
                                                           expmjd=expmjd_list,).transpose((2,0,1))

            if np.abs(dmag_arr).max() < self._dmag_cutoff:
                continue

            self._t_phot += time.time()-t_before_phot

            ############################
            # Process and output sources
            #
            t_before_out = time.time()
            mag_name_to_int = {'u':0, 'g':1, 'r':2, 'i':3, 'z':4, 'y':5}
            for i_obs, obs in enumerate(obs_valid):
                obshistid = obs.OpsimMetaData['obsHistID']

                obs_mag = obs.bandpass
                actual_i_mag = mag_name_to_int[obs_mag]
                assert mag_names[actual_i_mag] == obs_mag

                # only include those sources which fall on a detector for this pointing
                valid_chip_name_list, valid_obj = chip_name_dict[i_obs]

                actually_valid_sources = np.where(np.abs(dmag_arr[i_obs][actual_i_mag][valid_obj]) >= self._dmag_cutoff)
                if len(actually_valid_sources[0]) == 0:
                    continue

                valid_sources = chunk[valid_obj]
                local_column_cache = {}
                local_column_cache['deltaMagAvro'] = OrderedDict([('delta_%smag' % mag_names[i_mag], dmag_arr[i_obs][i_mag][valid_obj])
                                                                  for i_mag in range(len(mag_names))])
                local_column_cache['chipName'] = valid_chip_name_list

                # only include those sources for which np.abs(delta_mag) >= self._dmag_cutoff
                # this is technically only selecting sources that differ from the quiescent
                # magnitude by at least self._dmag_cutoff.  If a source changes from quiescent_mag+dmag
                # to quiescent_mag, it will not make the cut
                valid_sources = valid_sources[actually_valid_sources]
                for mag in mag_names:
                    local_column_cache['deltaMagAvro']['delta_%smag' % mag] = \
                    local_column_cache['deltaMagAvro']['delta_%smag' % mag][actually_valid_sources]

                local_column_cache['chipName'] = local_column_cache['chipName'][actually_valid_sources]

                i_star = 0
                cat = cat_list[i_obs]
                for valid_chunk, chunk_map in cat.iter_catalog_chunks(query_cache=[valid_sources], column_cache=local_column_cache):

                    if obshistid not in output_data_cache:
                        output_data_cache[obshistid] = {}


                    data_tag = '%d_%d' % (obs.OpsimMetaData['obsHistID'], i_chunk)

                    for col_name in ('uniqueId', 'raICRS', 'decICRS', 'mag', 'mag_uncertainty', 'dmag', 'chipName', 'varParamStr'):
                        if col_name not in output_data_cache[obshistid]:
                            if col_name == 'chipName' or col_name == 'varParamStr':
                                output_data_cache[obshistid][col_name] = list(valid_chunk[chunk_map[col_name]].astype(str))
                            else:
                                output_data_cache[obshistid][col_name] = list(valid_chunk[chunk_map[col_name]])
                        else:
                            if col_name == 'chipName' or col_name == 'varParamStr':
                                output_data_cache[obshistid][col_name] += list(valid_chunk[chunk_map[col_name]].astype(str))
                            else:
                                output_data_cache[obshistid][col_name] += list(valid_chunk[chunk_map[col_name]])

                    ct_to_write += len(valid_chunk[chunk_map['uniqueId']])
                    print('ct_to_write %d' % ct_to_write)
                    if ct_to_write >= 10000:
                        # self.output_to_hdf5(out_file, output_data_cache)
                        ct_to_write = 0
                        output_data_cache = {}

                    #print star_obj
                #if i_chunk > 10:
                #    exit()
            self._t_out += time.time()-t_before_out

        # if len(output_data_cache)>0:
        #    self.output_to_hdf5(out_file, output_data_cache)

        # out_file.close()
        self._t_mlt += photometry_catalog._total_t_MLT
        self._t_param_lc += photometry_catalog._total_t_param_lc
        self._t_apply_var += photometry_catalog._total_t_apply_var
        return len(obs_valid)
