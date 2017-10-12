import numpy as np
import time
from lsst.sims.utils import findHtmid, trixelFromHtmid
from lsst.sims.utils import angularSeparation, ObservationMetaData
from lsst.sims.catUtils.utils import _baseLightCurveCatalog
from lsst.sims.coordUtils import _chipNameFromRaDecLSST

from lsst.sims.catalogs.decorators import compound
from lsst.sims.photUtils import BandpassDict
from lsst.sims.catUtils.mixins import VariabilityStars, AstrometryStars
from lsst.sims.catUtils.mixins import CameraCoordsLSST, PhotometryBase
from lsst.sims.catUtils.mixins import ParametrizedLightCurveMixin

__all__ = ["AvroGenerator"]

class StellarVariabilityCatalog(VariabilityStars, AstrometryStars, PhotometryBase,
                                CameraCoordsLSST, _baseLightCurveCatalog):
    column_outputs = ['uniqueId', 'raICRS', 'decICRS',
                      'mag','mag_uncertainty', 'dmag', 'chipName']

    default_formats = {'f':'%.4g'}

    _d_ct = 0

    @compound('sigma_lsst_u', 'sigma_lsst_g', 'sigma_lsst_r',
              'sigma_lsst_i', 'sigma_lsst_z', 'sigma_lsst_y')
    def get_lsst_photometric_uncertainties(self):
        if not hasattr(self, 'lsstBandpassDict'):
            self.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()
        return self._magnitudeUncertaintyGetter(['lsst_u', 'lsst_g', 'lsst_r',
                                                 'lsst_i', 'lsst_z', 'lsst_y'],
                                                ['u', 'g', 'r', 'i', 'z', 'y'],
                                                'lsstBandpassDict')

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

        delta = self._variabilityGetter(self.get_lsst_magnitudes._colnames)
        magnitudes += delta

        return magnitudes

    @compound('mag', 'mag_uncertainty', 'dmag')
    def get_avroPhotometry(self):
        mag = self.column_by_name('lsst_%s' % self.obs_metadata.bandpass)
        mag_unc = self.column_by_name('sigma_lsst_%s' % self.obs_metadata.bandpass)
        dmag = self.column_by_name('%smag' % self.obs_metadata.bandpass) - mag

        return np.array([mag, mag_unc, dmag])


class AvroGenerator(object):

    def __init__(self, obs_list):
        plm = ParametrizedLightCurveMixin()
        plm.load_parametrized_light_curves()
        self.obs_list = np.array(obs_list)
        htmid_level = 7
        self.htmid_list = []
        for obs in self.obs_list:
            htmid = findHtmid(obs.pointingRA, obs.pointingDec, htmid_level)
            self.htmid_list.append(htmid)
        self.htmid_list = np.array(self.htmid_list)
        self.unq_htmid_list = np.unique(self.htmid_list)
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
              (len(self.obs_list), len(self.unq_htmid_list)))

    def alerts_from_db(self, dbobj):
        for i_h, htmid in enumerate(self.unq_htmid_list):
            print('processing %d --- %d of %d' % (htmid, i_h, len(self.unq_htmid_list)))
            t_start = time.time()
            self._process_htmid(htmid, dbobj)
            print("that took %e hours" % ((time.time()-t_start)/3600.0))


    def _process_htmid(self, htmid, dbobj, radius=1.75):
        valid_dexes = np.where(self.htmid_list == htmid)
        print('valid_dexes %s ' % str(valid_dexes))
        obs_valid = self.obs_list[valid_dexes]
        center_trixel = trixelFromHtmid(htmid)
        center_ra, center_dec = center_trixel.get_center()

        ra_list = []
        dec_list = []
        cat_list = []
        expmjd_list = []
        for obs in obs_valid:
            ra_list.append(obs.pointingRA)
            dec_list.append(obs.pointingDec)
            cat = StellarVariabilityCatalog(dbobj, obs_metadata=obs)
            cat.lsstBandpassDict =  self.bp_dict
            cat_list.append(cat)
            expmjd_list.append(obs.mjd.TAI)

        expmjd_list = np.array(expmjd_list)
        sorted_dex = np.argsort(expmjd_list)

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

        print('chunking')
        i_chunk = 0
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

            for i_obs, obs in enumerate(obs_valid):
                chip_name_list = _chipNameFromRaDecLSST(chunk['raJ2000'],
                                                        chunk['decJ2000'],
                                                        pm_ra=pmra,
                                                        pm_dec=pmdec,
                                                        parallax=px,
                                                        v_rad=vrad,
                                                        obs_metadata=obs)

                valid = np.where(np.char.find(chip_name_list.astype(str), 'R')==0)
                for name in chip_name_list[valid]:
                    assert name is not None
                valid_sources = chunk[valid]
                cat = cat_list[i_obs]
                i_star = 0
                for star_obj in cat.iter_catalog(query_cache=[valid_sources]):
                    pass
                    #print star_obj
                if i_chunk > 4:
                    exit()
