import numpy as np
from lsst.sims.utils import findHtmid, trixelFromHtmid
from lsst.sims.utils import angularSeparation, ObservationMetaData
from lsst.sims.catUtils.utils import _baseLightCurveCatalog
from lsst.sims.coordUtils import _chipNameFromRaDecLSST

from lsst.sims.catalogs.decorators import compound
from lsst.sims.photUtils import BandpassDict
from lsst.sims.catUtils.mixins import VariabilityStars, AstrometryStars
from lsst.sims.catUtils.mixins import PhotometryStars, CameraCoordsLSST

__all__ = ["AvroGenerator"]

class StellarVariabilityCatalog(VariabilityStars, PhotometryStars, AstrometryStars,
                                CameraCoordsLSST, _baseLightCurveCatalog):
    column_outputs = ['uniqueId', 'raICRS', 'decICRS',
                      'mag','mag_uncertainty', 'chipName']

    @compound('sigma_native_u', 'sigma_native_g', 'sigma_native_r',
              'sigma_native_i', 'sigma_native_z', 'sigma_native_y')
    def get_avroPhotometricUncertainty(self):
        if not hasattr(self, '_avro_bp_dict'):
            self._avro_bp_dict = BandpassDict.loadTotalBandpassesFromFiles()
        return self._magnitudeUncertaintyGetter(['native_u', 'native_g', 'native_r',
                                                 'native_i', 'native_z', 'native_y'],
                                                ['u', 'g', 'r', 'i', 'z', 'y'],
                                                '_avro_bp_dict')

    @compound('native_u','native_g','native_r','native_i','native_z','native_y')
    def get_native_magnitudes(self):
        """
        getter for LSST stellar magnitudes
        """

        magnitudes = np.array([self.column_by_name('umag'),
                               self.column_by_name('gmag'),
                               self.column_by_name('rmag'),
                               self.column_by_name('imag'),
                               self.column_by_name('zmag'),
                               self.column_by_name('ymag')])

        delta = self._variabilityGetter(self.get_lsst_magnitudes._colnames)
        magnitudes += delta

        return magnitudes

    @compound('mag', 'mag_uncertainty')
    def get_avroPhotometry(self):
        mag = self.column_by_name('native_%s' % self.obs_metadata.bandpass)
        mag_unc = self.column_by_name('sigma_native_%s' % self.obs_metadata.bandpass)
        return np.array([mag, mag_unc])


class AvroGenerator(object):

    def __init__(self, obs_list):
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
        for htmid in self.unq_htmid_list:
            print('processing %d' % htmid)
            self._process_htmid(htmid, dbobj)

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
            cat._avro_bp_dict =  self.bp_dict
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

        for chunk in data_iter:
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
                    print star_obj

