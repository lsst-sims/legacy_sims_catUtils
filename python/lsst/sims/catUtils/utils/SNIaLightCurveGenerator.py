from __future__ import print_function
import numpy as np

from lsst.sims.catalogs.measures.instance import compound
from lsst.sims.catUtils.mixins import SNIaCatalog, PhotometryBase
from lsst.sims.catUtils.utils import _baseLightCurveCatalog
from lsst.sims.catUtils.utils import LightCurveGenerator


__all__ = ["SNIaLightCurveGenerator"]



class _sniaLightCurveCatalog(_baseLightCurveCatalog, SNIaCatalog, PhotometryBase):

    _suppressDimSN = False
    _sn_params_id_cache = None
    _sn_params_cache = None


    @compound('c', 'x1', 'x0', 't0')
    def get_snparams(self):
        id_list = self.column_by_name('uniqueId')

        if self._sn_params_id_cache is None or \
        len(id_list)!=len(self._sn_params_id_cache) or \
        (id_list!=self._sn_params_id_cache).any():

            print('generating params')

            vals = self.SNparamDistFromHost(self.column_by_name('redshift'),
                                            self.column_by_name('snid'),
                                            self.column_by_name('cosmologicalDistanceModulus'))

            self._sn_params_id_cache = id_list
            self._sn_params_cache = np.array([vals[:, 0], vals[:, 1], vals[:, 2], vals[:, 3]])

        return self._sn_params_cache


    @compound("lightCurveMag", "sigma_lightCurveMag")
    def get_lightCurvePhotometry(self):

        return np.array([self.column_by_name("mag"),
                         self.column_by_name("mag_err")])





class SNIaLightCurveGenerator(LightCurveGenerator):

    def __init__(self, *args, **kwargs):
        self._lightCurveCatalogClass = _sniaLightCurveCatalog
        self._filter_cat = None
        super(SNIaLightCurveGenerator, self).__init__(*args, **kwargs)


    class _filterCatalogClass(_sniaLightCurveCatalog):
        column_outputs = ["uniqueId", "t0"]

    def _filter_chunk(self, raw_chunk):

        if self._filter_cat is None:
            self._filter_cat = self._filterCatalogClass(self._catalogdb)
            self._filter_cat.suppressDimSN = False
            db_required_columns = self._filter_cat.db_required_columns()

        valid_chunk = []
        for star_obj, raw_obj in zip(self._filter_cat.iter_catalog(query_cache=[raw_chunk]), raw_chunk):
            if not np.isnan(star_obj[1]) and \
            star_obj[1] > self._mjd_min-self._filter_cat.maxTimeSNVisible and \
            star_obj[1] < self._mjd_max+self._filter_cat.maxTimeSNVisible:

                valid_chunk.append(raw_obj)

        if len(valid_chunk)==0:
            return None

        return np.core.records.fromrecords(valid_chunk, dtype=raw_chunk.dtype)

