from __future__ import print_function
import numpy as np

from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from  lsst.sims.catUtils.mixins import PhotometryStars, VariabilityStars
from lsst.sims.catalogs.measures.instance import InstanceCatalog, compound
from lsst.sims.utils import haversine

import time

class _lightCurveCatalog(InstanceCatalog, VariabilityStars, PhotometryStars):

    column_outputs = ["uniqueId", "raJ2000", "decJ2000",
                      "lightCurveMag", "sigma_lightCurveMag"]

    @compound("lightCurveMag", "sigma_lightCurveMag")
    def get_lightCurvePhotometry(self):
        if len(self.obs_metadata.bandpass) != 1:
            raise RuntimeError("_lightCurveCatalog cannot handle bandpass "
                               "%s" % str(self.obs_metadata.bandpass))

        return np.array([
                 self.column_by_name("lsst_%s" % self.obs_metadata.bandpass),
                 self.column_by_name("sigma_lsst_%s" % self.obs_metadata.bandpass)])


class _chunkIterWrapper(object):

    def __init__(self, data):
        self._data = data
        self._ix = 0

    def __iter__(self):
        return self

    def next(self):
        if self._ix < len(self._data):
            self._ix += 1
            return self._data[self._ix-1]
        else:
            raise StopIteration


class LightCurveGenerator(object):

    def __init__(self, catalogdb, opsimdb, opsimdriver="sqlite"):
        self._generator = ObservationMetaDataGenerator(database=opsimdb,
                                                       driver=opsimdriver)

        self._catalogdb = catalogdb


    def generate_light_curves(self, ra, dec, bandpass):
        obs_list = self._generator.getObservationMetaData(
                                     fieldRA=ra,
                                     fieldDec=dec,
                                     telescopeFilter=bandpass,
                                     boundLength=1.75)

        output_dict = {}

        if len(obs_list) == 0:
            print("No observations found matching your criterion")
            return None

        t_start = time.clock()
        print('starting light curve generation')

        tol = 1.0e-12
        obs_groups = []
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

        cat = None

        for grp in obs_groups:
            if len(grp) == 1:
                obs = obs_list[grp[0]]
                if cat is None:
                    cat = _lightCurveCatalog(self._catalogdb, obs_metadata=obs)
                else:
                    cat.obs_metadata = obs
                    cat._gamma_cache = {}

                for star_obj in cat.iter_catalog():
                    if star_obj[0] not in output_dict:
                        output_dict[star_obj[0]] = []

                    output_dict[star_obj[0]].append((obs.mjd.TAI,
                                                     star_obj[3],
                                                     star_obj[4]))

            else:
                dataCache = None
                for ix in grp:
                    obs = obs_list[ix]
                    if cat is None:
                        cat = _lightCurveCatalog(self._catalogdb, obs_metadata=obs)
                        db_required_columns = cat.db_required_columns()

                    if dataCache is None:
                        query_result = cat.db_obj.query_columns(colnames=cat._active_columns,
                                                               obs_metadata=obs,
                                                               constraint=cat.constraint,
                                                               chunk_size=100000)
                        dataCache = []
                        for chunk in query_result:
                           dataCache.append(chunk)

                    query_result = _chunkIterWrapper(dataCache)
                    cat.obs_metadata = obs
                    cat._gamma_cache = {}

                    for star_obj in cat.iter_catalog(query_result=query_result):
                        if star_obj[0] not in output_dict:
                            output_dict[star_obj[0]] = []

                        output_dict[star_obj[0]].append((obs.mjd.TAI,
                                                         star_obj[3],
                                                         star_obj[4]))

        print('that took %e' % (time.clock()-t_start))
        return output_dict
