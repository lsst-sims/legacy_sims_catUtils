from __future__ import print_function
import numpy as np

from lsst.sims.catUtils.utils import ObservationMetaDataGenerator
from  lsst.sims.catUtils.mixins import PhotometryStars, VariabilityStars
from lsst.sims.catalogs.measures.instance import InstanceCatalog, compound

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

        for ix, obs in enumerate(obs_list):
            print('%d of %d at %e' % (ix, len(obs_list), time.clock()-t_start))
            cat = _lightCurveCatalog(self._catalogdb, obs_metadata=obs)
            for star_obj in cat.iter_catalog():
                    if star_obj[0] not in output_dict:
                        output_dict[star_obj[0]] = []

                    output_dict[star_obj[0]].append((obs.mjd.TAI,
                                                     star_obj[3],
                                                     star_obj[4]))

        print('that took %e' % (time.clock()-t_start))
        return output_dict
