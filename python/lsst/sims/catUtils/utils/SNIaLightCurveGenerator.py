from __future__ import print_function
import numpy as np

from lsst.sims.catUtils.mixins import SNIaCatalog, PhotometryBase
from lsst.sims.catUtils.utils import _baseLightCurveCatalog
from lsst.sims.catUtils.utils import LightCurveGenerator

from lsst.sims.catUtils.supernovae import SNObject, SNUniverse
from lsst.sims.photUtils import PhotometricParameters, calcGamma
from lsst.sims.photUtils import Sed, calcMagError_m5, BandpassDict

import time

__all__ = ["SNIaLightCurveGenerator"]


class _sniaLightCurveCatalog(_baseLightCurveCatalog, SNIaCatalog, PhotometryBase):

    column_outputs = ["uniqueId", "snid", "raJ2000", "decJ2000",
                      "cosmologicalDistanceModulus", "redshift", "EBV"]

    _suppressDimSN = False


class SNIaLightCurveGenerator(LightCurveGenerator):
    """
    This class will find all of the OpSim pointings in a particular region
    of the sky in a particular filter and then return light curves for all
    of the supernovae observed in that region of sky.

    Input parameters:
    -----------------
    catalogdb is a CatalogDBObject instantiation connecting to the database
    of supernovae to be observed.

    opsimdb is the path to the OpSim database of observation.

    opsimdriver (optional; default 'sqlite') indicates the database driver to
    be used when connecting to opsimdb.
    """

    def __init__(self, *args, **kwargs):
        self.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()
        self._lightCurveCatalogClass = _sniaLightCurveCatalog
        self._filter_cat = None
        self._ax_cache = None
        self._bx_cache = None
        self._ax_bx_wavelen = None
        self.phot_params = PhotometricParameters()
        self.sn_universe = SNUniverse()
        self.sn_universe.suppressDimSN = False
        self.z_cutoff = 1.2
        super(SNIaLightCurveGenerator, self).__init__(*args, **kwargs)

    class _filterCatalogClass(_sniaLightCurveCatalog):
        column_outputs = ["uniqueId", "t0"]

    def _light_curves_from_query(self, cat_dict, query_result, grp):

        t_dict = {}
        gamma_dict = {}
        m5_dict = {}
        t_min = None
        t_max = None
        for bp_name in cat_dict:
            self.lsstBandpassDict[bp_name].sbTophi()

            # generate a 2-D numpy array containing MJDs, m5, and photometric gamma values
            # for each observation in the given bandpass
            raw_array = np.array([[obs.mjd.TAI, obs.m5[bp_name],
                                   calcGamma(self.lsstBandpassDict[bp_name],
                                             obs.m5[obs.bandpass],
                                             self.phot_params)]
                                  for obs in grp if obs.bandpass == bp_name]).transpose()

            if len(raw_array) > 0:

                t_dict[bp_name] = raw_array[0]

                m5_dict[bp_name] = raw_array[1]

                gamma_dict[bp_name] = raw_array[2]

                local_t_min = t_dict[bp_name].min()
                local_t_max = t_dict[bp_name].max()
                if t_min is None or local_t_min < t_min:
                    t_min = local_t_min

                if t_max is None or local_t_max > t_max:
                    t_max = local_t_max

        snobj = SNObject()

        cat = cat_dict[cat_dict.keys()[0]]  # does not need to be associated with a bandpass

        dummy_sed = Sed()

        for chunk in query_result:
            t_start_chunk = time.time()
            for sn in cat.iter_catalog(query_cache=[chunk]):
                sn_rng = self.sn_universe.getSN_rng(sn[1])
                sn_t0 = self.sn_universe.drawFromT0Dist(sn_rng)
                if sn[5] <= self.z_cutoff and np.isfinite(sn_t0) and \
                    sn_t0 < t_max + cat.maxTimeSNVisible and \
                    sn_t0 > t_min - cat.maxTimeSNVisible:

                    sn_c = self.sn_universe.drawFromcDist(sn_rng)
                    sn_x1 = self.sn_universe.drawFromx1Dist(sn_rng)
                    sn_x0 = self.sn_universe.drawFromX0Dist(sn_rng, sn_x1, sn_c, sn[4])

                    snobj.set(t0=sn_t0, c=sn_c, x1=sn_x1, x0=sn_x0, z=sn[5])

                    for bp_name in t_dict:
                        t_list = t_dict[bp_name]
                        m5_list = m5_dict[bp_name]
                        gamma_list = gamma_dict[bp_name]
                        bandpass = self.lsstBandpassDict[bp_name]
                        if len(t_list) == 0:
                            continue

                        if snobj.maxtime() >= t_list[0] and snobj.mintime() <= t_list[-1]:
                            active_dexes = np.where(np.logical_and(t_list >= snobj.mintime(),
                                                                   t_list <= snobj.maxtime()))

                            t_active = t_list[active_dexes]
                            m5_active = m5_list[active_dexes]
                            gamma_active = gamma_list[active_dexes]

                            if len(t_active) > 0:
                                wave_ang = bandpass.wavelen*10.0
                                mask = np.logical_and(wave_ang > snobj.minwave(),
                                                      wave_ang < snobj.maxwave())

                                wave_ang = wave_ang[mask]
                                snobj.set(mwebv=sn[6])
                                sn_ff_buffer = snobj.flux(time=t_active, wave=wave_ang)*10.0
                                flambda_grid = np.zeros((len(t_active), len(bandpass.wavelen)))
                                for ff, ff_sn in zip(flambda_grid, sn_ff_buffer):
                                    ff[mask] = np.where(ff_sn > 0.0, ff_sn, 0.0)

                                fnu_grid = flambda_grid*bandpass.wavelen* \
                                           bandpass.wavelen*dummy_sed._physParams.nm2m* \
                                           dummy_sed._physParams.ergsetc2jansky/dummy_sed._physParams.lightspeed

                                flux_list = \
                                (fnu_grid*bandpass.phi).sum(axis=1)*(bandpass.wavelen[1]-bandpass.wavelen[0])

                                mag_list = dummy_sed.magFromFlux(flux_list)

                                acceptable = np.where(np.isfinite(mag_list))
                                mag_error_list = calcMagError_m5(mag_list[acceptable], bandpass,
                                                                 m5_active[acceptable], self.phot_params,
                                                                 gamma=gamma_active[acceptable])

                                if len(acceptable) > 0:

                                    if sn[0] not in self.truth_dict:
                                        self.truth_dict[sn[0]] = {}
                                        self.truth_dict[sn[0]]['t0'] = sn_t0
                                        self.truth_dict[sn[0]]['x1'] = sn_x1
                                        self.truth_dict[sn[0]]['x0'] = sn_x0
                                        self.truth_dict[sn[0]]['c'] = sn_c
                                        self.truth_dict[sn[0]]['z'] = sn[5]
                                        self.truth_dict[sn[0]]['E(B-V)'] = sn[6]

                                    if sn[0] not in self.mjd_dict:
                                        self.mjd_dict[sn[0]] = {}
                                        self.bright_dict[sn[0]] = {}
                                        self.sig_dict[sn[0]] = {}

                                    if bp_name not in self.mjd_dict[sn[0]]:
                                        self.mjd_dict[sn[0]][bp_name] = []
                                        self.bright_dict[sn[0]][bp_name] = []
                                        self.sig_dict[sn[0]][bp_name] = []

                                for tt, mm, ee in zip(t_active[acceptable], mag_list[acceptable],
                                                      mag_error_list[0]):

                                    self.mjd_dict[sn[0]][bp_name].append(tt)
                                    self.bright_dict[sn[0]][bp_name].append(mm)
                                    self.sig_dict[sn[0]][bp_name].append(ee)

            print("chunk of ", len(chunk), " took ", time.time()-t_start_chunk)

