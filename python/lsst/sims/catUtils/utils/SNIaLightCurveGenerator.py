from __future__ import print_function
import numpy as np
import warnings

from lsst.sims.catUtils.mixins import SNIaCatalog, PhotometryBase
from lsst.sims.catUtils.utils import _baseLightCurveCatalog
from lsst.sims.catUtils.utils import LightCurveGenerator

from lsst.sims.catUtils.supernovae import SNObject, SNUniverse
from lsst.sims.photUtils import PhotometricParameters, calcGamma
from lsst.sims.photUtils import Sed, calcSNR_m5, BandpassDict

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

    Note: in this case, the method light_curves_from_pointings returns
    fluxes and flux errors in maggies.

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
        self._brightness_name = 'flux'
        super(SNIaLightCurveGenerator, self).__init__(*args, **kwargs)

    def light_curves_from_pointings(self, pointings, chunk_size=100000, lc_per_field=None):
        if lc_per_field is not None:
            warnings.warn("You have set lc_per_field in the SNIaLightCurveGenerator. "
                          "This will limit the number of candidate galaxies queried from the "
                          "CatSim database per field-of-view.  Because supernovae are randomly "
                          "populated in those galaxies, there is no guarantee that the galaxies "
                          "queried will have supernovae in them.  If the galaxies you actually "
                          "query do not host supernovae, you could get fewer light curves than "
                          "you expect.")

        return LightCurveGenerator.light_curves_from_pointings(self, pointings,
                                                               chunk_size=chunk_size,
                                                               lc_per_field=lc_per_field)

    def _get_query_from_group(self, grp, chunk_size, lc_per_field=None):
        """
        Override _get_query_from_group.  The probabilistic nature of SNe requires
        that we always actually do the query with lc_per_field=None (since we can't be
        guaranteed that any given galaxy, though it contains a SN, will have an SN that
        is going off during our time of interest).
        """
        return LightCurveGenerator._get_query_from_group(self, grp, chunk_size, lc_per_field=None)

    class _filterCatalogClass(_sniaLightCurveCatalog):
        column_outputs = ["uniqueId", "t0"]

    def _light_curves_from_query(self, cat_dict, query_result, grp, lc_per_field=None):

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

        n_actual_sn = 0  # how many SN have we actually delivered?

        for chunk in query_result:

            if lc_per_field is not None and n_actual_sn >= lc_per_field:
                break

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

                                acceptable = np.where(flux_list>0.0)

                                flux_error_list = flux_list[acceptable]/ \
                                                  calcSNR_m5(dummy_sed.magFromFlux(flux_list[acceptable]),
                                                             bandpass,
                                                             m5_active[acceptable], self.phot_params,
                                                             gamma=gamma_active[acceptable])

                                if len(acceptable) > 0:

                                    n_actual_sn += 1
                                    if lc_per_field is not None and n_actual_sn > lc_per_field:
                                        break

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

                                for tt, ff, ee in zip(t_active[acceptable], flux_list[acceptable],
                                                      flux_error_list[0]):

                                    self.mjd_dict[sn[0]][bp_name].append(tt)
                                    self.bright_dict[sn[0]][bp_name].append(ff/3631.0)
                                    self.sig_dict[sn[0]][bp_name].append(ee/3631.0)

            print("chunk of ", len(chunk), " took ", time.time()-t_start_chunk)

