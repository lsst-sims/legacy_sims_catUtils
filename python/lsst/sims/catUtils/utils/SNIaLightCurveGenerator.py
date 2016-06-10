from __future__ import print_function
import numpy as np

from lsst.sims.catalogs.measures.instance import compound
from lsst.sims.catUtils.mixins import SNIaCatalog, PhotometryBase
from lsst.sims.catUtils.utils import _baseLightCurveCatalog
from lsst.sims.catUtils.utils import LightCurveGenerator

from lsst.sims.catUtils.supernovae import SNObject
from lsst.sims.photUtils import PhotometricParameters, calcGamma
from lsst.sims.photUtils import Sed, calcMagError_m5

__all__ = ["SNIaLightCurveGenerator"]



class _sniaLightCurveCatalog(_baseLightCurveCatalog, SNIaCatalog, PhotometryBase):

    column_outputs = ["uniqueId", "raJ2000", "decJ2000",
                      "t0", "c", "x1", "x0", "redshift", "EBV"]

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
        self._ax_cache = None
        self._bx_cache = None
        self._ax_bx_wavelen = None
        self.phot_params = PhotometricParameters()
        super(SNIaLightCurveGenerator, self).__init__(*args, **kwargs)


    class _filterCatalogClass(_sniaLightCurveCatalog):
        column_outputs = ["uniqueId", "t0"]


    def _light_curves_from_query(self, cat, query_result, grp):

        bandpass = cat.lsstBandpassDict[grp[0].bandpass]

        t_list = np.array([obs.mjd.TAI for obs in grp])

        gamma_list = np.array([calcGamma(bandpass, obs.m5[obs.bandpass], self.phot_params)
                               for obs in grp])

        for chunk in query_result:
            for sn in cat.iter_catalog(query_cache=[chunk]):
                if np.isfinite(sn[3]) and \
                sn[3]<t_list[-1]+cat.maxTimeSNVisible and \
                sn[3]>t_list[0]-cat.maxTimeSNVisible:

                    snobj = SNObject()
                    snobj.set(t0=sn[3], c=sn[4], x1=sn[5], x0=sn[6], z=sn[7])

                    if snobj.maxtime()>=t_list[0] and snobj.mintime()<=t_list[-1]:
                        active_dexes = np.where(np.logical_and(t_list>=snobj.mintime(),
                                                               t_list<=snobj.maxtime()))

                        t_active = t_list[active_dexes]
                        m5_active = np.array([obs.m5[obs.bandpass] for obs in grp])[active_dexes]
                        gamma_active = gamma_list[active_dexes]

                        if len(t_active)>0:
                            wave_ang = bandpass.wavelen * 10.0
                            mask1 = wave_ang > snobj.minwave()
                            mask2 = wave_ang < snobj.maxwave()
                            mask = mask1 & mask2
                            wave_ang = wave_ang[mask]
                            flambda_grid = snobj.flux(time=t_active, wave=wave_ang)*10.0

                            flambda = np.zeros(len(bandpass.wavelen))*np.nan
                            mag_list = []
                            for ix, ff in enumerate(flambda_grid):
                                flambda[mask] = ff
                                flambda = np.where(flambda > 0., flambda, 0.)
                                ss = Sed(wavelen=bandpass.wavelen, flambda=flambda)

                                if self._ax_bx_wavelen is None or \
                                len(self._ax_bx_wavelen) != len(bandpass.wavelen) or \
                                (self._ax_bx_wavelen!=bandpass.wavelen).any():

                                    self._ax_bx_wavelen = bandpass.wavelen
                                    self._ax_cache, self._bx_cache = ss.setupCCMab()

                                ss.addCCMDust(a_x=self._ax_cache, b_x=self._bx_cache, ebv=sn[8])
                                flux = ss.calcFlux(bandpass)
                                mag = ss.magFromFlux(flux)
                                mag_list.append(mag)

                            mag_list = np.array(mag_list)

                            acceptable = np.where(np.isfinite(mag_list))
                            mag_error_list = calcMagError_m5(mag_list[acceptable], bandpass,
                                                             m5_active[acceptable], self.phot_params,
                                                             gamma=gamma_active[acceptable])

                            if sn[0] not in self.mjd_dict:
                                self.mjd_dict[sn[0]] = []
                                self.mag_dict[sn[0]] = []
                                self.sig_dict[sn[0]] = []

                            for tt, mm, ee in zip(t_active[acceptable], mag_list[acceptable],
                                                  mag_error_list[0]):

                                self.mjd_dict[sn[0]].append(tt)
                                self.mag_dict[sn[0]].append(mm)
                                self.sig_dict[sn[0]].append(ee)
