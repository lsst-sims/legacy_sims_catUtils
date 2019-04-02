import numpy as np
import os

from lsst.sims.utils import radiansFromArcsec
from lsst.sims.catUtils.mixins import MLTflaringMixin
from lsst.sims.photUtils import PhotometricParameters
from lsst.sims.photUtils import BandpassDict
from lsst.sims.photUtils import Sed

quiet_mag = 17.0
dummy_sed = Sed()

phot_params = PhotometricParameters(nexp=1, exptime=30)

mlt_model = MLTflaringMixin()
mlt_model.photParams = phot_params
mlt_model.lsstBandpassDict = BandpassDict.loadTotalBandpassesFromFiles()
mlt_model._actually_calculated_columns = list(['lsst_%s' % bp
                                               for bp in 'ugrizy'])

fname = os.path.join(os.environ['SIMS_DATA_DIR'],'catUtilsData',
                     'mlt_shortened_lc_171012.npz')

assert os.path.isfile(fname)

lc_curve_time = {}
t_max_dex_dict = {}
with np.load(fname) as lc_file:
    for name in lc_file.keys():
        if 'time' in name:
            lc_curve_time[name] = lc_file[name]
        else:
            dex = np.argmax(lc_file[name])
            t_max_dex_dict[name] = dex

lc_name_list = ['early_active', 'early_inactive',
                'mid_active', 'mid_inactive',
                'late_active']

milliarcsec = radiansFromArcsec(0.001)
ebv_grid = np.arange(0.01, 7.01, 0.01)
parallax_grid = np.arange(0.01*milliarcsec,50.0*milliarcsec, milliarcsec)

n_steps = len(parallax_grid)

out_cache = {}
out_cache['parallax_grid(mas)'] = parallax_grid/milliarcsec

for lc_name_root in lc_name_list:
    for lc_dex in range(4):
        lc_name = lc_name_root + '_%d' % lc_dex
        print(lc_name)
        t_max_dex = t_max_dex_dict[lc_name+'_u']

        for bp in 'grizy':
            t = t_max_dex_dict[lc_name+'_%s' % bp]
            assert t == t_max_dex

        t_max = lc_curve_time[lc_name+'_time'][t_max_dex]+mlt_model._survey_start

        params = {}
        params['lc'] = np.array([lc_name]*n_steps)
        params['t0'] = np.array([0.0]*n_steps)

        q_mags = {}
        for bp in 'ugrizy':
            q_mags[bp] = quiet_mag*np.ones(n_steps, dtype=float)

        for ebv_val in ebv_grid:
            ebv = np.array([ebv_val]*n_steps, dtype=float)

            dmag_base = mlt_model.applyMLTflaring(np.where(np.array(n_steps*[True])),
                                                  params, t_max,
                                                  parallax=parallax_grid, ebv=ebv,
                                                  quiescent_mags = q_mags)

            mag_base = dmag_base+quiet_mag
            flux1 = dummy_sed.fluxFromMag(mag_base)
            flux0 = dummy_sed.fluxFromMag(quiet_mag)
            dflux = flux1-flux0
            for ii, bp in enumerate('ugrizy'):
                dflux_name = '%s_ebv_%.2f_%s' % (lc_name, ebv_val, bp)
                assert dflux_name not in out_cache
                out_cache[dflux_name] = dflux[ii]

np.savez('mlt_dflux_lookup.npz', **out_cache)
