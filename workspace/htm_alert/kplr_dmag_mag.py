import numpy as np
import time
import os
import multiprocessing as mpc
from lsst.sims.catUtils.mixins import ParametrizedLightCurveMixin
from lsst.sims.catUtils.mixins import create_variability_cache

def get_dmag(lc_id_list, out_file_name, seed, lock, v_cache):
    rng = np.random.RandomState(seed)
    plc = ParametrizedLightCurveMixin()
    v_cache = create_variability_cache()
    plc.load_parametrized_light_curves(variability_cache=v_cache)


    dt = (60.0*10.0)/(3600.0*24.0)
    mjd_base = np.arange(0.0, 180.0+2.0*dt, dt)
    dmag_list = np.zeros(len(lc_id_list), dtype=float)
    for i_lc, lc_int in enumerate(lc_id_list):
        dmag_max = -1.0
        for mjd_start in np.arange(59580.0, 59580.0+3654.0+181.0,180.0):
            mjd_arr = mjd_base+mjd_start
            d_mjd = rng.random_sample(len(mjd_arr))*0.9*dt
            mjd_arr += d_mjd
            q_flux, d_flux = plc._calc_dflux(lc_int, mjd_arr,
                                             variability_cache=v_cache)

            d_flux_max = np.abs(d_flux/q_flux).max()
            dmag_max_local = 2.5*np.log10(1.0+d_flux_max)

            if dmag_max_local>dmag_max:
                dmag_max = dmag_max_local

        dmag_list[i_lc] = dmag_max

    lock.acquire()
    with open(out_file_name, 'a') as out_file:
        for lc_id, dmag in zip(lc_id_list, dmag_list):
            out_file.write('%d %e' % (lc_id, dmag))
    lock.release()

if __name__ == "__main__":


    out_dir = os.path.join('/astro', 'store', 'pogo4', 'danielsf', 'kplr_dmag')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    if not os.path.isdir(out_dir):
        raise RuntimeError('%s is not a dir' % out_dir)

    out_file_name = os.path.join(out_dir,'kpr_dmag_171204.txt')
    if os.path.exists(out_file_name):
        raise RuntimeError('%s exists' % out_file_name)

    n_proc = 10

    plc = ParametrizedLightCurveMixin()
    v_cache = create_variability_cache()
    plc.load_parametrized_light_curves(variability_cache=v_cache)

    lc_id_list = []
    for i_p in range(n_proc):
        lc_id_list.append([])
    i_p = 0
    for lc_id in v_cache['_PARAMETRIZED_LC_MODELS']:
        lc_id_list[i_p].append(lc_id)
        i_p += 1
        if i_p>=n_proc:
            i_p = 0

    lock = mpc.Lock()
    p_list = []
    for i_p in range(n_proc):
        p = mpc.Process(target=get_dmag,
                        args=[lc_id_list[i_p],, out_file_name,
                              101+i_p, lock, v_cache])
        p.start()
        p_list.append(p)

    for p in p_list:
        p.join()

    print('all done')
