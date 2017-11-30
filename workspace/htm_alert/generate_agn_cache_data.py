import numpy as np
import os
import json
import multiprocessing as mproc
import argparse
import time

from lsst.sims.catUtils.mixins import ExtraGalacticVariabilityModels


def simulate_agn(galid_list, param_dict_list, output_dir, log_file, lock):
    evm = ExtraGalacticVariabilityModels()

    duration = 3654.0
    survey_start = 59580.0

    tot_steps = 0
    for param_dict in param_dict_list:
        tau = param_dict['pars']['agn_tau']
        dt = tau/100.0
        local_steps = int(np.ceil(duration/dt))
        tot_steps += local_steps

    t_start = time.time()
    n_simulated_steps = 0
    for galid, param_dict in zip(galid_list, param_dict_list):

        tau = param_dict['pars']['agn_tau']
        dt = tau/100.0
        mjd_arr = survey_start + np.arange(0.0, duration+5.0*dt, dt)

        params = {}
        for key in param_dict['pars']:
            params[key]=np.array([param_dict['pars'][key]])

        dmag_arr = evm.applyAgn(np.array([[0]]), params, mjd_arr)
        assert dmag_arr.shape == (6, 1, len(mjd_arr))

        out_name = os.path.join(output_dir,'agn_%d_lc.txt' % galid)
        with open(out_name, 'w') as out_file:
            out_file.write('# ')
            for mag_name in ('u', 'g', 'r', 'i', 'z', 'y'):
                out_file.write('%e ' % param_dict['pars']['agn_sf%s' % mag_name])
            out_file.write('\n')
            for i_time in range(len(mjd_arr)):
                out_file.write('%.6f ' % mjd_arr[i_time])
                for i_filter in range(6):
                    out_file.write('%.4f ' % (dmag_arr[i_filter][0][i_time]))
                out_file.write('\n')

        n_simulated_steps += len(mjd_arr)
        lock.acquire()
        with open(log_file, 'a') as out_file:
            elapsed = (time.time()-t_start)/3600.0
            projected =tot_steps*elapsed/n_simulated_steps
            out_file.write('simulated %d -- %d -- tau %e hrs elapsed %e proj %e\n' %
                           (galid, len(mjd_arr),tau,elapsed,projected))
        lock.release()


if __name__== "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--n_proc', type=int, default=4)
    parser.add_argument('--out_dir', type=str, default=None)
    parser.add_argument('--log_file', type=str, default=None)
    args= parser.parse_args()

    if args.out_dir is None:
        raise RuntimeError("must specify output dir")

    if not os.path.isdir(args.out_dir):
        raise RuntimeError('%s is not a dir' % args.out_dir)

    if args.log_file is None:
        raise RuntimeError("must specify log file")

    dtype = np.dtype([('galid', int), ('varParamStr', str, 300)])
    agn_data_file_name = 'agn_var_param_str.txt'
    if not os.path.exists(agn_data_file_name):
        raise RuntimeError('%s does not exist' % agn_data_file_name)

    data = np.genfromtxt(agn_data_file_name, delimiter=';',
                         dtype=dtype)


    galid_list = []
    param_dict_list = []
    steps_list = []
    for galid, varparamstr in zip(data['galid'], data['varParamStr']):
        galid_list.append(galid)
        param_dict = json.loads(varparamstr)
        param_dict_list.append(param_dict)
        dt = param_dict['pars']['agn_tau']/100.0
        steps_list.append(365.25/dt)

    steps_list = np.array(steps_list)
    galid_list = np.array(galid_list)
    param_dict_list = np.array(param_dict_list)

    sorted_dex = np.argsort(-1.0*steps_list)

    galid_list = galid_list[sorted_dex]
    param_dict_list = param_dict_list[sorted_dex]
    steps_list = steps_list[sorted_dex]

    n_steps = np.zeros(args.n_proc, dtype=int)
    galid_grid = []
    param_dict_grid = []
    for galid, param_dict, steps in zip(galid_list, param_dict_list, steps_list):
        if len(galid_grid) < args.n_proc:
            galid_grid.append([])
            galid_grid[-1].append(galid)
            param_dict_grid.append([])
            param_dict_grid[-1].append(param_dict)
            n_steps[len(galid_grid)-1] += steps
        else:
            min_dex = np.argmin(n_steps)
            galid_grid[min_dex].append(galid)
            param_dict_grid[min_dex].append(param_dict)
            n_steps[min_dex] += steps

    assert len(galid_grid) == args.n_proc
    assert len(param_dict_grid) == args.n_proc

    lock = mproc.Lock()
    for i_proc in range(args.n_proc):
        p = mproc.Process(target=simulate_agn,
                          args=(galid_grid[i_proc],
                                param_dict_grid[i_proc],
                                args.out_dir,
                                args.log_file,
                                lock))

        p.start()
