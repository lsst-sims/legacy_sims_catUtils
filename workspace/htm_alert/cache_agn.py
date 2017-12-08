import h5py
import numpy as np
import os
import argparse
import time

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--in_dir', type=str, default=None)
    parser.add_argument('--out_dir', type=str, default=None)
    parser.add_argument('--t_step', type=int, default=365,
                        help="The length of time (in days) in which to "
                        "group light curves")


    log_file_name = 'agn_cache_log.txt'
    if os.path.exists(log_file_name):
        os.unlink(log_file_name)

    args = parser.parse_args()
    if args.in_dir is None:
        raise RuntimeError("you must specify in_dir")
    if args.out_dir is None:
        raise RuntimeError("you must specify out_dir")

    print('cacheing AGN in pid %d' % os.getpid())

    t_start_arr = np.arange(59580.0, 59580.0+3654.0, args.t_step)

    lc_dtype = np.dtype([('mjd', float), ('du', float), ('dg', float),
                         ('dr', float), ('di', float), ('dz', float)])

    list_of_all_files = os.listdir(args.in_dir)
    t_start = time.time()
    created_file_names = set()
    for i_file, file_name in enumerate(list_of_all_files):
        if file_name.startswith('agn') and file_name.endswith('lc.txt'):
            full_name = os.path.join(args.in_dir, file_name)
            galid = int(file_name.split('_')[1])
            with open(full_name, 'r') as in_file:
                first_line = in_file.readline().strip()
                params = first_line.split()
                assert params[0] == '#'
                normalizing_factors = []
                for i_filter in range(6):
                    normalizing_factors.append(float(params[i_filter+1]))

            lc_data = np.genfromtxt(full_name, dtype=lc_dtype)

            duration = lc_data['mjd'][-1] - lc_data['mjd'][0]
            dt = lc_data['mjd'][1]-lc_data['mjd'][0]
            for i_t, t_start in enumerate(t_start_arr):
                t_name = int(np.floor(t_start))
                out_name = os.path.join(args.out_dir, 'agn_lcs_%d.hdf5' % t_name)
                mjd_min = float(t_name)-5.0*dt
                mjd_max = t_start+args.t_step+5.0*dt
                if i_t == 0:
                    valid = np.where(lc_data['mjd']<mjd_max)
                elif i_t == (len(t_start_arr)-1):
                    valid = np.where(lc_data['mjd']>mjd_min)
                else:
                    valid = np.where(np.logical_and(lc_data['mjd']>mjd_min,
                                                    lc_data['mjd']<mjd_max))

                if out_name not in created_file_names:
                    if os.path.exists(out_name):
                        raise RuntimeError('%s exists already' % out_name)
                    created_file_names.add(out_name)

                with h5py.File(out_name, 'a') as out_file:
                    out_file.create_dataset('%d_norm' % galid, data=normalizing_factors)
                    out_file.create_dataset('%d_mjd' % galid, data=lc_data['mjd'][valid])
                    out_file.create_dataset('%d_du' % galid, data=lc_data['du'][valid])

        with open(log_file_name,'a') as out_file:
            elapsed = (time.time()-t_start)/3600.0
            out_file.write('did %d of %d in %.2e hours\n' %
                           (i_file, len(list_of_all_files), elapsed))
    print('all done')
