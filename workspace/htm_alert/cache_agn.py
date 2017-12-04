import numpy as np
import os
import argparse
import time

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--in_dir', type=str, default=None)
    parser.add_argument('--out_dir', type=str, default=None)
    parser.add_argument('--t_step', type=float, default=365.25,
                        help="The length of time (in days) in which to "
                        "group light curves")


    args = parser.parse_args()
    if args.in_dir is None:
        raise RuntimeError("you must specify in_dir")
    if args.out_dir is None:
        raise RuntimeError("you must specify out_dir")

    print('cacheing AGN in pid %d' % os.getpid())

    lc_dtype = np.dtype([('mjd', float), ('du', float), ('dg', float),
                         ('dr', float), ('di', float), ('dz', float)])

    galid_list = []
    normalizing_factors = {}
    lc_data = {}
    list_of_all_files = os.listdir(args.in_dir)
    for file_name in list_of_all_files:
        if file_name.startswith('agn') and file_name.endswith('lc.txt'):
            full_name = os.path.join(args.in_dir, file_name)
            galid = int(file_name.split('_')[1])
            galid_list.append(galid)
            with open(full_name, 'r') as in_file:
                first_line = in_file.readline().strip()
                params = first_line.split()
                assert params[0] == '#'
                normalizing_factors[galid] = []
                for i_filter in range(6):
                    normalizing_factors[galid].append(float(params[i_filter+1]))

            raw_data = np.genfromtxt(full_name, dtype=lc_dtype)
            lc_data[galid] = {}
            lc_data[galid]['mjd'] = raw_data['mjd']
            lc_data[galid]['du'] = raw_data['du']

    t_start = 59580.0
    duration = 3654.0
    n_t_steps = np.floor(duration/args.t_step).astype(int)
    print('read in all agn')
    time.sleep(30)
    exit(1)

    for i_t, t_min in enumerate(np.arange(t_start, t_start+duration, args.t_step)):
        t_max = t_min + args.t_step
        output_lc = {}
        print('time range %e %e' % (t_min, t_max))

        for galid, file_name in zip(galid_list, lc_file_list):
            data = np.genfromtxt(file_name, dtype=lc_dtype)
            dt = data['mjd'][1] - data['mjd'][0]
            valid = np.where(np.logical_and(data['mjd']>t_min-10.0*dt,
                                            data['mjd']<t_max+10.0*dt))

            mjd = data['mjd'][valid]
            du = data['du'][valid]
            output_lc['%d_mjd' % galid] = mjd
            output_lc['%d_du' % galid] = du
            output_lc['%d_norms' % galid] = normalizing_factors[galid]

        out_file_name = os.path.join(args.out_dir, 'agn_lc_cache_%.1f_%.1f.npz' % (t_min, t_max))
        with open(out_file_name, 'wb') as file_handle:
            np.savez(file_handle, **output_lc)
