import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
import os
import pickle

def get_rate(fname, cut=None):
    assert os.path.isfile(fname)

    with open(fname, 'rb') as in_file:
        alert_dict = pickle.load(in_file)

    mjd_arr = -1.0*np.ones(len(alert_dict), dtype=float)
    for ii, name in enumerate(alert_dict):
        if cut is not None:
            if alert_dict[name][2] != cut:
                continue
        mjd_arr[ii] = alert_dict[name][0]

    del alert_dict
    valid = np.where(mjd_arr>0.0)
    return mjd_arr[valid]

    mjd_arr = np.round(mjd_arr).astype(int)

    raw_date_arr, raw_ct_arr = np.unique(mjd_arr, return_counts=True)

    date_arr = np.arange(59580, 59580+3653, dtype=int)
    ct_arr = np.zeros(len(date_arr), dtype=int)
    for ii, date in enumerate(date_arr):
        if date in raw_date_arr:
            dex = np.where(raw_date_arr==date)
            ct_arr[ii] = raw_ct_arr[dex]

    print('processed %s' % fname)
    return date_arr, ct_arr

if __name__ == "__main__":
    region_list = ['region_1', 'region_2', 'region_3', 'region_4']

    region_dict={}
    region_dict['region_1'] = 'ecliptic; 9.7e4 visits'
    region_dict['region_2'] = 'south galactic pole; 9.7e4 visits'
    region_dict['region_3'] = 'galactic center; 3.4e4 visits'
    region_dict['region_4'] = 'galactic anti-center; 450 deg^2; 2.1e4 visits'

    plt.figure(figsize=(20,20))
    for i_reg, region_name in enumerate(region_list):
        plt.subplot(2,2,i_reg+1)
        star_file = '%s_stars.pickle' % region_name
        assert os.path.isfile(star_file)
        agn_file = '%s_agn.pickle' % region_name
        assert os.path.isfile(agn_file)
        sne_file = '%s_sne.pickle' % region_name
        assert os.path.isfile(sne_file)


        handle_list = []
        name_list = []      
        for i_fig, (file_name, color, pop, cut) in \
        enumerate(zip([star_file, star_file, star_file,
                       star_file, agn_file, sne_file],
                      ['r', 'r', 'r', 'r', 'b', 'g'],
                      ['stars (all)', 'M dwarfs', 'RRLy', 'kplr', 'agn', 'sne'],
                      [None, 2, 3, 1, None, None])):

            linestyle='-'
            if pop == 'M dwarfs' or pop=='RRLy' or pop=='kplr':
                linestyle='--'
            date_arr = get_rate(file_name, cut=cut)

            date_arr = np.round(date_arr).astype(int)
            n_events = len(date_arr)
            i_yr = []
            avg_rate = []
            for ii in range(10):
                valid = np.where(np.logical_and(date_arr>59580+ii*365,
                                          date_arr<=59580+(ii+1)*365))
                i_yr.append(ii+1)
                valid_date = date_arr[valid]
                #print(pop,ii,len(valid_date))
                if len(valid_date)>0:
                    avg_rate.append(len(valid_date)/(len(np.unique(valid_date))))
                else:
                    avg_rate.append(0.0)

            hh, = plt.plot(i_yr, avg_rate, linestyle=linestyle, linewidth=3)
            handle_list.append(hh)
            name_list.append('%s (%.2e)' % (pop, n_events))
        plt.xticks(fontsize=20)
        plt.xlabel('year', fontsize=20)
        plt.yticks(fontsize=20)
        plt.yscale('log')
        plt.ylim((1,2000000))
        plt.ylabel('alerts/day', fontsize=20)
        plt.title(region_dict[region_name],
                      fontsize=20)
        plt.legend(handle_list, name_list, loc=0, fontsize=20)

        #plt.legend([star_handle, agn_handle, sne_handle],
        #           ['stars', 'agn', 'SNe'], fontsize=20, loc=0)

    plt.tight_layout()
    plt.savefig('plots/all_regions_avg_rate.png')
    plt.close()
