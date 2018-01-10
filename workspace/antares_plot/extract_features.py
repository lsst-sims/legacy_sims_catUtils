import argparse
import os
import scipy
import sqlite3
import numpy as np
import time

import FeatureExtraction

def get_lc(cursor, band, unq):
    lc_cmd = 'SELECT m.TAI, q.flux, a.dflux, a.SNR '
    lc_cmd += 'FROM alert_data AS a '
    lc_cmd += 'INNER JOIN metadata AS m ON a.obshistId=m.obshistId '
    lc_cmd += 'INNER JOIN quiescent_flux AS q ON q.uniqueId=a.uniqueId '
    lc_cmd +=  'AND q.band=m.band '
    lc_cmd += 'WHERE a.uniqueId=%d AND m.band=%d' % (unq, band)
    lc_data = cursor.execute(lc_cmd).fetchall()
    time_arr = np.zeros(len(lc_data), dtype=float)
    flux_arr = np.zeros(len(lc_data), dtype=float)
    sig_arr = np.zeros(len(lc_data), dtype=float)
    for i_f in range(len(lc_data)):
        time_arr[i_f] = lc_data[i_f][0]
        flux_arr[i_f] = lc_data[i_f][1] + lc_data[i_f][2]
        sig_arr[i_f] = flux_arr[i_f]/lc_data[i_f][3]
    sorted_dex = np.argsort(time_arr)
    time_arr = time_arr[sorted_dex]
    flux_arr = flux_arr[sorted_dex]
    sig_arr = sig_arr[sorted_dex]

    return time_arr, flux_arr, sig_arr


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--data_dir', type=str, default=None)
    parser.add_argument('--prefix', type=str, default=None)
    parser.add_argument('--out_file', type=str, default=None)
    parser.add_argument('--limit', type=int, default=None)

    args = parser.parse_args()
    if args.data_dir is None:
        raise RuntimeError("Must specify data_dir")
    if args.prefix is None:
        raise RuntimeError("Must specify prefix")
    if args.out_file is None:
        raise RuntimeError("Must specify out_file")

    list_of_files = os.listdir(args.data_dir)

    obj_ct = 0
    t_start = time.time()
    with open(args.out_file, 'w') as out_file:
        for file_name in list_of_files:
            if args.limit is not None and obj_ct>=args.limit:
                break
            if not file_name.endswith('sqlite.db'):
                continue
            if not file_name.startswith(args.prefix):
                continue
            full_name = os.path.join(args.data_dir, file_name)
            assert os.path.exists(full_name)
            elapsed = (time.time()-t_start)/3600.0
            print('connecting to %s -- obj %d elapsed %.2e hrs' % (file_name, obj_ct, elapsed))
            with sqlite3.connect(full_name) as connection:
                cursor = connection.cursor()
                unique_id_cmd = 'SELECT uniqueId FROM quiescent_flux WHERE band=1'
                unq_list = cursor.execute(unique_id_cmd).fetchall()
                for unq_val in unq_list:
                    unq = unq_val[0]
                    g_time, g_flux, g_sig = get_lc(cursor, 1, unq)
                    i_time, i_flux, i_sig = get_lc(cursor, 3, unq)
                    if g_flux.max()-g_flux.min()<1.0e-30:
                        continue
                    if i_flux.max()-i_flux.min()<1.0e-30:
                        continue
                    if len(g_flux)<3 or len(i_flux)<3:
                        continue
                    g_entropy = FeatureExtraction.entropy(g_flux, g_sig)
                    i_entropy = FeatureExtraction.entropy(i_flux, i_sig)
                    g_hlr = FeatureExtraction.hlratio(g_flux)
                    i_hlr = FeatureExtraction.hlratio(i_flux)
                    g_quart = FeatureExtraction.quartile_range(g_flux)
                    i_quart = FeatureExtraction.quartile_range(i_flux)
                    (g_skew, g_kurt,
                     g_stdevmean) = FeatureExtraction.skewness_and_kurtosis(g_flux)
                    (i_skew, i_kurt,
                     i_stdevmean) = FeatureExtraction.skewness_and_kurtosis(i_flux)
                    g_mad = FeatureExtraction.median_absolute_deviation(g_flux)
                    i_mad = FeatureExtraction.median_absolute_deviation(i_flux)
                    g_k = FeatureExtraction.stetson_k(g_flux, g_sig)
                    i_k = FeatureExtraction.stetson_k(i_flux, i_sig)
                    g_eta = FeatureExtraction.von_neumann_ratio(g_flux)
                    i_eta = FeatureExtraction.von_neumann_ratio(i_flux)
                    g_w, g_p_val = scipy.stats.shapiro(g_flux)
                    i_w, i_p_val = scipy.stats.shapiro(i_flux)

                    try:
                        (g_period_sigma,
                         g_period,
                         g_period_snr,
                         g_fap) = FeatureExtraction.periodic_features(g_time,
                                                                      g_flux,
                                                                      g_sig)
                        (i_period_sigma,
                         i_period,
                         i_period_snr,
                         i_fap) = FeatureExtraction.periodic_features(i_time,
                                                                      i_flux,
                                                                      i_sig)
                    except np.linalg.linalg.LinAlgError:
                        continue

                    if g_period_sigma < i_period_sigma:
                        period = g_period
                        period_sigma = g_period_sigma
                        period_snr = g_period_snr
                        fap = g_fap
                    else:
                        period = i_period
                        period_sigma = i_period_sigma
                        period_snr = i_period_snr
                        fap = i_fap

                    feature_vec = np.array([g_k, g_eta, g_w, g_hlr, g_kurt, g_entropy,
                                            g_mad, g_stdevmean, g_quart, g_skew,
                                            i_k, i_eta, i_w, i_hlr, i_kurt, i_entropy,
                                            i_mad, i_stdevmean, i_quart, i_skew,
                                            period, period_sigma, period_snr, fap])

                    if len(np.where(np.isnan(feature_vec))[0]) != 0:
                        continue
                    if len(np.where(np.isinf(feature_vec))[0]) != 0:
                        continue

                    out_file.write('%e %e %e %e %e %e %e %e %e %e ' %
                                   (g_k, g_eta, g_w, g_hlr, g_kurt, g_entropy,
                                    g_mad, g_stdevmean, g_quart, g_skew))
                    out_file.write('%e %e %e %e %e %e %e %e %e %e '%
                                   (i_k, i_eta, i_w, i_hlr, i_kurt, i_entropy,
                                    i_mad, i_stdevmean, i_quart, i_skew))
                    out_file.write('%e %e %e %e %d %d\n'%
                                   (period, period_sigma, period_snr, fap, len(g_time), len(i_time)))

                    obj_ct += 1
                    if args.limit is not None and obj_ct>=args.limit:
                        break
