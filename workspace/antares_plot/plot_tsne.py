import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

n_bins = 100

dtype = np.dtype([('x', float), ('y', float)])

rrly_data = np.genfromtxt('rrly_features_180108_tsne_features.txt', dtype=dtype)
eb_data = np.genfromtxt('eb_features_180109_tsne_features.txt', dtype=dtype)

plt.figsize = (30,30)

counts, xbins, ybins = np.histogram2d(rrly_data['x'], rrly_data['y'], bins=n_bins)

flat_counts = counts.flatten()
valid = np.where(flat_counts>0.0)
flat_counts = flat_counts[valid]
sorted_dex = np.argsort(-1.0*flat_counts)
sorted_counts = flat_counts[sorted_dex]
count_sum = np.cumsum(sorted_counts)
tot_sum = count_sum[-1]


one_sig_dex = np.argmin(np.abs(count_sum-0.68*tot_sum))
two_sig_dex = np.argmin(np.abs(count_sum-0.95*tot_sum))
three_sig_dex = np.argmin(np.abs(count_sum-0.99*tot_sum))

one_sig_val = sorted_counts[one_sig_dex]
two_sig_val = sorted_counts[two_sig_dex]
three_sig_val = sorted_counts[three_sig_dex]

print(tot_sum,sorted_counts[-1],two_sig_val,one_sig_val)

plt.contour(counts.transpose(),[two_sig_val, one_sig_val],
            extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()],
            colors='r')


counts, xbins, ybins = np.histogram2d(eb_data['x'], eb_data['y'], bins=n_bins)

flat_counts = counts.flatten()
valid = np.where(flat_counts>0.0)
flat_counts = flat_counts[valid]
sorted_dex = np.argsort(-1.0*flat_counts)
sorted_counts = flat_counts[sorted_dex]
count_sum = np.cumsum(sorted_counts)
tot_sum = count_sum[-1]


one_sig_dex = np.argmin(np.abs(count_sum-0.68*tot_sum))
two_sig_dex = np.argmin(np.abs(count_sum-0.95*tot_sum))
three_sig_dex = np.argmin(np.abs(count_sum-0.99*tot_sum))

one_sig_val = sorted_counts[one_sig_dex]
two_sig_val = sorted_counts[two_sig_dex]
three_sig_val = sorted_counts[three_sig_dex]

plt.contour(counts.transpose(),[two_sig_val, one_sig_val],
            extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()],
            colors='b')

plt.savefig('test_tsne.png')
