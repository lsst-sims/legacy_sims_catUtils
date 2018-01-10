import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np

def plot_contour(x,y,colors):
    counts, xbins, ybins = np.histogram2d(x, y, bins=100)

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

    one_sig_entries = np.where(flat_counts>=one_sig_val)
    one_sig_sum = flat_counts[one_sig_entries].sum()
    three_sig_entries = np.where(flat_counts>=three_sig_val)
    three_sig_sum = flat_counts[three_sig_entries].sum()

    print(tot_sum,sorted_counts[-1],two_sig_val,one_sig_val)
    print(one_sig_sum/tot_sum,three_sig_sum/tot_sum)

    plt.contour(counts.transpose(),[three_sig_val],
                extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()],
                colors=colors)

    plt.contourf(counts.transpose(), [one_sig_val, sorted_counts.max()],
                 extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()],
                 colors=colors, alpha=0.5)

if __name__ == "__main__":

    dtype = np.dtype([('x', float), ('y', float)])

    rrly_data = np.genfromtxt('rrly_features_180108_tsne_features_100k.txt', dtype=dtype)
    eb_data = np.genfromtxt('eb_features_180109_tsne_features_100k.txt', dtype=dtype)

    plt.figsize = (30,30)

    plot_contour(rrly_data['x'], rrly_data['y'], colors='r')
    plot_contour(eb_data['x'], eb_data['y'], colors='b')

    plt.savefig('test_tsne_100k.png')
