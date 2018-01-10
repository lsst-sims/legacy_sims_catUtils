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
plt.contour(counts.transpose(),
            extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()],
            colors='r')
   
counts, xbins, ybins = np.histogram2d(eb_data['x'], eb_data['y'], bins=n_bins)
plt.contour(counts.transpose(),
            extent=[xbins.min(), xbins.max(), ybins.min(), ybins.max()],
            colors='b')

plt.savefig('test_tsne.png')
