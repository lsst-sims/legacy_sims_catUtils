"""
Test the simulation of a SNIa light curve using snObject at specified ra, dec,
and specific SN (SALT2) parameters. 
"""
import os
import numpy as np
import matplotlib.pyplot as plt

from lsst.sims.photUtils.Photometry import PhotometryStars, Sed, Bandpass
from lsst.sims.photUtils.Photometry import PhotometryBase
from snObject import SNObject

import utils_for_test as tu
# from unittest import assertAlmostEqual
ra = 204.
dec = -30.
SN = SNObject(ra, dec)
SN.set(x0=1.847e-6, x1=0.1, c=0., z=0.2)

bandPassList = ['u', 'g', 'r', 'i', 'z', 'y']
photometry = PhotometryBase()
photometry.loadBandPassesFromFiles(bandPassList)
lsstbands = photometry.bandPassList
l = []
for time in np.arange(-20., 50., 1.0):
    t = time*np.ones(len(bandPassList))
    t.tolist()
    x = SN.lsstbandmags(lsstbands, time=time)
    # y = SNCosmoSN.bandmag(band=sncosmobands, time=t, magsys='ab')
    e = [time]
    e += x.tolist()
    # e += y.tolist()
    l.append(e)
header = "time(mjd) u g r i z y"
lc = np.array(l)
np.savetxt('testData/lc.dat', lc, fmt='%10.6f', header=header)

# The 'testData/standard_lc.dat' file was written using the following command
#     for a SN at ra  = 204., -30. , with x0 = 1.846e-6, x1=0.1, c=0., z=0.2
#     using the SALT2-extended model
# np.savetxt('testData/standard_lc.dat', np.array(l), fmt='%10.6f', header=header)
std_lc = np.loadtxt('testData/standard_lc.dat')
np.testing.assert_allclose(std_lc, lc)
#ass
