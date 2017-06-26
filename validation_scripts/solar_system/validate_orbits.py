"""
This script will validate the MBA tables on fatboy by querying a random
sample of asteroids and using the orbit propagation and ephemeride
generation code in sims_movingObjects to verify that the position parametrized
on fatboy is within 5 milliarcseconds of the position calculated directly by
pyoorb.

Running this script will require installing and setting up sims_movingObjects.

You will also need to compile the coordinate transformation code coordtrans.c
by running
    gcc -lm -o coord coordtrans.c
"""

from subprocess import check_output

def equatorial_from_ecliptic(lon, lat):
    """
    Convert ecliptic coordinates to RA, Dec.
    Not vectorized; must do coordinate pairs one-at-a-time.

    Parameters:
    -----------
    lon is ecliptic longitude in degrees

    lat is ecliptic latitude in degrees

    Returns:
    -------
    RA in degrees

    Dec in degrees
    """

    try:
        result = check_output(["./coord", "-e", "%.6f" % lon, "%.6f" % lat])
    except OSError:
        print('\nTrying to call ./coord resulted in an exception\n'
              'Are you sure you have compiled coordtrans.c?\n'
              '    gcc -lm -o coord coordtrans.c\n')
        raise
    result = result.split('\n')
    ra_dec = result[1]
    ra_dec = ra_dec.split(':')[1]
    ra_dec = ra_dec.strip().split()
    return 15.0*float(ra_dec[0]), float(ra_dec[1])

from lsst.sims.catUtils.baseCatalogModels import MBAObj
from lsst.sims.utils import ObservationMetaData

try:
    # if you are on UW campus/VPN
    mba_db = MBAObj(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
                    port=1433, driver='mssql+pymssql')
except:
    # if you need to use the SSH tunnel to connect to fatboy
    mba_db = MBAObj()

ra, dec = equatorial_from_ecliptic(213.0, 1.1)

obs = ObservationMetaData(mjd=60121.67,
                          pointingRA=ra,
                          pointingDec=dec,
                          boundLength=1.75,
                          boundType='circle')

colnames = ['raJ2000', 'decJ2000']

import time

t_start = time.time()
results = mba_db.query_columns(colnames=colnames,
                               obs_metadata=obs,
                               chunk_size=10000)

line_ct = 0
for chunk in results:
    print chunk
    line_ct += len(chunk)
print 'got %d in %e' % (line_ct, time.time()-t_start)
