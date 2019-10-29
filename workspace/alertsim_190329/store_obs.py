"""
Get the bare minimum parameters from the OpSim database.
Store them in an hdf5 file
"""
import os
import sqlite3
import h5py
import numpy as np

import time

opsim_db = '/Users/danielsf/physics/lsst_150412'
opsim_db = os.path.join(opsim_db, 'Development/garage/OpSimData')
opsim_db = os.path.join(opsim_db,'minion_1016_sqlite.db')
assert os.path.isfile(opsim_db)

t_start = time.time()
query = "SELECT ditheredRA, ditheredDec, expMJD, obsHistID, rotSkyPos, "
query += "filter, fiveSigmaDepth FROM Summary "
query += "GROUP BY obsHistID ORDER BY obsHistID"

filt_to_int = {}
for ii, bp in enumerate('ugrizy'):
    filt_to_int[bp] = ii

conn = sqlite3.connect(opsim_db)
c = conn.cursor()
results = c.execute(query).fetchall()

ra = []
dec = []
mjd = []
obs_id = []
rotsky = []
bp = []
m5 = []

t_0 = time.time()
for rr in results:
    ra.append(float(rr[0]))
    dec.append(float(rr[1]))
    mjd.append(float(rr[2]))
    obs_id.append(int(rr[3]))
    rotsky.append(float(rr[4]))
    bp.append(filt_to_int[rr[5]])
    m5.append(float(rr[6]))

with h5py.File('data/obs_params.h5', 'w') as out_file:
    out_file.create_dataset('ra', data=np.degrees(np.array(ra)))
    out_file.create_dataset('dec', data=np.degrees(np.array(dec)))
    out_file.create_dataset('mjd', data=np.array(mjd))
    out_file.create_dataset('obsHistID', data=np.array(obs_id))
    out_file.create_dataset('rotSkyPos', data=np.degrees(np.array(rotsky)))
    out_file.create_dataset('filter', data=np.array(bp))
    out_file.create_dataset('m5', data = np.array(m5))

print('that took %e hrs' % ((time.time()-t_start)/3600.0))
