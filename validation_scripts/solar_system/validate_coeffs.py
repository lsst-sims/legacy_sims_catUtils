"""
This script will loop through the tables of MBA orbit coefficients, verifying
that every object in our model has a set of coefficients assigned for every
time period simulated.
"""

from lsst.sims.catalogs.db import DBObject

try:
    # if you are on the University of Washington campus or VPN
    db = DBObject(database='LSSTSSM', host='fatboy.phys.washington.edu',
                  port=1433, driver='mssql+pymssql')
except:
    # try connecting through an SSH tunnel
    db = DBObject(database='LSSTSSM', host='localhost', port=51433,
                  driver='mssql_pymssql')

table_name_list = db.get_table_names()

import numpy as np

# first, get a list of all of the object IDs in the MBA model
dtype = np.dtype([('id', int)])
results = db.execute_arbitrary('SELECT ssmid from C14_MBA_name_map',
                               dtype=dtype)
id_list = results['id']

print id_list
print len(id_list)

dtype = np.dtype([('id', int), ('mjd_start', float), ('mjd_end', float)])

import time
t_start = time.time()
query = 'SELECT ssmid, mjd_start, mjd_end FROM C14_MBA_59580 '
query += 'ORDER BY mjd_start' % id_list[0]

prev_end = np.ones(len(id_list)+1, dtype=float)*59580.0

tol = 1.0e-9

results = db.get_arbitrary_chunk_iterator(query, dtype=dtype, chunk_size=100000)
row_ct = 0
overhead = time.time()-t_start
t_start = time.time()
for chunk in results:
    n_chunk = len(chunk)
    sub_divisions = 0
    is_okay = False
    while not is_okay:
        sub_divisions += 1
        is_okay = True
        for i_sub in range(sub_divisions):
            i_start = n_chunk*i_sub/sub_divisions
            i_end = n_chunk*(i_sub+1)/sub_divisions
            unq = np.unique(chunk['id'][i_start:i_end])
            if len(unq) != i_end-i_start:
                is_okay = False
                break

    for i_sub in range(sub_divisions):
        i_start = n_chunk*i_sub/sub_divisions
        i_end = n_chunk*(i_sub+1)/sub_divisions
        id_vec = chunk['id'][i_start:i_end]
        mjd_start_vec = chunk['mjd_start'][i_start:i_end]
        mjd_end_vec = chunk['mjd_end'][i_start:i_end]

        invalid = np.where(np.abs(mjd_start_vec-prev_end[id_vec])>tol)
        assert len(invalid[0]) == 0
        prev_end[id_vec] = mjd_end_vec

    row_ct += n_chunk
    elapsed = time.time()-t_start
    expected = overhead + 2.09e8*elapsed/row_ct
    expected = expected/(24.0*60.0*60.0)
    print row_ct,elapsed+overhead, expected, sub_divisions

invalid = np.where(np.abs(prev_end-59610.0)>tol)
assert len(invalid[0]) == 1

print 'that took ',time.time()-t_start
