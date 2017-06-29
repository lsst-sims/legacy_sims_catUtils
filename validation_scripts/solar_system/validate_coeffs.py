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

dtype = np.dtype([('id', int), ('start', float), ('end', float)])

import time
t_start = time.time()
query = 'SELECT ssmid, mjd_start, mjd_end FROM C14_MBA_59580 '
query += 'ORDER BY mjd_start' % id_list[0]

prev_end = np.zeros(len(id_list)+1)

tol = 1.0e-6

results = db.get_arbitrary_chunk_iterator(query, dtype=dtype, chunk_size=10000)
row_ct = 0
for chunk in results:
    for line in chunk:
        row_ct += 1
        if prev_end[line['id']]<1.0:
            assert np.abs(line['start']-59580.0)<tol
        else:
            assert np.abs(line['start']-prev_end[line['id']])<tol
        prev_end[line['id']] = line['end']
        if row_ct%1000==0:
            elapsed = time.time()-t_start
            expected = 2.09e8*elapsed/row_ct
            expected = expected/(24.0*60.0*60.0)
            print row_ct,elapsed, expected

invalid = np.where(np.abs(prev_end-59610.0)>tol)
assert len(invalid[0]) == 1

print 'that took ',time.time()-t_start
