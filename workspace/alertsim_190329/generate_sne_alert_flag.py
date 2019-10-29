from lsst.sims.catalogs.db import DBObject
db = DBObject(database='LSST', host='localhost',
              port=51432, driver='mssql+pymssql')

import numpy as np

rng = np.random.RandomState(7153323)
dtype = np.dtype([('id', int), ('htmid', int)])
query = 'SELECT id, htmid FROM galaxy'

data_iter = db.get_arbitrary_chunk_iterator(query, dtype=dtype,
                                            chunk_size=10000)

with open('galaxy_sne_flag.txt', 'w') as out_file:
    for chunk in data_iter:
        flag_vals = rng.randint(0,10,size=len(chunk))
        for hh, ii, ff in zip(chunk['htmid'], chunk['id'], flag_vals):
            out_file.write('%d;%d;%d\n' % (hh, ii, ff))
