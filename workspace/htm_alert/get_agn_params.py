from lsst.sims.catalogs.db import DBObject

db = DBObject(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
              port=1433, driver='mssql+pymssql')

query = 'SELECT galid, varParamStr '
query += 'FROM galaxy WHERE varParamStr IS NOT NULL'

import numpy as np
dtype = np.dtype([('galid', int), ('varParamStr', str, 300)])

data_iter = db.get_arbitrary_chunk_iterator(query,
                                            dtype=dtype, chunk_size=10000)

with open('agn_var_param_str.txt', 'w') as out_file:
    for chunk in data_iter:
        for gal in chunk:
            out_file.write('%d;%s\n' %
            (gal['galid'],gal['varParamStr']))
