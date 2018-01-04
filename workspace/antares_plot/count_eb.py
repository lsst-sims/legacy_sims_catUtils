from EclipsingBinaryModule import OBAFGKObj
from lsst.utils import getPackageDir
import json
import os

eb_cat = os.path.join(getPackageDir('sims_catUtils'),
                      'workspace', 'antares_plot', 'data',
                      'villanova_eb_catalog.txt')

if not os.path.exists(eb_cat):
    raise RuntimeError('%s does not exist' % eb_cat)

eb_set = set()
with open(eb_cat, 'r') as in_file:
    for line in in_file:
        if line[0] == '#':
            continue
        params = line.strip().split()
        eb_set.add(int(params[0]))

db = OBAFGKObj(database='LSSTCATSIM', host='fatboy.phys.washington.edu',
               port=1433, driver='mssql+pymssql')

htmid = 2683

data_iter = db.query_columns_htmid(colnames=['simobjid', 'htmid', 'varParamStr'],
                                   chunk_size=10000,
                                   htmid=htmid,
                                   constraint='imag<24.0')

total_ct = 0
eb_ct = 0
for chunk in data_iter:
    for star in chunk:
        total_ct += 1
        param_dict = json.loads(star['varParamStr'])
        kid = param_dict['p']['lc']
        if kid in eb_set:
            eb_ct += 1
    print('eb %.2e tot %.2e' % (eb_ct, total_ct))
print('eb %.2e tot %.2e' % (eb_ct, total_ct))
