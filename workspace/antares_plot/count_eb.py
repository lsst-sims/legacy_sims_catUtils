from lsst.sims.catUtils.baseCatalogModels import StarBase
from lsst.sims.catUtils.utils import StellarAlertDBObjMixin
import json
import os

class OBAFGKObj(StellarAlertDBObjMixin, StarBase):
    objid = 'obafgkstars'
    tableid = 'StarOBAFGKForceseek'
    objectTypeId = 14
    doRunTest = True
    #These types should be matched to the database.
    #: Default map is float.  If the column mapping is the same as the column name, None can be specified
    columns = [('id','simobjid', int),
               ('raJ2000', 'ra*PI()/180.'),
               ('decJ2000', 'decl*PI()/180.'),
               ('glon', 'gal_l*PI()/180.'),
               ('glat', 'gal_b*PI()/180.'),
               ('magNorm', '(-2.5*log(flux_scale)/log(10.)) - 18.402732642'),
               ('properMotionRa', '(mura/(1000.*3600.))*PI()/180.'),
               ('properMotionDec', '(mudecl/(1000.*3600.))*PI()/180.'),
               ('parallax', 'parallax*PI()/648000000.'),
               ('galacticAv', 'CONVERT(float, ebv*3.1)'),
               ('radialVelocity', 'vrad'),
               ('variabilityParameters', 'varParamStr', str, 256),
               ('sedFilename', 'sedfilename', str, 40)]

eb_cat = os.path.join('/Users', 'danielsf', 'physics', 'lsst_160212',
                      'Development', 'CatSimLC', 'data',
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
