from lsst.sims.catUtils.baseCatalogModels import StarBase
from lsst.sims.catUtils.utils import StellarAlertDBObjMixin

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
