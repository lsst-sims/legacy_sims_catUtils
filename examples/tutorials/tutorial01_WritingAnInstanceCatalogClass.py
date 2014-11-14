"""
This tutorial shows how to write getters so that an InstanceCatalog
daughter class can calculate new column names
"""

import numpy
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import *
from lsst.sims.catalogs.measures.instance import InstanceCatalog, cached, compound

class TutorialCatalog(InstanceCatalog):
    column_outputs = ['raJ2000', 'decJ2000', 'sum', 'difference', 'quotient',
                      'df1', 'df2']

    default_columns = [('df1', 3.5, float), ('df2', 'default', (str,7))]

    @cached
    def get_quotient(self):
        ra = self.column_by_name('raJ2000')
        dec = self.column_by_name('decJ2000')
        return ra/dec
    
    @compound('sum', 'difference')
    def get_sumAndDifference(self):
        ra = self.column_by_name('raJ2000')
        dec = self.column_by_name('decJ2000')
        return numpy.array([ra+dec, ra-dec])

myDB = CatalogDBObject.from_objid('allstars')
obs_metadata = ObservationMetaData(unrefractedRA=220.0, unrefractedDec=19.0,
                                   boundType='circle', boundLength=0.1)

cat = TutorialCatalog(myDB, obs_metadata=obs_metadata)
cat.write_catalog('tutorial_on_getters.txt')
