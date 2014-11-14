"""
This tutorial will introduce how we use mixins to compartmentalize InstanceCatalog functionality.

"""

import numpy
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData, haversine
from lsst.sims.catUtils.baseCatalogModels import *
from lsst.sims.catalogs.measures.instance import InstanceCatalog, cached, compound

class ExampleMixin(object):

    @cached
    def get_raPlusOne(self):
        rr = numpy.degrees(self.column_by_name('raJ2000'))
        return rr+1.0
    
    @compound('sum', 'difference')
    def get_math(self):
        rr = self.column_by_name('raJ2000')
        dd = self.column_by_name('decJ2000')
        return numpy.array([rr+dd, rr-dd])

class TutorialCatalog(InstanceCatalog, ExampleMixin):
    column_outputs = ['raJ2000', 'decJ2000', 'raPlusOne', 'sum', 'difference']

    transformations = {'raJ2000':numpy.degrees, 'decJ2000':numpy.degrees,
                       'sum':numpy.degrees, 'difference':numpy.degrees}

    catalog_type = 'tutorial_catalog'

myDB = CatalogDBObject.from_objid('allstars')
obs_metadata = ObservationMetaData(unrefractedRA=220.0, unrefractedDec=19.0,
                                   boundType='circle', boundLength=0.1, 
                                   mjd=52000.0)

cat = TutorialCatalog(myDB, obs_metadata=obs_metadata)
cat.write_catalog('tutorial_mixin_catalog.txt')

obs_metadata = ObservationMetaData(unrefractedRA=120.0, unrefractedDec=-5.0,
                                   boundType='circle', boundLength=0.1, 
                                   mjd=52000.0)

cat = myDB.getCatalog('tutorial_catalog', obs_metadata=obs_metadata)
cat.write_catalog('tutorial_mixin_get_catalog.txt')
