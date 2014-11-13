"""
This tutorial will show an example of incorporating actual mixins from the stack
"""

import numpy
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData, haversine
from lsst.sims.catUtils.baseCatalogModels import *
from lsst.sims.catalogs.measures.instance import InstanceCatalog, cached
from lsst.sims.coordUtils import AstrometryStars
from lsst.sims.photUtils import PhotometryStars

class TutorialCatalog(InstanceCatalog, AstrometryStars, PhotometryStars):
    column_outputs = ['raJ2000', 'decJ2000', 'lsst_u', 'raObserved', 'decObserved',
                      'shift']

    transformations = {'raJ2000':numpy.degrees, 'decJ2000':numpy.degrees,
                      'raObserved':numpy.degrees, 'decObserved':numpy.degrees,
                       'shift':numpy.degrees}

    catalog_type = 'tutorial_catalog'

    @cached
    def get_shift(self):
        r0 = self.column_by_name('raJ2000')
        d0 = self.column_by_name('decJ2000')
        r1 = self.column_by_name('raObserved')
        d1 = self.column_by_name('decObserved')
        
        return haversine(r0, d0, r1, d1)

myDB = CatalogDBObject.from_objid('allstars')
obs_metadata = ObservationMetaData(unrefractedRA=220.0, unrefractedDec=19.0,
                                   boundType='circle', boundLength=0.1, 
                                   mjd=52000.0)

cat = TutorialCatalog(myDB, obs_metadata=obs_metadata)
cat.write_catalog('tutorial_on_actual_getters.txt')

obs_metadata = ObservationMetaData(unrefractedRA=120.0, unrefractedDec=-5.0,
                                   boundType='circle', boundLength=0.1, 
                                   mjd=52000.0)

cat = myDB.getCatalog('tutorial_catalog', obs_metadata=obs_metadata)
cat.write_catalog('tutorial_on_get_catalog.txt')
