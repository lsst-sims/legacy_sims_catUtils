"""
This is a version of tutorial02_UsingMixins.py without the running commentary
"""

import numpy
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import *
from lsst.sims.catalogs.measures.instance import InstanceCatalog, cached, compound

def radiansToArcSec(value):
    """
    An example unit transformation that converts radians into arc seconds
    """
    return 3600.0*numpy.degrees(value)

class ExampleMixin(object):
    """
    An example mixin that provides getters to TutorialCatalog
    """

    def get_raInArcSec(self):
        """
        Returns RA in radians.  Will be converted to arcseconds by transformations = {}
        """
        return self.column_by_name('raJ2000')

    @cached
    def get_raPlusOneRadian(self):
        rr = self.column_by_name('raJ2000')
        return rr+1.0

    @compound('sum', 'difference')
    def get_math(self):
        rr = self.column_by_name('raJ2000')
        dd = self.column_by_name('decJ2000')
        return numpy.array([rr+dd, rr-dd])

class TutorialCatalog(InstanceCatalog, ExampleMixin):
    """
    An example InstanceCatalog that relies on ExampleMixin to provide getters for some
    of its columns
    """
    column_outputs = ['raJ2000', 'decJ2000', 'raPlusOneRadian', 'sum', 'difference',
                     'raInArcSec']

    #Recall that all angles are manipulated as radians inside the code.
    #Therefore, to get outputs in degrees, we must define transformations
    #for the columns we want to transform.
    #
    #Note that 'raPlusOneRadian' is not converted and will thus be written
    #in radians.
    transformations = {'raJ2000':numpy.degrees, 'decJ2000':numpy.degrees,
                       'sum':numpy.degrees, 'difference':numpy.degrees,
                       'raInArcSec':radiansToArcSec}

    #This is the key value that needs to be passed to CatalogDBObject.getCatalog()
    #in order to instantiate a TutorialCatalog
    catalog_type = 'tutorial_catalog'





myDB = CatalogDBObject.from_objid('allstars')
obs_metadata = ObservationMetaData(unrefractedRA=220.0, unrefractedDec=19.0,
                                   boundType='circle', boundLength=0.1,
                                   mjd=52000.0)

#First just write a catalog the way we are used to
cat = TutorialCatalog(myDB, obs_metadata=obs_metadata)
cat.write_catalog('tutorial_mixin_catalog.txt')


#Now use CatalogDBObject.getCatalog() to write a catalog (using a different
#ObservationMetaData)
obs_metadata = ObservationMetaData(unrefractedRA=120.0, unrefractedDec=-5.0,
                                   boundType='circle', boundLength=0.1,
                                   mjd=52000.0)

cat = myDB.getCatalog('tutorial_catalog', obs_metadata=obs_metadata)
cat.write_catalog('tutorial_mixin_get_catalog.txt')
