"""
This tutorial will introduce how we use mixins to compartmentalize InstanceCatalog functionality.

Not every InstanceCatalog will have the same outputs.  Even those that do have the same outputs
may require different getters (one examples that comes to mind is doing photometry on stars and
galaxies; for galaxies, one must consider flux from the bulge, the disk, and the agn, all of which
are stored separately in the database; for stars there is only one source of flux to consider).
Therefore, we try as much as possible to define commonly used getters in mixins.  Mixins are
classes that define methods explicitly for the purpose of being inherited by other classes
(though, occasionally you will find a mixin with classmethods that are meant to be called
independently).  Thus mixins should not have an __init__() (since they will rely on InstanceCatalog's
__init__()).  Indeed, in the context of InstanceCatalogs, they should only contain getters and
methods in support of getters that InstanceCatalog daughter classes may need.

Below, we show an example of a mixin and how it interacts with an example InstanceCatalog.
The class TutorialCatalog calls for the columns 'raPlusOneRadian', 'sum', and 'difference.'
These columns do not exist in the database and TutorialCatalog does not contain getters defining
them.  The class ExampleMixin does contain getters defining these columns.  By making
TutorialCatalog inherit from ExampleMixin, we pass these getters on to TutorialCatalog and
thus allow it to compute the desired columns.

Two other bits of functionality are introduced below:

1) Transformations:
Sometimes you will want to define a unit transformation so the data in your written Catalog is in
different units than the data being manipulated by the code.  The most obvious example is that
CatSim policy is to store all angles in radians while they are being manipulated.  However, RA and Dec
are stored in degrees in most databases.  Therefore, RA and Dec are converted into radians when
they are passed into the code (that will be discussed in tutorial04), manipulated as radians,
and then converted back to degrees before being written out to a catalog.

Transformations of this sort can be handled by the InstanceCatalog member variable transformations.
transformations is a dict.  The keys of this dict are the names of columns to be converted.
The values of this dict are methods to be called on those columns.  For example, if your InstanceCatalog
class has

transformations = {'raJ2000':numpy.degrees}

Then, the values in 'raJ2000' will be passed through numpy.degrees() before being written to the
catalog.  We illustrate this method by converting all of the columns in TutorialCatalog to degrees
and by additionally writing out the contrived case of 'raInArcSec,' which converts the value of
'raJ2000' into arc seconds.

2) CatalogDBObject.getCatalog()
The second piece of functionality introduced below is the method CatalogDBObject.getCatalog().
Analogous to CatalogDBObject.from_objid() (see tutorial00 and tutorial04), CatalogDBObject.getCatalog()
allows the user to take a CatalogDBObject and immediately convert it into a catalog.  This is
accomplished by

cat = CatalogDBobject.getCatalog(catalogType)

where catalogType is a string corresponding to the value of the member varialbe catalog_type in
the desired InstanceCatalog daughter class.
"""

import numpy
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData, haversine
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

    #Recall that all angled are manipulated as radians inside the code.
    #Therefore, to get outputs in degrees, we must define transformations
    #for the columns we want to transform.
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
