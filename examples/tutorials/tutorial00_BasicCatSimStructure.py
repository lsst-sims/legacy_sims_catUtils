"""
This script shows how all the pieces of CatSim fit together without much
attempt at explanation

The general point of the CatSim software is to produce a catalog of objects
visible through a given telescope at a given time from a given position.

The objects are read in from a database which is handled by the CatalogDBObject
class.

The parameters of the telescope pointing (direction, date, telescope location)
are handled by the ObservationMetaData class.

This information is combined and output as a catalog by the InstanceCatalogClass.

Below we walk through a cartoon example, instantiating each of these classes one
at a time.
"""

##############CatalogDBObject

"""
CatalogDBObject is a class that connects our python code to database tables.
It is defined in

sims_catalogs_generation/python/lsst/sims/catalogs/generation/db/dbConnection.py

Different daughter classes of this method have been written to specific tables
in the fatboy databases.  CatalogDBObject contains a class method from_objid()
which allows you to instantiate these daughter classes by referring to their
objid member variable.

The daughter classes insantiated below are defined in

sims_catUtils/python/lsst/sims/catUtils/baseCatalogModels/GalaxyModels.py
sims_catUtils/python/lsst/sims/catUtils/baseCatalogModels/StarModels.py

There is also a file in that directory that defines an interface to a table
of solar system objects and one that defines an interface to the Opsim 3.61 run.
"""

from lsst.sims.catalogs.generation.db import CatalogDBObject
from lsst.sims.catUtils.baseCatalogModels import *

myGalaxyDB = CatalogDBObject.from_objid('galaxyTiled')
myStarDB = CatalogDBObject.from_objid('allstars')

##############ObservationMetaData

"""
The ObservationMetaData class defines a particular pointing of the telescope.
The InstanceCatalog classes below will end up querying our CatalogDBObject
and using an ObservationMetaData object to constrain that query so that
the written catalog only contains objects that are visible from the telescope
at a given time, in a given direction, and within a given field of view size.

The ObservationMetaData class is defined in

sims_catalogs_generation/python/lsst/sims/catalogs/generation/db/ObservationMetaData.py

Generally, an ObservationMetaData is defined by an unrefractedRA, unrefractedDec (these
are the above-atmosphere directions of the pointing; in degrees), a boundType
(either 'circle' or 'box') defining the shape of the field of view, a boundLength
(either a float or a numpy array) defining the size of the field of view in degrees,
and an MJD defining the date of the observation (though this is optional if you do not
require astrometric calculations to be done on the catalog).

You can also specify a site which is an instantiation of the Site class defined in

sims_utils/python/lsst/sims/utils/Site.py

which characterizes the telescope's location.  This defaults to the LSST site.

There are other optional arguments, mostly related to interfacing the catalog with
PhoSim.  See the class's docstring for more detailed information.
"""

from lsst.sims.catalogs.generation.db import ObservationMetaData

obs_metadata = ObservationMetaData(unrefractedRA = 220.0,
                                   unrefractedDec = 19.0,
                                   boundType = 'circle',
                                   boundLength = 0.2,
                                   mjd = 52000.0)

##############InstanceCatalog

"""
The InstanceCatalog class is defined in

sims_catalogs_measures/python/lsst/sims/catalogs/measures/instance/InstanceCatalog.py

The InstanceCatalog (or daughter classes thereof) define what data should be output
to the catalog (i.e. do you just want ra and dec, or do you also want magnitudes integrated
over your telescope's pass bands, the name of the camera chip that actually sees the object,
etc.)  Tutorials 01 and 02 will show how the InstanceCatalog actually gets and processes
this data.  Below, we will just demonstrate the user interface.

The daughter classes of InstanceCatalog used below are defined in

sims_catUtils/python/lsst/sims/catUtils/exampleCatalogDefinitions/refCatalogExamples.py

This portion of the script will result in two output files

star_example.txt
galaxy_example.txt

being written to the current working directory (the one contains only stars; the other
contains only galaxies)
"""

from lsst.sims.catUtils.exampleCatalogDefinitions import RefCatalogGalaxyBase, \
                                                  RefCatalogStarBase

myStarCat = RefCatalogStarBase(myStarDB, obs_metadata=obs_metadata)
myStarCat.write_catalog('star_example.txt')

myGalaxyCat = RefCatalogGalaxyBase(myGalaxyDB, obs_metadata=obs_metadata)
myGalaxyCat.write_catalog('galaxy_example.txt')

##############alternate ObservationMetaData

"""
Above we used an ObservationMetaData object with a circular bound on the field of view.
Below we try a square field of view, just so you can see that the results actually
come out as advertised.

This portion of the script will write

star_example_square.txt

to the current working directory
"""
squareObsMetadata = ObservationMetaData(unrefractedRA = 220.0,
                                       unrefractedDec = 19.0,
                                       boundType = 'box',
                                       boundLength = 0.3,
                                       mjd = 52000.0)

myStarCat = RefCatalogStarBase(myStarDB, obs_metadata=squareObsMetadata)
myStarCat.write_catalog('star_example_square.txt')
