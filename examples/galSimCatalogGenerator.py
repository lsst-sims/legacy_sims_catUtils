"""
This file shows how to generate a GalSim input catalog named galSim_example.txt.
galSimInterpreter.py will read that file in using the GalSimInterpreter class (see
galSimInterpeter.py) and turn them into FITS images.

This somewhat convoluted API (writing the catalog to a file and then reading it
in with a separate script that uses GalSim to produce an image) is made necessary
by the fact that we cannot install GalSim using the Anaconda version of Python on
Macs

stackoverflow.com/questions/23771608/trouble-installing-galsim-on-osx-with-anaconda

To wit: if the user has a system python, SCons will try to install galsim using that
python even if the user tells SCons to use the Anaconda python.  Therefore, we cannot
use the stack-supplied python to run GalSim directly for the user.
"""

from __future__ import with_statement
import os
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData
from lsst.sims.catUtils.exampleCatalogDefinitions.galSimCatalogExamples import \
        GalSimGalaxies

from lsst.sims.catUtils.baseCatalogModels import *

#select an OpSim pointing
obsMD = CatalogDBObject.from_objid('opsim3_61')
obs_metadata = obsMD.getObservationMetaData(88625744, 0.01, makeCircBounds = True)

#grab a database of galaxies (in this case, galaxy bulges)
gals = CatalogDBObject.from_objid('galaxyBulge')

#now append a bunch of objects with 2D sersic profiles to our output file
galaxy_galSim = GalSimGalaxies(gals, obs_metadata=obs_metadata)
galaxy_galSim.write_catalog('galSim_example.txt')
