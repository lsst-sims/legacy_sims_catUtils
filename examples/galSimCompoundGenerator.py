import os
import galsim
from lsst.sims.catalogs.generation.db import radiansToArcsec
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import *
from lsst.sims.catUtils.galSimInterface import *

class testGalSimStars(GalSimStars):
    band_pass_names = ['u','g']
    PSF = ExampleOpticalPSF()

class testGalSimGalaxies(GalSimGalaxies):
    band_pass_names = ['u', 'g']

class testGalSimAgn(GalSimAgn):
    band_pass_names = ['u', 'g']
    PSF = ExampleOpticalPSF()

#if you want to use the actual LSST camera
#from lsst.obs.lsstSim import LsstSimMapper

#select an OpSim pointing
obsMD = OpSim3_61DBObject()
obs_metadata = obsMD.getObservationMetaData(88625744, 0.05, makeCircBounds = True)

#grab a database of galaxies (in this case, galaxy bulges)
stars = CatalogDBObject.from_objid('allstars')

#now append a bunch of objects with 2D sersic profiles to our output file
stars_galSim = testGalSimStars(stars, obs_metadata=obs_metadata)

#If you want to use the LSST camera, uncomment the line below.
#You can similarly assign any camera object you want here, as long
#as you do it immediately after instantiating GalSimGalaxies
#stars_galSim.camera = LsstSimMapper().camera

catName = 'galSim_compound_examle.txt'
stars_galSim.write_catalog(catName, chunk_size=100)

print 'done with stars'

bulges = CatalogDBObject.from_objid('galaxyBulge')
bulge_galSim = testGalSimGalaxies(bulges, obs_metadata=obs_metadata)
bulge_galSim.copyGalSimInterpreter(stars_galSim)
bulge_galSim.write_catalog(catName, write_header=False,
                            write_mode='a')

print 'done with bulges'

disks = CatalogDBObject.from_objid('galaxyDisk')
disk_galSim = testGalSimGalaxies(disks, obs_metadata=obs_metadata)
disk_galSim.copyGalSimInterpreter(bulge_galSim)
disk_galSim.write_catalog(catName, write_header=False, write_mode='a')

print 'done with disks'

agn = CatalogDBObject.from_objid('galaxyAgn')
agn_galSim = testGalSimAgn(agn, obs_metadata=obs_metadata)
agn_galSim.copyGalSimInterpreter(disk_galSim)
agn_galSim.write_catalog(catName, write_header=False, write_mode='a')

agn_galSim.write_images(nameRoot='compound')
