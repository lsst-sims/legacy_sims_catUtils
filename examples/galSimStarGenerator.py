import os
import galsim
from lsst.sims.catalogs.generation.db import radiansToArcsec
from lsst.sims.catalogs.measures.instance import InstanceCatalog
from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import *
from lsst.sims.catUtils.galSimInterface import *

class ExampleOpticalPSF(PSFbase):

    def _getPSF(self, x_pupil=None, y_pupil=None):
        #psf = galsim.OpticalPSF(lam_over_diam=radiansToArcsec(7.5e-8))
        psf = galsim.Gaussian(sigma=0.1)
        return psf

class testGalSimStars(GalSimStars):
    band_pass_names = ['u','g']
    
    PSF = ExampleOpticalPSF()

#if you want to use the actual LSST camera
#from lsst.obs.lsstSim import LsstSimMapper

#select an OpSim pointing
obsMD = OpSim3_61DBObject()
obs_metadata = obsMD.getObservationMetaData(88625744, 0.1, makeCircBounds = True)

#grab a database of galaxies (in this case, galaxy bulges)
stars = CatalogDBObject.from_objid('allstars')

#now append a bunch of objects with 2D sersic profiles to our output file
stars_galSim = testGalSimStars(stars, obs_metadata=obs_metadata)

#If you want to use the LSST camera, uncomment the line below.
#You can similarly assign any camera object you want here, as long
#as you do it immediately after instantiating GalSimGalaxies
#stars_galSim.camera = LsstSimMapper().camera

print obs_metadata.unrefractedRA, obs_metadata.unrefractedDec
stars_galSim.write_catalog('galSim_star_example.txt', chunk_size=100)
stars_galSim.write_images(nameRoot='star')
