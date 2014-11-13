"""
This script reads in the galSim_example.txt catalog created by
galSimCatalogGenerator.py and uses it to draw FITS files for each
of the detectors defined in that catalog.

See the documentation at the top of galSimCatalogGenerator.py for an
explanation of why the two scripts must be separate.

Run this script using the python for which you have GalSim installed.

Be sure to run it in a shell in which all of the LSST stack environment
variables have been set, otherwise you will not have access to all of the
photUtils functionality needed by the GalSimInterpreter.
"""

import os
from lsst.sims.catUtils.exampleCatalogDefinitions import GalSimInterpreter

#specify a bandpass through which to observe the galaxies
bandPass = os.path.join(os.getenv('THROUGHPUTS_DIR'),'baseline','total_g.dat')

gs = GalSimInterpreter()

#read in our catalog of galaxy bulges
gs.readCatalog('galSim_example.txt')

#write the images to files of the name galsimTest_detectorName.fits
name = 'galsimTest_'
gs.drawCatalog(bandPass=bandPass, fileNameRoot=name)
