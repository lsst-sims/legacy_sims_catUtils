#If you do not have GalSim installed for the version of python you use to 
#run the stack, you need to stop here and copy the code below into 
#a different script and run it using the version of python for which
#you do have GalSim installed.  Be sure to do this inside a shell in which
#the LSST environment variables have been set.  That will allow you to still
#import the photUtils functionality needed to make the GalSimInterpreter work

exit()
#specify a bandpass through which to observe the galaxies
bandPass = os.path.join(os.getenv('THROUGHPUTS_DIR'),'baseline','total_g.dat')

gs = GalSimInterpreter()

#read in our catalog of galaxy bulges
gs.readCatalog('galsim_example.txt')

#write the images to files of the name galsimTest_detectorName.fits
name = 'galsimTest_'
gs.drawCatalog(bandPass=bandPass, fileNameRoot=name)
