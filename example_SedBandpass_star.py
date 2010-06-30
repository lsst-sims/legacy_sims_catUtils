import lsst.sims.catalogs.measures.photometry.Sed as Sed
import lsst.sims.catalogs.measures.photometry.Bandpass as Bandpass
import numpy as n

# let's treat this star as if it were many stars, with many fluxnorms and calculate the magnitude
# in a single bandpass from those many stars

# pretend we got a list (or array of strings) of star SED names
starskeys = []
for i in range(0, 10):
    starskeys.append('exampleSED.dat')

# pretend we got fluxnorm from somewhere
fluxnorm_min = 3.81e-12
fluxnorm_max = 1.2e-13
fluxnormstep = (fluxnorm_max - fluxnorm_min) / float(len(starskeys))
fluxnorm = n.arange(fluxnorm_min, fluxnorm_max, fluxnormstep)

# and similarly, that we got an array of EBV values
ebvmin = 0
ebvmax = 2
ebvstep = (ebvmax - ebvmin)/float(len(starskeys))
ebv = n.arange(ebvmin, ebvmax, ebvstep)

# okay, instantiate the bandpass used to calculate the magnitude 
rband = Bandpass()
rband.readThroughput("exampleBandpass.dat")

# instantiate stars into dictionary
stars = {}
rootdir = "./"
for key in starskeys:
    stars[key] = Sed()
    stars[key].readSED_flambda(rootdir + key)

# we know the stars have the same wavelength array (stars[i].wavelen)
# (because they were the same SED!)  BUT, in general we may not know this is true.
# so check .. (this is not strictly necessary, but does speed up some of the calculations
#  below. so, for many stars, do it.)

wavelen_base = stars[starskeys[0]].wavelen
for key in starskeys:
    if stars[key].needResample(wavelen_match = wavelen_base):
        stars[key].resampleSED(wavelen_match=wavelen_base)

# create the a,b arrays for all the stars (because we resampled the stars onto the
#  same wavelength range we can just calculate a/b once, and this is slow)

a, b = stars[starskeys[0]].setupCCMab()

# pretend we want to read mags into an array .. you could just as easily put it into a
# dictionary or list, with small variations in the code
mags = n.empty(len(starskeys), dtype='float') 
for i in range(len(starskeys)):
    # make a copy of the original SED *if* you want to 'reuse' the SED for multiple magnitude
    # calculations with various fluxnorms and dust applications (otherwise just use the object
    # you instantiated above)
    tmpstar = Sed(wavelen=stars[starskeys[i]].wavelen, flambda=stars[starskeys[i]].flambda)
    tmpstar.multiplyFluxNorm(fluxnorm[i])
    tmpstar.addCCMDust(a, b, ebv=ebv[i])
    mags[i] = tmpstar.calcMag(rband)


# show results
print "#sedname      fluxnorm     ebv    magnitude"
for i in range(len(starskeys)):
    print "%s %.5g %.5f %.5f" %(starskeys[i], fluxnorm[i], ebv[i], mags[i])
    
