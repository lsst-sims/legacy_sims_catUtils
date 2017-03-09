"""
This script will simply generate cartoon throughputs so that we can
test whether or not the GalSim interface can be pointed to alternate
throughputs.
"""

from builtins import zip
from __future__ import with_statement
import numpy

bandpassNameList = ['x', 'y', 'z']
centerList = [300.0, 700.0, 1000.0]
widthList = [100.0, 300.0, 200.0]
wavelen = numpy.arange(100.0, 2000.0, 0.2)

lineCenterList = [400.0, 600.0, 800.0]
lineWidthList = [10.0, 250, 10.0]
lineAmpList = [0.5, 0.3, 0.6]

fakeAtmo = numpy.array([1.0-(ww-100.0)*0.0002 for ww in wavelen])

for amp, center, width in zip(lineAmpList, lineCenterList, lineWidthList):
    line = numpy.array([1.0-amp*numpy.exp(-0.5*numpy.power((ww-center)/width,2))
                            for ww in wavelen])
    fakeAtmo = fakeAtmo * line

with open('fakeAtmo.dat', 'w') as outputFile:
    for ww, ff in zip(wavelen, fakeAtmo):
        outputFile.write('%e %e\n' % (ww, ff));


fakeSky = numpy.array([1.0e-15 + 2.0e-15*numpy.exp(ww*0.0003) for ww in wavelen])
with open('fakeSky.dat', 'w') as outFile:
    for ww, ff in zip(wavelen, fakeSky):
        outFile.write('%e %e\n' % (ww, ff))

fakeM1 = numpy.array([numpy.exp(-(ww-100.0)*0.0001) for ww in wavelen])
with open('fakeM1.dat', 'w') as outputFile:
    for ww, ff in zip(wavelen, fakeM1):
        outputFile.write('%e %e\n' % (ww, ff))
    
fakeM2 = numpy.array([1.0-numpy.power(0.0004*(ww-1000.0),2) for ww in wavelen])
with open('fakeM2.dat', 'w') as outputFile:
    for ww, ff in zip(wavelen, fakeM2):
        outputFile.write('%e %e\n' % (ww, ff))

for bandpass, center, width in zip(bandpassNameList, centerList, widthList):
    throughput = numpy.array([numpy.exp(-0.5*numpy.power((ww-center)/width,2))
                              for ww in wavelen])

    filterName = 'fakeFilter_' + bandpass + '.dat'
    with open(filterName, 'w') as outFile:
        for ww, ff in zip(wavelen, throughput):
            outFile.write('%e %e\n' % (ww, ff))

    total = throughput * fakeAtmo
    total = total * fakeM1
    total = total * fakeM2
    
    totalName = 'fakeTotal_' + bandpass + '.dat'
    with open(totalName, 'w') as outFile:
        for ww, ff in zip(wavelen, total):
            outFile.write('%e %e\n' % (ww, ff))
