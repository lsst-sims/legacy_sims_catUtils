"""
This script will create cartoon SEDs so that we can test our ability
to re-direct the GalSim catalogs to alternate SEDs.
"""

from __future__ import with_statement
import os
import numpy

sedNameList = ['fakeSed1.dat', 'fakeSed2.dat', 'fakeSed3.dat']
wavelen = numpy.arange(100.0, 2000.0, 0.2)

functionList = [
               lambda ww: ww*0.1+100.0*numpy.exp(-0.5*numpy.power((ww-500.0)/100.0, 2)),
               lambda ww: ww*0.3+(ww*0.1)*numpy.sin(ww/400.0),
               lambda ww: 1.0-numpy.power(0.0005*(ww-500.0),2)
               ]

lineCenterList = [
                 [300.0, 500.0, 700.0],
                 [400.0, 600.0, 800.0],
                 [150.0, 800.0, 1500.0]
                 ]

lineWidthList = [
                [20.0, 30.0, 10.0],
                [50.0, 40.0, 5.0],
                [10.0, 20.0, 30.0]
                ]

lineAmpList = [
              [0.5, 0.4, 0.3],
              [0.6, 0.1, 0.2],
              [0.7, 0.5, 0.7]
              ]


for sedName, centerList, widthList, ampList, fn in \
             zip(sedNameList, lineCenterList, lineWidthList, lineAmpList, functionList):

    flux = numpy.array([fn(ww) for ww in wavelen])
    for center, width, amp in zip(centerList, widthList, ampList):
        flux *= numpy.array([1.0-amp*numpy.exp(-0.5*numpy.power((ww-center)/width,2))
                             for ww in wavelen])

    with open(sedName, 'w') as outFile:
        for ww, ff in zip(wavelen, flux):
            outFile.write('%e %e\n' % (ww, ff))
