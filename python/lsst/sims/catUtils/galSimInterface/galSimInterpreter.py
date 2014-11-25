"""
This script contains all of the methods needed to turn a GalSimCatalog
(see exampleCatalogExamples/galSimCatalogExamples.py) into FITS images

There is the class GalSimInterpreter, which can be imported into other scripts
(assuming you are running a stack-enabled Python) and a main function which uses
it to create FITS images.
"""

import os
import numpy
import galsim

__all__ = ["GalSimInterpreter"]

class GalSimDetector(object):
    """
    This class stores information about individual detectors for use by the GalSimInterpreter
    """

    def __init__(self, name=None, xCenter=None, yCenter=None,
                 xMin=None, xMax=None, yMin=None, yMax=None,
                 plateScale=None):
        """
        param [in] name is a string denoting the name of the detector (this should be the
        same name that will be returned by the astrometry method findChipName()

        param [in] xCenter is the x pupil coordinate of the center of the detector in arcseconds

        param [in] yCenter is the y pupil coordinate of the cneter of the detector in arcseconds

        param [in] xMin, xMax, yMin, yMax are the corresponding minimum and maximum values of the
        pupil coordinates on this detector in arcseconds

        param [in] plateScale in arcseconds per pixel on this detector
        """

        self.name = name
        self.xCenter = xCenter
        self.yCenter = yCenter
        self.xMin = xMin
        self.xMax = xMax
        self.yMin = yMin
        self.yMax = yMax
        self.plateScale = plateScale

class GalSimInterpreter(object):
    """
    This class will read in a GalSimCatalog from a file and write
    FITS images for each of the detectors specified in the catalog.

    The basic work flow is:

    myInterpreter = GalSimInterpreter()

    myInterpreter.readCatalog('myCatalog.dat')

    myInterpreter.drawCatalog(fileNameRoot='myFitsImages', bandpass='/path/to/bandpass/file.dat')

    This will produce a series of FITS files named myFitsImages_R_0_0_S_0_0.fits, etc., one for each
    detector in the catalog.
    """

    def __init__(self):

        #in case we want to draw images using the Fourier transform
        self.bigfft = galsim.GSParams(maximum_fft_size=10000)

        self.data = None
        self.detectors = []

        #this is a list of which detectors which objects fall on
        #(in case an object is near the edge of a detector and some
        #of the light from the object falls on that detector)
        self.chipsImpinged = None

    def readCatalog(self, catalogFile):
        """
        A method to read in a catalog file

        param [in] catalogFile is a string denoting the name of the catalog file to be read

        This method will store the objects in the catalog in self.data and the detectors
        in the camera (also described by the catalog file) in self.detectors

        For an example of the type of catalog file desired, look at the catalog classes in
        exampleCatalogDefinitions/galSimCatalogExamples.py
        """

        #these are the datatypes associated with each of the columns
        dataTypes={
                   'x_pupil':float, 'y_pupil':float, 'chipName':(str, 126),
                   'galSimSedName':(str,200), 'positionAngle': float,
                   'majorAxis':float, 'minorAxis':float,
                   'sindex':float, 'halfLightRadius':float, 'galSimType':(str,126)
                   }

        #this is the datatype used for any extraneous columns that make it into the catalog
        defaultType = (str,200)

        cat = open(catalogFile,'r')
        lines = cat.readlines()
        cat.close()

        #read in the header data from the catalog
        for line in lines:
            if line[0] != '#':
                #we have now read in all of the header data from the catalog
                break

            line = line.replace("#","").strip()
            line = numpy.array(line.split(';'))

            #read in the columns and construct a numpy.dtype from them
            if line[0] == 'galSimType':
                dtype = numpy.dtype(
                        [(ww, dataTypes[ww]) if ww in dataTypes else (ww, defaultType) for ww in line]
                        )

            #read in the data describing a detector
            if line[0] == 'detector':
                name = line[1]
                xCenter = float(line[2])
                yCenter = float(line[3])
                xMin = float(line[4])
                xMax = float(line[5])
                yMin = float(line[6])
                yMax = float(line[7])
                plateScale = float(line[8])

                detector = GalSimDetector(name=name, xCenter=xCenter, yCenter=yCenter,
                                          xMin=xMin, xMax=xMax, yMin=yMin, yMax=yMax,
                                          plateScale=plateScale)

                self.detectors.append(detector)

        #now read in the data from the catalogFile
        self.data = numpy.genfromtxt(catalogFile, dtype=dtype, delimiter=';')

        #now go through each object in self.data and check to see if it is likely to
        #illumine a detector other than the one it is assigned to
        self.chipsImpinged = []
        for entry in self.data:
            chipList = []
            chipList.append(entry['chipName'])
            for detector in self.detectors:
                if detector.name not in chipList and self._doesObjectImpingeOnChip(obj=entry, detector=detector):
                    chipList.append(detector.name)

            self.chipsImpinged.append(chipList)

    def _doesObjectImpingeOnChip(self, obj=None, detector=None):
        """
        Compare an object to a detector and determine whether or not that object will cast any
        light on that detector (in case the object is near the edge of a detector and will cast some
        incidental light onto it).

        param [in] obj the row from self.data describing the astronomical object

        param [in] detector a GalSimDetector object describing the detector

        returns True if the object does illuminate the detector; False if not
        """

        #12 November 2014
        #The method below is not terribly clever.  It triples the half light radius of the object,
        #maps it to an ellipse with the same major-to-minor axis ratio as the actual object,
        #then assumes that the object is a circle with this (enlarged) half light radius.  Any detector
        #that is within this radius of the object is considered illuminated by the object.
        #This was meant to be a very safe estimate (choosing detectors that are farther from the object
        #than they probably need to be).  Something more clever and targeted is welcome

        if obj['x_pupil'] >= detector.xMin and \
           obj['x_pupil'] <= detector.xMax:

            isBetweenX = True
        else:
            isBetweenX = False

        if obj['y_pupil'] >= detector.yMin and \
           obj['y_pupil'] <= detector.yMax:

            isBetweenY = True
        else:
            isBetweenY = False

        if isBetweenX and isBetweenY:
            #the object falls on the detector directly
            return True

        radius = 3.0*obj['halfLightRadius']
        ratio = obj['minorAxis']/obj['majorAxis']
        majorAxis = radius/numpy.sqrt(ratio)

        #check if light from the object bleed across any of the detector's boundaries
        if isBetweenY:
            if obj['x_pupil'] <= detector.xMin and detector.xMin - obj['x_pupil'] < majorAxis:
                return True

            if obj['x_pupil'] >= detector.xMax and obj['x_pupil'] - detector.xMax < majorAxis:
                return True

        if isBetweenX:
            if obj['y_pupil'] <= detector.yMin and detector.yMin - obj['y_pupil'] < majorAxis:
                return True

            if obj['y_pupil'] >= detector.yMax and obj['y_pupil'] - detector.yMax <majorAxis:
                return True

        #see if light from the object bleeds through any of the detector's corners
        for xx in [detector.xMin, detector.xMax]:
            for yy in [detector.yMin, detector.yMax]:
                distance = numpy.sqrt(numpy.power(xx - obj['x_pupil'],2) + \
                           numpy.power(yy - obj['y_pupil'],2))

                if distance < majorAxis:
                    return True

        return False

    def drawCatalog(self, fileNameRoot=None, bandPass=None):
        """
        Loop through the detectors read in by readCatalog and create FITS images for each

        param [in] fileNameRoot is a string containing the root of the names of the resulting FITS
        files.  The output files will be named fileNameRoot_detectorName.fits, i.e.
        fileNameRoot_R_0_0_S_1_2.fits, etc.

        param [in] bandPass is a string containing the name of the file which contains the filter
        bandpass over which you want to integrate the flux.  This file should be two columns: wavelength
        in nanometers and throughput from 0.0 to 1.0
        """

        #read in the bandpass
        bp = galsim.Bandpass(bandPass)
        for detector in self.detectors:
            self.drawImage(fileNameRoot=fileNameRoot, bandPass=bp, detector=detector)

    def drawImage(self, fileNameRoot=None, bandPass=None, detector=None):
        """
        Draw the FITS file associated with a specific detector

        param [in] fileNameRoot is a string containing the root of the names of the resulting fits
        files.  The output files will be named fileNameRoot_detectorName.fits, i.e.
        fileNameRoot_R_0_0_S_1_2.fits, etc.


        param [in] bandPass is galsim.bandpass object denoting the bandpass over which to integrate flux

        param [in] detector is a GalSimDetector object indicating the detector for which to draw the FITS file
        """

        #format the name of the detector to add to the name of the FITS file
        detectorName = detector.name
        detectorName = detectorName.replace(',','_')
        detectorName = detectorName.replace(':','_')
        detectorName = detectorName.replace(' ','_')

        fileName = fileNameRoot+detectorName+'.fits'

        #set the size of the image
        nx = int((detector.xMax - detector.xMin)/detector.plateScale)
        ny = int((detector.yMax - detector.yMin)/detector.plateScale)
        image = galsim.Image(nx,ny)

        #draw all of the objects who illumine this detector
        drawn = 0
        if self.data is not None:
            for (ii, entry) in enumerate(self.data):
                if detector.name in self.chipsImpinged[ii]:
                    drawn += 1
                    if entry['galSimType'] == 'galaxy':
                        self.drawGalaxy(entry=entry, image=image, detector=detector, bandPass=bandPass)
                    else:
                        print "Apologies: the GalSimInterpreter does not yet have a method to draw "
                        print entry['galSimType']
                        print "objects\n"

        if drawn>0:
            image.write(fileName)

    def drawGalaxy(self, entry=None, image=None, detector=None, bandPass=None):
        """
        Draw the image of a galaxy.

        param [in] entry is an element from self.data containing the information on an astronomical object

        param [in] image is a galsim Image object into which we will draw this object

        param [in] detector is a GalSimDetector object denoting the detectror with which the image is associated

        param [in] bandPass is a galsim.bandpass object denoting the bandpass over which to integrate flux
        """

        #create a Sersic profile
        obj = galsim.Sersic(n=entry['sindex'], half_light_radius=entry['halfLightRadius'],
                            gsparams=self.bigfft)

        #turn the Sersic profile into an ellipse
        obj = obj.shear(q=entry['minorAxis']/entry['majorAxis'], beta=entry['positionAngle']*galsim.radians)

        #by default, galsim draws objects at the center of the image;
        #this will shift the object into its correct position
        dx=entry['x_pupil']-detector.xCenter
        dy=entry['y_pupil']-detector.yCenter
        obj = obj.shift(dx, dy)

        #declare a spectrum() function for galSim to use
        spectrum = galsim.SED(entry['galSimSedName'],
                                 flux_type='flambda')

        #convolve the Sersic profile with the spectrum
        obj = obj*spectrum

        #add this object to the image, integrating over the bandPass
        #Note: by specifying method='real_space', we are asking GalSim
        #to directly integrate over the pixels.  Other options are to use
        #method='fft' and integrate using a Fourier transform (though this method can
        #be finnicky for large images) or method='phot' in which case photons are drawn
        #from the SED and shot at the chip (a la phoSim)
        image = obj.drawImage(bandpass=bandPass, scale=detector.plateScale, image=image,
                                  add_to_image=True, method='real_space')


def main():
    """
    This method reads in the galSim_example.txt catalog created by
    galSimCatalogGenerator.py and uses it to draw FITS files for each
    of the detectors defined in that catalog.

    See the documentation at the top of galSimCatalogGenerator.py for an
    explanation of why the two scripts must be separate.

    Run this script using the python for which you have GalSim installed.

    Be sure to run it in a shell in which all of the LSST stack environment
    variables have been set, otherwise you will not have access to all of the
    photUtils functionality needed by the GalSimInterpreter.
    """

    try:
        #specify a bandpass through which to observe the galaxies
        bandPass = os.path.join(os.getenv('THROUGHPUTS_DIR'),'baseline','total_g.dat')
    except AttributeError:
        print "You must set the environment variable THROUGPUTS_DIR to point"
        print "to wherever the LSST throughputs are stored on your machine\n"
        print "probably something like $LSST_HOME/yourOS/throughputs/version/"
        exit()

    gs = GalSimInterpreter()

    #read in our catalog of galaxy bulges
    gs.readCatalog('galSim_example.txt')

    #write the images to files of the name galsimTest_detectorName.fits
    name = 'galsimTest_'
    gs.drawCatalog(bandPass=bandPass, fileNameRoot=name)

if __name__ == "__main__":
    main()
