import os
import numpy
import galsim
from lsst.sims.photUtils import Sed, Bandpass

__all__ = ["GalSimInterpreter"]

class GalSimDetector(object):
    def __init__(self, name=None, xCenter=None, yCenter=None,
                 xMin=None, xMax=None, yMin=None, yMax=None,
                 plateScale=None):

        self.name = name
        self.xCenter = xCenter
        self.yCenter = yCenter
        self.xMin = xMin
        self.xMax = xMax
        self.yMin = yMin
        self.yMax = yMax
        self.plateScale = plateScale

class GalSimInterpreter(object):

    def __init__(self):
        self.imsimband = Bandpass()
        self.imsimband.imsimBandpass()
        self.bigfft = galsim.GSParams(maximum_fft_size=10000)
        self.sedDir = os.getenv('SIMS_SED_LIBRARY_DIR')
        self.data = None
        self.detectors = []
        self.chipsImpinged = None

    def readCatalog(self, catalogFile):
        dataNeeded = ['x_pupil', 'y_pupil', 'magNorm', 'sedFilepath',
                      'redshift', 'positionAngle',
                      'galacticAv', 'galacticRv', 'internalAv', 'internalRv',
                      'majorAxis', 'minorAxis', 'sindex', 'halfLightRadius']

        dataTypes={
                   'x_pupil':float, 'y_pupil':float, 'magNorm':float, 'chipName':(str, 126),
                   'sedFilepath':(str,126), 'redshift':float, 'positionAngle': float,
                   'galacticAv':float, 'galacticRv':float, 'internalAv':float,
                   'internalRv':float, 'majorAxis':float, 'minorAxis':float,
                   'sindex':float, 'halfLightRadius':float, 'galSimType':(str,126)
                   }

        defaultType = (str,126)

        cat = open(catalogFile,'r')
        lines = cat.readlines()
        cat.close()

        for line in lines:
            if line[0] != '#':
                break

            line = line.replace("#","").strip()
            line = numpy.array(line.split(';'))

            if line[0] == 'galSimType':
                dtype = numpy.dtype(
                        [(ww, dataTypes[ww]) if ww in dataTypes else (ww, defaultType) for ww in line]
                        )

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

        print 'n_detectors ',len(self.detectors)
        self.data = numpy.genfromtxt(catalogFile, dtype=dtype, delimiter=';')

        self.chipsImpinged = []
        for entry in self.data:
            chipList = []
            chipList.append(entry['chipName'])
            for detector in self.detectors:
                if detector.name not in chipList and self._doesObjectImpingeOnChip(obj=entry, detector=detector):
                    chipList.append(detector.name)

            self.chipsImpinged.append(chipList)

    def _doesObjectImpingeOnChip(self, obj=None, detector=None):
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
            return True

        radius = 3.0*obj['halfLightRadius']
        ratio = obj['minorAxis']/obj['majorAxis']
        majorAxis = radius/numpy.sqrt(ratio)

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

        for xx in [detector.xMin, detector.xMax]:
            for yy in [detector.yMin, detector.yMax]:
                distance = numpy.sqrt(numpy.power(xx - obj['x_pupil'],2) + \
                           numpy.power(yy - obj['y_pupil'],2))

                if distance < majorAxis:
                    return True

        return False

    def drawCatalog(self, fileNameRoot=None, bandPass=None):
        for detector in self.detectors:
            self.drawImage(fileNameRoot=fileNameRoot, bandPass=bandPass, detector=detector)

    def drawImage(self, fileNameRoot=None, bandPass=None, detector=None):

        self.bandPass = galsim.Bandpass(bandPass)

        detectorName = detector.name
        detectorName = detectorName.replace(',','_')
        detectorName = detectorName.replace(':','_')
        detectorName = detectorName.replace(' ','_')

        fileName = fileNameRoot+detectorName+'.fits'

        nx = int((detector.xMax - detector.xMin)/detector.plateScale)
        ny = int((detector.yMax - detector.yMin)/detector.plateScale)
        image = galsim.Image(nx,ny)

        drawn = 0
        if self.data is not None:
            for (ii, entry) in enumerate(self.data):
                if detector.name in self.chipsImpinged[ii]:
                    drawn += 1
                    if entry['galSimType'] == 'galaxy':
                        self.drawGalaxy(entry=entry, image=image, detector=detector)
                    else:
                        print entry['galSimType']

        if drawn>0:
            image.write(fileName)

    def drawGalaxy(self, entry=None, image=None, detector=None):

        sedFile = os.path.join(self.sedDir,entry['sedFilepath'])

        obj = galsim.Sersic(n=entry['sindex'], half_light_radius=entry['halfLightRadius'],
                            gsparams=self.bigfft)

        obj = obj.shear(q=entry['minorAxis']/entry['majorAxis'], beta=entry['positionAngle']*galsim.radians)

        dx=entry['x_pupil']-detector.xCenter
        dy=entry['y_pupil']-detector.yCenter

        obj = obj.shift(dx, dy)

        sed = Sed()
        sed.readSED_flambda(sedFile)
        fNorm = sed.calcFluxNorm(entry['magNorm'], self.imsimband)
        sed.multiplyFluxNorm(fNorm)
        a_int, b_int = sed.setupCCMab()
        sed.addCCMDust(a_int, b_int, A_v=entry['internalAv'], R_v=entry['internalRv'])
        sed.redshiftSED(entry['redshift'], dimming=False)
        sed.addCCMDust(a_int, b_int, A_v=entry['galacticAv'], R_v=entry['galacticRv'])

        spectrum = galsim.SED(spec = lambda ll:
                                 numpy.interp(ll, sed.wavelen, sed.flambda),
                                 flux_type='flambda')

        obj = obj*spectrum

        image = obj.drawImage(bandpass=self.bandPass, scale=detector.plateScale, image=image,
                                  add_to_image=True, method='real_space')

