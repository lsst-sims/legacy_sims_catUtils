"""
The catalog examples in this file will write catalog files that can be read
by galSimInterpreter.py (to be written), which will use galSim to turn them
into an image.
"""

import numpy
import os
from lsst.sims.catalogs.measures.instance import InstanceCatalog, cached, is_null
from lsst.sims.coordUtils import CameraCoords, AstrometryGalaxies
from lsst.sims.photUtils import EBVmixin, Sed, CosmologyWrapper, Bandpass
import lsst.afw.cameraGeom.testUtils as camTestUtils
import lsst.afw.geom as afwGeom
from lsst.afw.cameraGeom import PUPIL, PIXELS, FOCAL_PLANE

#if you want to use the actual LSST camera:
#from lsst.obs.lsstSim import LsstSimMapper

__all__ = ["GalSimGalaxies"]

def radiansToArcsec(value):
    return 3600.0*numpy.degrees(value)

class GalSimBase(InstanceCatalog, CameraCoords):

    cannot_be_null = ['galSimSedName']

    column_outputs = ['galSimType', 'chipName', 'x_pupil', 'y_pupil', 'galSimSedName',
                      'majorAxis', 'minorAxis', 'sindex', 'halfLightRadius',
                      'positionAngle']

    transformations = {'x_pupil':radiansToArcsec,
                       'y_pupil':radiansToArcsec}

    default_formats = {'S':'%s', 'f':'%.9g', 'i':'%i'}

    delimiter = ';'

    objectCounter = 0
    sedDirName = 'galSimSedDir'

    camera = camTestUtils.CameraWrapper().camera

    #if you want to use the actual LSST camera
    #camera = LsstSimMapper().camera

    def get_sedFilepath(self):
        #copied from the phoSim catalogs
        return numpy.array([self.specFileMap[k] if self.specFileMap.has_key(k) else None
                         for k in self.column_by_name('sedFilename')])

    def get_galSimSedName(self):
        if os.path.exists(self.sedDirName):
            if not os.path.isdir(self.sedDirName):
                os.unlink(self.sedDirName)
        if not os.path.exists(self.sedDirName):
            os.mkdir(self.sedDirName)

        returnNames = []
        directory = os.path.join(os.getcwd(),self.sedDirName)
        return self.calculateGalSimSed(directory)

    def write_header(self, file_handle):
        """
        Overwrite the write_header method from InstanceCatalog because, in order to run GalSim,
        we need to print information about the detectors into the header of the catalog file
        """

        for dd in self.camera:
            cs = dd.makeCameraSys(PUPIL)
            centerPupil = self.camera.transform(dd.getCenter(FOCAL_PLANE),cs).getPoint()
            centerPixel = dd.getCenter(PIXELS).getPoint()

            translationPixel = afwGeom.Point2D(centerPixel.getX()+1, centerPixel.getY()+1)
            translationPupil = self.camera.transform(
                                    dd.makeCameraPoint(translationPixel, PIXELS), cs).getPoint()

            plateScale = numpy.sqrt(numpy.power(translationPupil.getX()-centerPupil.getX(),2)+
                                    numpy.power(translationPupil.getY()-centerPupil.getY(),2))/numpy.sqrt(2.0)
            xmax = None
            xmin = None
            ymax = None
            ymin = None
            for corner in dd.getCorners(FOCAL_PLANE):
                pt = self.camera.makeCameraPoint(corner, FOCAL_PLANE)
                pp = self.camera.transform(pt, cs).getPoint()
                if xmax is None or pp.getX() > xmax:
                    xmax = pp.getX()
                if xmin is None or pp.getX() < xmin:
                    xmin = pp.getX()
                if ymax is None or pp.getY() > ymax:
                    ymax = pp.getY()
                if ymin is None or pp.getY() < ymin:
                    ymin = pp.getY()

            xcenter = 3600.0*numpy.degrees(centerPupil.getX())
            ycenter = 3600.0*numpy.degrees(centerPupil.getY())
            xmin = 3600.0*numpy.degrees(xmin)
            xmax = 3600.0*numpy.degrees(xmax)
            ymin = 3600.0*numpy.degrees(ymin)
            ymax = 3600.0*numpy.degrees(ymax)
            plateScale = 3600.0*numpy.degrees(plateScale)

            file_handle.write('#detector;%s;%f;%f;%f;%f;%f;%f;%f\n' %
                             (dd.getName(), xcenter, ycenter, xmin, xmax, ymin, ymax, plateScale))

        InstanceCatalog.write_header(self, file_handle)

class GalSimGalaxies(GalSimBase, AstrometryGalaxies, EBVmixin):
    catalog_type = 'galsim_galaxy'
    default_columns = [('galacticAv', 0.1, float),
                       ('galSimType', 'galaxy', (str,6))]

    def calculateGalSimSed(self, outputDirectory):
        actualSEDnames = self.column_by_name('sedFilepath')
        redshift = self.column_by_name('redshift')
        internalAv = self.column_by_name('internalAv')
        internalRv = self.column_by_name('internalRv')
        galacticAv = self.column_by_name('galacticAv')
        galacticRv = self.column_by_name('galacticRv')
        magNorm = self.column_by_name('magNorm')

        sedDir = os.getenv('SIMS_SED_LIBRARY_DIR')

        #for setting magNorm
        imsimband = Bandpass()
        imsimband.imsimBandpass()

        cosmology = CosmologyWrapper()

        outputNames=[]

        for (sedName, zz, iAv, iRv, gAv, gRv, norm) in \
            zip(actualSEDnames, redshift, internalAv, internalRv, galacticAv, galacticRv, magNorm):

            if is_null(sedName):
                outputNames.append('None')
            else:
                #load the SED of the object
                sedFile = os.path.join(sedDir, sedName)
                sed = Sed()
                sed.readSED_flambda(sedFile)

                #normalize the SED
                fNorm = sed.calcFluxNorm(norm, imsimband)
                sed.multiplyFluxNorm(fNorm)

                #apply dust extinction (internal)
                a_int, b_int = sed.setupCCMab()
                sed.addCCMDust(a_int, b_int, A_v=iAv, R_v=iRv)

                #apply redshift; do not dim the SED; that should be left to the CosmologyWrapper
                sed.redshiftSED(zz, dimming=False)

                #apply cosmological distance modulus
                luminosity = sed.calcFlux(imsimband)
                distanceModulus = cosmology.distanceModulus(zz)
                flux = luminosity * numpy.power(10.0, -0.4*distanceModulus)
                sed.multiplyFluxNorm(flux)

                #apply dust extinction (galactic)
                a_int, b_int = sed.setupCCMab()
                sed.addCCMDust(a_int, b_int, A_v=gAv, R_v=gRv)

                individualName = ('galSimSed_galaxy_%d.dat' % self.objectCounter)
                name = os.path.join(outputDirectory,individualName)
                outputNames.append(name)
                self.objectCounter += 1

                sed.writeSED(name)

        return numpy.array(outputNames)
