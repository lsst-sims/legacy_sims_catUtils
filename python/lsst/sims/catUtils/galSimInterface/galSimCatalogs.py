"""
The catalog examples in this file will write catalog files that can be read
by galSimInterpreter.py (to be written), which will use galSim to turn them
into an image.

Because of the bug documented here

stackoverflow.com/questions/23771608/trouble-installing-galsim-on-osx-with-anaconda

this code must do all of the processing of the object SEDs.  We cannot assume
that the python which will convert the SEDs into images will have access to the
stack.  Therefore, the catalogs will write the processed (redshifted, dust extincted, etc.)
SEDs to a scratch directory which the image-creating script will read from.

See get_galSimSedName() in the class GalSimBase and calculateGalSimSed in its daughter
classes to see how this is implemented.
"""

import numpy
import os
import eups
from lsst.sims.catalogs.generation.db import radiansToArcsec
from lsst.sims.catalogs.measures.instance import InstanceCatalog, cached, is_null
from lsst.sims.coordUtils import CameraCoords, AstrometryGalaxies
from lsst.sims.catUtils.galSimInterface import GalSimInterpreter, GalSimDetector
from lsst.sims.photUtils import EBVmixin, Sed, Bandpass
import lsst.afw.cameraGeom.testUtils as camTestUtils
import lsst.afw.geom as afwGeom
from lsst.afw.cameraGeom import PUPIL, PIXELS, FOCAL_PLANE

__all__ = ["GalSimGalaxies"]

class GalSimBase(InstanceCatalog, CameraCoords):
    """
    This is the base class from which all GalSim catalog classes will inherit.
    It sets the column outputs and getters.

    Daughter classes of this class must set self.calculateGalSimSed() which accepts
    the name of a directory.  This method will loop through the objects in the catalog,
    calculate their SEDs (corrected for redshift, dust extinction, magNorm, etc.)
    and write them to the specified directory.  calculateGalSimSed() will then return
    a numpy array containing the full names (absolute path) of the SEDs written for
    each object in the catalog.  The catalog will contain these filenames so that
    GalSimInterpreter (defined in examples/galSimInterpreter.py) can find the SEDs when
    it is writing FITS images from the catalog.

    This class also assigns a placeholder camera to the catalog.  If the user wishes
    to use a different camera, they can assign it directly at run time after instantiating
    the GalSimCatalog class (see examples/galSimCatalogGenerator.py for a demonstration).
    """

    cannot_be_null = ['sedFilepath']

    column_outputs = ['galSimType', 'uniqueId', 'chipName', 'x_pupil', 'y_pupil', 'sedFilepath',
                      'majorAxis', 'minorAxis', 'sindex', 'halfLightRadius',
                      'positionAngle','fitsFiles']

    transformations = {'x_pupil':radiansToArcsec,
                       'y_pupil':radiansToArcsec,
                       'halfLightRadius':radiansToArcsec}

    default_formats = {'S':'%s', 'f':'%.9g', 'i':'%i'}

    delimiter = ';'

    band_pass_names = ['u','g','r','i','z','y']
    band_pass_directory = eups.productDir('throughputs')
    band_pass_root = 'baseline/total'

    #Consulting the file sed.py in GalSim/galsim/ it appears that GalSim expects
    #its SEDs to ultimately be in units of ergs/nm so that, when called, they can
    #be converted to photons/nm (see the function __call__() and the assignment of
    #self._rest_photons in the __init__() of galsim's sed.py file).  Thus, we need
    #to read in our SEDs, normalize them, and then multiply by the exposure time
    #and the effective area to get from ergs/s/cm^2/nm to ergs/nm.  
    #
    #The gain parameter should convert between photons and ADU (so: it is the
    #traditional definition of "gain" -- electrons per ADU -- multiplied by the
    #quantum efficiency of the detector).  Because we fold the quantum efficiency
    #of the detector into our total_[u,g,r,i,z,y].dat bandpass files
    #(see the readme in the THROUGHPUTS_DIR/baseline/), we only need to multiply
    #by the electrons per ADU gain.
    #
    #These parameters can be set to different values by redefining them in daughter 
    #classes of this class.
    #
    effective_area = numpy.pi*(6.5*100.0/2.0)**2 #copied from Sed.py in sims_photUtils
    exposure_time = 15.0 #copied from Sed.py in sims_photUtils
    gain = 2.3 #copied from Sed.py in sims_photUtils

    #This is just a place holder.  If you want to assign a different camera,
    #you can do so immediately after instantiating this class
    camera = camTestUtils.CameraWrapper().camera

    detectorImages = {}
    objectHasBeenDrawn = []
    uniqueSeds = {}
    catalogHasBeenWritten = False

    def get_sedFilepath(self):
        #copied from the phoSim catalogs
        return numpy.array([self.specFileMap[k] if self.specFileMap.has_key(k) else None
                         for k in self.column_by_name('sedFilename')])

    def _calculateGalSimSeds(self):
        """
        Apply any physical corrections to the objects' SEDS (redshift them, apply dust, etc.).
        Return a list of Sed objects containing the SEDS
        """
        
        sedList = []
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

        outputNames=[]

        for (sedName, zz, iAv, iRv, gAv, gRv, norm) in \
            zip(actualSEDnames, redshift, internalAv, internalRv, galacticAv, galacticRv, magNorm):

            if is_null(sedName):
                sedList.append(None)
            else:
                if sedName in self.uniqueSeds:
                    sed = Sed(wavelen=self.uniqueSeds[sedName].wavelen,
                              flambda=self.uniqueSeds[sedName].flambda,
                              fnu=self.uniqueSeds[sedName].fnu,
                              name=self.uniqueSeds[sedName].name)
                else:
                    #load the SED of the object
                    sed = Sed()
                    sedFile = os.path.join(sedDir, sedName)
                    sed.readSED_flambda(sedFile)
                    sedCopy = Sed(wavelen=sed.wavelen, flambda=sed.flambda,
                                  fnu=sed.fnu, name=sed.name)
                    self.uniqueSeds[sedName] = sedCopy

                #normalize the SED
                fNorm = sed.calcFluxNorm(norm, imsimband)
                sed.multiplyFluxNorm(fNorm*self.exposure_time*self.effective_area)

                #apply dust extinction (internal)
                if iAv != 0.0 and iRv != 0.0:
                    a_int, b_int = sed.setupCCMab()
                    sed.addCCMDust(a_int, b_int, A_v=iAv, R_v=iRv)

                #13 November 2014
                #apply redshift; there is no need to apply the distance modulus from
                #sims/photUtils/CosmologyWrapper; I believe magNorm takes that into account
                #
                #also: no need to apply dimming for the same reason
                if zz != 0.0:
                    sed.redshiftSED(zz, dimming=False)

                #apply dust extinction (galactic)
                a_int, b_int = sed.setupCCMab()
                sed.addCCMDust(a_int, b_int, A_v=gAv, R_v=gRv)
                sedList.append(sed)

        return sedList

    @cached
    def get_fitsFiles(self):
        objectNames = self.column_by_name('uniqueId')
        xPupil = self.column_by_name('x_pupil')
        yPupil = self.column_by_name('y_pupil')
        halfLight = self.column_by_name('halfLightRadius')
        minorAxis = self.column_by_name('minorAxis')
        majorAxis = self.column_by_name('majorAxis')
        positionAngle = self.column_by_name('positionAngle')
        sindex = self.column_by_name('sindex')
        sedList = self._calculateGalSimSeds()

        output = []
        for (name, xp, yp, hlr, minor, major, pa, ss, sn) in \
            zip(objectNames, xPupil, yPupil, halfLight, minorAxis, majorAxis, positionAngle,
            sedList, sindex):
            
            if name in self.objectHasBeenDrawn:
                raise RuntimeError('Trying to draw %s more than once' % str(name))
            self.objectHasBeenDrawn.append(name)
            chipsString, chipsList = self.galSimInterpreter.findAllChips(xPupil=xp, yPupil=yp,
                                                                         minorAxis=minor, majorAxis=major,
                                                                         halfLightRadius=hlr)
            output.append(chipsString)
            self.galSimInterpreter.drawObject(galSimType=self.galsim_type, detectorList=chipsList,
                                              sindex=sn, minorAxis=minor,
                                              majorAxis=major, positionAngle=pa, halfLightRadius=hlr,
                                              x_pupil=xp, y_pupil=yp, sed=ss)
        
        if hasattr(self,'galSimInterpreter'):
            self.galSimInterpreter.compressObjectList()
        return numpy.array(output)

    def _getBandPasses(self):
        bandPassFiles = []
        for bpn in self.band_pass_names:
            name = self.band_pass_directory+'/'+self.band_pass_root+'_'+bpn+'.dat'
            bandPassFiles.append(name)

        return bandPassFiles, self.band_pass_names

    def write_header(self, file_handle):
        """
        Overwrite the write_header method from InstanceCatalog because, in order to run GalSim,
        we need to print information about the detectors into the header of the catalog file
        """

        detectors = []
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

            xCenter = 3600.0*numpy.degrees(centerPupil.getX())
            yCenter = 3600.0*numpy.degrees(centerPupil.getY())
            xMin = 3600.0*numpy.degrees(xmin)
            xMax = 3600.0*numpy.degrees(xmax)
            yMin = 3600.0*numpy.degrees(ymin)
            yMax = 3600.0*numpy.degrees(ymax)
            plateScale = 3600.0*numpy.degrees(plateScale)

            file_handle.write('#detector;%s;%f;%f;%f;%f;%f;%f;%f\n' %
                             (dd.getName(), xCenter, yCenter, xMin, xMax, yMin, yMax, plateScale))

            detector = GalSimDetector(name=dd.getName(), xCenter=xCenter, yCenter=yCenter,
                                      xMin=xMin, yMin=yMin, xMax=xMax, yMax=yMax,
                                      plateScale=plateScale)
        
            detectors.append(detector)

        bandPassFiles, bandPassNames = self._getBandPasses()
        
        self.galSimInterpreter = GalSimInterpreter(detectors=detectors, bandPassNames=bandPassNames,
                                                   bandPassFiles=bandPassFiles, gain=self.gain)
        InstanceCatalog.write_header(self, file_handle)

    def write_catalog(self, *args, **kwargs):
        InstanceCatalog.write_catalog(self, *args, **kwargs)
        self.catalogHasBeenWritten = True

    def write_images(self, isTest=False, nameRoot=None):
        if self.catalogHasBeenWritten is False:
            print "Cannot write GalSim images until you write the GalSim catalog"
            return

        self.galSimInterpreter.writeImages(isTest=isTest, nameRoot=nameRoot)

class GalSimGalaxies(GalSimBase, AstrometryGalaxies, EBVmixin):
    """
    This is a GalSimCatalog class for galaxy components (i.e. objects that are shaped
    like Sersic profiles).

    See the docstring in GalSimBase for explanation of how this class should be used.
    """

    catalog_type = 'galsim_galaxy'
    galsim_type = 'galaxy'
    default_columns = [('galacticAv', 0.1, float),
                       ('galSimType', 'galaxy', (str,6))]


