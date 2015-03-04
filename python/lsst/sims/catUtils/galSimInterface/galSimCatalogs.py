"""
The catalog classes in this file use the InstanceCatalog infrastructure to construct
FITS images for each detector-filter combination on a simulated camera.  This is done by
instantiating the class GalSimInterpreter.  This GalSimInterpreter is the class which
actually generates the FITS images.  As the GalSim Instance Catalogs are iterated over,
each object in the catalog is passed to tothe GalSimInterpeter, which adds the object
to the appropriate FITS images.  The user can then write the images to disk by calling
the write_images method in the GalSim Instance Catalog.

Objects are passed to the GalSimInterpreter by the get_fitsFiles getter function, which
adds a column to the InstanceCatalog indicating which detectors' FITS files contain each image.

Note: because each GalSim Instance Catalog has its own GalSimInterpreter, it means
that each GalSimInterpreter will only draw FITS images containing one type of object
(whatever type of object is contained in the GalSim Instance Catalog).  If the user
wishes to generate FITS images containing multiple types of object, the method
copyGalSimInterpreter allows the user to pass the GalSimInterpreter from one
GalSim Instance Catalog to another (so, the user could create a GalSim Instance
Catalog of stars, generate that InstanceCatalog, then create a GalSim Instance Catalog
of galaxies, pass the GalSimInterpreter from the star catalog to this new catalog,
and thus create FITS images that contain both stars and galaxies).
"""

import numpy
import os
import eups
import copy
from lsst.sims.utils import radiansToArcsec
from lsst.sims.catalogs.measures.instance import InstanceCatalog, cached, is_null
from lsst.sims.coordUtils import CameraCoords, AstrometryGalaxies, AstrometryStars
from lsst.sims.catUtils.galSimInterface import GalSimInterpreter, GalSimDetector
from lsst.sims.photUtils import EBVmixin, Sed, Bandpass
import lsst.afw.cameraGeom.testUtils as camTestUtils
import lsst.afw.geom as afwGeom
from lsst.afw.cameraGeom import PUPIL, PIXELS, FOCAL_PLANE

__all__ = ["GalSimGalaxies", "GalSimAgn", "GalSimStars"]

class GalSimBase(InstanceCatalog, CameraCoords):
    """
    The catalog classes in this file use the InstanceCatalog infrastructure to construct
    FITS images for each detector-filter combination on a simulated camera.  This is done by
    instantiating the class GalSimInterpreter.  GalSimInterpreter is the class which
    actually generates the FITS images.  As the GalSim Instance Catalogs are iterated over,
    each object in the catalog is passed to tothe GalSimInterpeter, which adds the object
    to the appropriate FITS images.  The user can then write the images to disk by calling
    the write_images method in the GalSim Instance Catalog.

    Objects are passed to the GalSimInterpreter by the get_fitsFiles getter function, which
    adds a column to the InstanceCatalog indicating which detectors' FITS files contain each image.

    Note: because each GalSim Instance Catalog has its own GalSimInterpreter, it means
    that each GalSimInterpreter will only draw FITS images containing one type of object
    (whatever type of object is contained in the GalSim Instance Catalog).  If the user
    wishes to generate FITS images containing multiple types of object, the method
    copyGalSimInterpreter allows the user to pass the GalSimInterpreter from one
    GalSim Instance Catalog to another (so, the user could create a GalSim Instance
    Catalog of stars, generate that InstanceCatalog, then create a GalSim Instance Catalog
    of galaxies, pass the GalSimInterpreter from the star catalog to this new catalog,
    and thus create FITS images that contain both stars and galaxies; see galSimCompoundGenerator.py
    in the examples/ directory of sims_catUtils for an example).

    This class (GalSimBase) is the base class for all GalSim InstanceCatalogs.  Daughter
    classes of this calss need to behave like ordinary InstanceCatalog daughter classes
    with the following exceptions:

    1) If they re-define column_outputs, they must be certain to include the column
    'fitsFiles', as the getter for this column (defined in this class) calls all of the
    GalSim image generation infrastructure

    2) Daughter classes of this class must define a member variable galsim_type that is either
    'sersic' or 'pointSource'.  This variable tells the GalSimInterpreter how to draw the
    object (to allow a different kind of image profile, define a new method in the GalSimInterpreter
    class similar to drawPoinSource and drawSersic)

    3) The variables bandpass_names (a list of the form ['u', 'g', 'r', 'i', 'z', 'y']),
    bandpass_directory, and bandpass_root should be defined to tell the GalSim Instance Catalog
    where to find the files defining the bandpasses to be used for these FITS files.
    The GalSim Instance Catalog will look for bandpass files in files with the names

    for bpn in bandpass_names:
        name = self.bandpass_directory+'/'+self.bandpass_root+'_'+bpn+'.dat'

    4) Telescope parameters such as exposure time, area, and gain are controled by the
    GalSim InstanceCatalog member variables exposure_time (in s), effective_area (in cm^2),
    and gain (in photons per ADU)

    Daughter classes of GalSimBase will generate both FITS images for all of the detectors/filters
    in their corresponding cameras and Instance Catalogs listing all of the objects
    contained in those images.  The catalog is written using the normal write_catalog()
    method provided for all InstanceClasses.  The FITS files are drawn using the write_images()
    method that is unique to GalSim Instance Catalogs.  The FITS file will be named something like:

    DetectorName_FilterName.fits

    (a typical LSST fits file might be R_0_0_S_1_0_y.fits)

    Note: If you call write_images() before calling write_catalog(), nothing will happen.
    Objects are only added to the GalSimInterpreter in the course of writing
    the Instance Catalog.
    """

    #This is sort of a hack; it prevents findChipName in coordUtils from dying
    #if an object lands on multiple science chips.
    allow_multiple_chips = True

    #There is no point in writing things to the InstanceCatalog that do not have SEDs and/or
    #do not land on any detectors
    cannot_be_null = ['sedFilepath', 'fitsFiles']

    column_outputs = ['galSimType', 'uniqueId', 'chipName', 'x_pupil', 'y_pupil', 'sedFilepath',
                      'majorAxis', 'minorAxis', 'sindex', 'halfLightRadius',
                      'positionAngle','fitsFiles']

    transformations = {'x_pupil':radiansToArcsec,
                       'y_pupil':radiansToArcsec,
                       'halfLightRadius':radiansToArcsec}

    default_formats = {'S':'%s', 'f':'%.9g', 'i':'%i'}

    #This is used as the delimiter because the names of the detectors printed in the fitsFiles
    #column contain both ':' and ','
    delimiter = ';'

    #default variables telling the InstanceCatalog where to find bandpass data
    bandpass_names = ['u','g','r','i','z','y']
    bandpass_directory = eups.productDir('throughputs')
    bandpass_root = 'baseline/total'

    #This member variable will define a PSF to convolve with the sources.
    #See the classes PSFbase, ExampleGaussianPSF, and ExampleOpticalPSF in
    #galSimUtilities.py for more information
    PSF = None

    #This member variable can store a GalSim noise model instantiation
    #which will be applied to the FITS images by calling add_noise()
    noise = None

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

    #This is just a place holder for the camera object associated with the Instance Catalog.
    #If you want to assign a different camera, you can do so immediately after instantiating this class
    camera = camTestUtils.CameraWrapper().camera


    uniqueSeds = {} #a cache for un-normalized SED files, so that we do not waste time on I/O

    hasBeenInitialized = False

    galSimInterpreter = None #the GalSimInterpreter instantiation for this catalog
                             #This class is either passed in from another catalog using
                             #copyGalSimInterpreter, or initialized in the write_header method

    totalDrawings = 0
    totalObjects = 0

    def _initializeGalSimCatalog(self):
        """
        Initializes an empy list of objects that have already been drawn to FITS images.
        We do not want to accidentally draw an object twice.

        Objects are stored based on their uniqueId values.
        """
        self.objectHasBeenDrawn = []
        self.initializeGalSimInterpreter()
        self.hasBeenInitialized = True

    def get_sedFilepath(self):
        """
        Maps the name of the SED as stored in the database to the file stored in
        sims_sed_library
        """
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
                    #we have already read in this file; no need to do it again
                    sed = Sed(wavelen=self.uniqueSeds[sedName].wavelen,
                              flambda=self.uniqueSeds[sedName].flambda,
                              fnu=self.uniqueSeds[sedName].fnu,
                              name=self.uniqueSeds[sedName].name)
                else:
                    #load the SED of the object
                    sed = Sed()
                    sedFile = os.path.join(sedDir, sedName)
                    sed.readSED_flambda(sedFile)

                    flambdaCopy = copy.deepcopy(sed.flambda)

                    #If the SED is zero inside of the bandpass, GalSim raises an error.
                    #This sets a minimum flux value of 1.0e-30 so that the SED is never technically
                    #zero inside of the bandpass.
                    sed.flambda = numpy.array([ff if ff>1.0e-30 else 1.0e-30 for ff in flambdaCopy])
                    sed.fnu = None

                    #copy the unnormalized file to uniqueSeds so we don't have to read it in again
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


    def get_fitsFiles(self):
        """
        This getter returns a column listing the names of the detectors whose corresponding
        FITS files contain the object in question.  The detector names will be separated by a '//'

        This getter also passes objects to the GalSimInterpreter to actually draw the FITS
        images.
        """
        objectNames = self.column_by_name('uniqueId')
        xPupil = self.column_by_name('x_pupil')
        yPupil = self.column_by_name('y_pupil')
        halfLight = self.column_by_name('halfLightRadius')
        minorAxis = self.column_by_name('minorAxis')
        majorAxis = self.column_by_name('majorAxis')
        positionAngle = self.column_by_name('positionAngle')
        sindex = self.column_by_name('sindex')

        #correct the SEDs for redshift, dust, etc.  Return a list of Sed objects as defined in
        #sims_photUtils/../../Sed.py
        sedList = self._calculateGalSimSeds()

        if self.hasBeenInitialized is False:
            #start keeping track of the names of the objects that have already been drawn
            self._initializeGalSimCatalog()

        output = []
        for (name, xp, yp, hlr, minor, major, pa, ss, sn) in \
            zip(objectNames, xPupil, yPupil, halfLight, minorAxis, majorAxis, positionAngle,
            sedList, sindex):

            if ss is None or name in self.objectHasBeenDrawn:
                #do not draw objects that have no SED or have already been drawn
                output.append(None)
                if name in self.objectHasBeenDrawn:
                    #15 December 2014
                    #This should probably be an error.  However, something is wrong with
                    #the SQL on fatboy such that it does return the same objects more than
                    #once (at least in the case of stars).  Yusra is currently working to fix
                    #the problem.  Until then, this will just warn you that the same object
                    #appears twice in your catalog and will refrain from drawing it the second
                    #time.
                    print 'Trying to draw %s more than once ' % str(name)

            else:

                self.objectHasBeenDrawn.append(name)

                #actually draw the object
                detectorsString = self.galSimInterpreter.drawObject(galSimType=self.galsim_type,
                                                  sindex=sn, minorAxis=minor,
                                                  majorAxis=major, positionAngle=pa, halfLightRadius=hlr,
                                                  xPupil=xp, yPupil=yp, sed=ss)

                output.append(detectorsString)

        return numpy.array(output)

    def _getBandpasses(self):
        """
        Create a list of paths to the files containing bandpass data.

        returns the list of paths and the list self.bandpass_names (which are just
        tags identifying each bandpass a la ['u', 'g', 'r', 'i', 'z', 'y'])
        """
        bandpassFiles = []
        for bpn in self.bandpass_names:
            name = self.bandpass_directory+'/'+self.bandpass_root+'_'+bpn+'.dat'
            bandpassFiles.append(name)

        return bandpassFiles, self.bandpass_names

    def copyGalSimInterpreter(self, otherCatalog):
        """
        Copy the camera, GalSimInterpreter, from another GalSim Instance Catalog
        so that multiple types of object (stars, AGN, galaxy bulges, galaxy disks, etc.)
        can be drawn on the same FITS files.

        @param [in] otherCatalog is another GalSim Instance Catalog that already has
        an initialized GalSimInterpreter

        See galSimCompoundGenerator.py in the examples/ directory of sims_catUtils for
        an example of how this is used.
        """
        self.camera = otherCatalog.camera
        self.galSimInterpreter = otherCatalog.galSimInterpreter

        #set the PSF to the current PSF; in this way, compound FITS files do not
        #have to have the same PSF for all types of objects (though I'm not sure
        #that is actually physical)
        self.galSimInterpreter.setPSF(self.PSF)

    def write_header(self, file_handle):

        if not self.hasBeenInitialized:
            self._initializeGalSimCatalog()

        for detector in self.galSimInterpreter.detectors:
            file_handle.write('#detector;%s;%f;%f;%f;%f;%f;%f;%f\n' %
                                 (detector.name, detector.xCenter, detector.yCenter, detector.xMin,
                                  detector.xMax, detector.yMin, detector.yMax, detector.plateScale))

        InstanceCatalog.write_header(self, file_handle)
    

    def initializeGalSimInterpreter(self):
        """
        Overwrite the write_header method from InstanceCatalog.

        This adds information about the camera to the InstanceCatalog header.
        For each detector it writes:

        detector name
        center coordinates in arc seconds
        xmin and xmax in arc seconds
        ymin and ymax in arc seconds
        plateScale (arc seconds per pixel)

        This method will also call _getBandpasses to construct the paths to
        the files containing the bandpass data and construct the GalSimInterpreter
        if it is still None
        """

        if self.galSimInterpreter is None:

            #This list will contain instantiations of the GalSimDetector class
            #(see galSimInterpreter.py), which stores detector information in a way
            #that the GalSimInterpreter will understand
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

                detector = GalSimDetector(name=dd.getName(), xCenter=xCenter, yCenter=yCenter,
                                          xMin=xMin, yMin=yMin, xMax=xMax, yMax=yMax,
                                          plateScale=plateScale)

                detectors.append(detector)

            bandpassFiles, bandpassNames = self._getBandpasses()

            self.galSimInterpreter = GalSimInterpreter(detectors=detectors, bandpassNames=bandpassNames,
                                                       bandpassFiles=bandpassFiles, gain=self.gain)

            self.galSimInterpreter.setPSF(PSF=self.PSF)

    def add_noise(self):
        """
        Adds the noise model stored in self.noise to the images stored
        in the GalSimInterpreter
        """

        if self.noise is not None:
            self.galSimInterpreter.addNoise(noiseWrapper = self.noise, obs_metadata=self.obs_metadata)

    def write_images(self, nameRoot=None):
        """
        Writes the FITS images associated with this InstanceCatalog.

        Cannot be called before write_catalog is called.

        @param [in] nameRoot is an optional string prepended to the names
        of the FITS images.  The FITS images will be named

        nameRoot_DetectorName_FilterName.fits

        (e.g. myImages_R_0_0_S_1_1_y.fits for an LSST-like camera with
        nameRoot = 'myImages')
        """
        namesWritten = self.galSimInterpreter.writeImages(nameRoot=nameRoot)
        return namesWritten

class GalSimGalaxies(GalSimBase, AstrometryGalaxies, EBVmixin):
    """
    This is a GalSimCatalog class for galaxy components (i.e. objects that are shaped
    like Sersic profiles).

    See the docstring in GalSimBase for explanation of how this class should be used.
    """

    catalog_type = 'galsim_galaxy'
    galsim_type = 'sersic'
    default_columns = [('galacticAv', 0.1, float),
                       ('galSimType', 'sersic', (str,6))]

class GalSimAgn(GalSimBase, AstrometryGalaxies, EBVmixin):
    """
    This is a GalSimCatalog class for AGN.

    See the docstring in GalSimBase for explanation of how this class should be used.
    """
    catalog_type = 'galsim_agn'
    galsim_type = 'pointSource'
    default_columns = [('galacticAv', 0.1, float),
                      ('galSimType', 'pointSource', (str,11)),
                      ('majorAxis', 0.0, float),
                      ('minorAxis', 0.0, float),
                      ('sindex', 0.0, float),
                      ('positionAngle', 0.0, float),
                      ('halfLightRadius', 0.0, float),
                      ('internalAv', 0.0, float),
                      ('internalRv', 0.0, float)]

class GalSimStars(GalSimBase, AstrometryStars, EBVmixin):
    """
    This is a GalSimCatalog class for stars.

    See the docstring in GalSimBase for explanation of how this class should be used.
    """
    catalog_type = 'galsim_stars'
    galsim_type = 'pointSource'
    default_columns = [('galacticAv', 0.1, float),
                      ('galSimType', 'pointSource', (str,11)),
                      ('internalAv', 0.0, float),
                      ('internalRv', 0.0, float),
                      ('redshift', 0.0, float),
                      ('majorAxis', 0.0, float),
                      ('minorAxis', 0.0, float),
                      ('sindex', 0.0, float),
                      ('positionAngle', 0.0, float),
                      ('halfLightRadius', 0.0, float)]
