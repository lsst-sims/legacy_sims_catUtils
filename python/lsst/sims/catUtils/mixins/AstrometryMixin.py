import numpy
import ctypes
import math
import palpy as pal
import lsst.afw.geom as afwGeom
from lsst.afw.cameraGeom import PUPIL, PIXELS, FOCAL_PLANE
from lsst.afw.cameraGeom import SCIENCE
from lsst.sims.catalogs.measures.instance import compound
from lsst.sims.utils import haversine, arcsecFromRadians, radiansFromArcsec, \
                            galacticFromEquatorial, sphericalFromCartesian, \
                            cartesianFromSpherical

from lsst.sims.coordUtils.AstrometryUtils import appGeoFromICRS, observedFromAppGeo
from lsst.sims.coordUtils.AstrometryUtils import observedFromICRS, calculatePupilCoordinates
from lsst.sims.coordUtils.CameraUtils import findChipName, calculatePixelCoordinates
from lsst.sims.coordUtils.CameraUtils import calculateFocalPlaneCoordinates

__all__ = ["AstrometryBase", "AstrometryStars", "AstrometryGalaxies",
           "CameraCoords"]

class AstrometryBase(object):
    """Collection of astrometry routines that operate on numpy arrays"""

    @compound('glon','glat')
    def get_galactic_coords(self):
        """
        Getter for galactic coordinates, in case the catalog class does not provide that

        Reads in the ra and dec from the data base and returns columns with galactic
        longitude and latitude.

        All angles are in radians
        """
        ra=self.column_by_name('raJ2000')
        dec=self.column_by_name('decJ2000')

        glon, glat = galacticFromEquatorial(ra,dec)

        return numpy.array([glon,glat])


    @compound('x_pupil','y_pupil')
    def get_pupilFromSky(self):
        """
        Take an input RA and dec from the sky and convert it to coordinates
        in the pupil.
        """

        raObj = self.column_by_name('raObserved')
        decObj = self.column_by_name('decObserved')

        return calculatePupilCoordinates(raObj, decObj, epoch=self.db_obj.epoch,
                                         obs_metadata=self.obs_metadata)

class CameraCoords(AstrometryBase):
    """Methods for getting coordinates from the camera object"""
    camera = None
    allow_multiple_chips = False #this is a flag which, if true, would allow
                                 #findChipName to return objects that land on
                                 #multiple chips; only the first chip would be
                                 #written to the catalog



    def get_chipName(self):
        """Get the chip name if there is one for each catalog entry"""
        xPupil, yPupil = (self.column_by_name('x_pupil'), self.column_by_name('y_pupil'))
        return findChipName(xPupil=xPupil, yPupil=yPupil, camera=self.camera,
                            allow_multiple_chips=self.allow_multiple_chips)

    @compound('xPix', 'yPix')
    def get_pixelCoordinates(self):
        """Get the pixel positions (or nan if not on a chip) for all objects in the catalog"""
        if not self.camera:
            raise RuntimeError("No camera defined.  Cannot calculate pixel coordinates")
        chipNames = self.column_by_name('chipName')
        xPupil, yPupil = (self.column_by_name('x_pupil'), self.column_by_name('y_pupil'))

        return calculatePixelCoordinates(xPupil = xPupil, yPupil = yPupil, chipNames=chipNames,
                                         camera=self.camera)

    @compound('xFocalPlane', 'yFocalPlane')
    def get_focalPlaneCoordinates(self):
        """Get the focal plane coordinates for all objects in the catalog."""
        xPupil, yPupil = (self.column_by_name('x_pupil'), self.column_by_name('y_pupil'))

        return calculateFocalPlaneCoordinates(xPupil = xPupil, yPupil = yPupil, camera=self.camera)

class AstrometryGalaxies(AstrometryBase):
    """
    This mixin contains a getter for the corrected RA and dec which ignores parallax and proper motion
    """

    @compound('raPhoSim','decPhoSim')
    def get_phoSimCoordinates(self):
        ra = self.column_by_name('raJ2000')
        dec = self.column_by_name('decJ2000')
        return observedFromICRS(ra, dec, includeRefraction = False, obs_metadata=self.obs_metadata,
                                epoch=self.db_obj.epoch)


    @compound('raObserved','decObserved')
    def get_observedCoordinates(self):
        """
        convert mean coordinates in the International Celestial Reference Frame
        to observed coordinates
        """
        ra = self.column_by_name('raJ2000')
        dec = self.column_by_name('decJ2000')
        return observedFromICRS(ra, dec, obs_metadata=self.obs_metadata, epoch=self.db_obj.epoch)


class AstrometryStars(AstrometryBase):
    """
    This mixin contains a getter for the corrected RA and dec which takes account of proper motion and parallax
    """

    def observedStellarCoordinates(self, includeRefraction = True):
        """
        Getter which converts mean coordinates in the International Celestial
        Reference Frame to observed coordinates.
        """

        #TODO
        #are we going to store proper motion in raw radians per year
        #or in sky motion = cos(dec) * (radians per year)
        #PAL asks for radians per year inputs

        pr = self.column_by_name('properMotionRa') #in radians per year
        pd = self.column_by_name('properMotionDec') #in radians per year
        px = self.column_by_name('parallax') #in radians
        rv = self.column_by_name('radialVelocity') #in km/s; positive if receding
        ra = self.column_by_name('raJ2000')
        dec = self.column_by_name('decJ2000')

        return observedFromICRS(ra, dec, pm_ra = pr, pm_dec = pd, parallax = px, v_rad = rv,
                     includeRefraction = includeRefraction, obs_metadata=self.obs_metadata,
                     epoch=self.db_obj.epoch)


    @compound('raPhoSim','decPhoSim')
    def get_phoSimCoordinates(self):
        return self.observedStellarCoordinates(includeRefraction = False)

    @compound('raObserved','decObserved')
    def get_observedCoordinates(self):
        return self.observedStellarCoordinates()
