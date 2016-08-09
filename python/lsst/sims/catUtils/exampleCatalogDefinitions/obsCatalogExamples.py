"""Instance Catalog"""
import numpy
import lsst.utils
from lsst.sims.catalogs.definitions import InstanceCatalog
from lsst.sims.catUtils.mixins import AstrometryStars, CameraCoords, PhotometryStars
from lsst.obs.lsstSim.utils import loadCamera

__all__ = ["ObsStarCatalogBase"]

class ObsStarCatalogBase(InstanceCatalog, AstrometryStars, PhotometryStars, CameraCoords):
    comment_char = ''
    camera = loadCamera(lsst.utils.getPackageDir('obs_lsstSim'))
    catalog_type = 'obs_star_cat'
    column_outputs = ['uniqueId', 'raObserved', 'decObserved', 'lsst_r', 'sigma_lsst_r',
                      'chipName', 'xPix', 'yPix']
    default_formats = {'S':'%s', 'f':'%.8f', 'i':'%i'}
    transformations = {'raObserved':numpy.degrees, 'decObserved':numpy.degrees}
