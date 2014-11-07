"""
The catalog examples in this file will write catalog files that can be read
by galSimInterpreter.py (to be written), which will use galSim to turn them
into an image.
"""

import numpy
from lsst.sims.catalogs.measures.instance import InstanceCatalog, cached
from lsst.sims.coordUtils import CameraCoords, AstrometryGalaxies
from lsst.sims.photUtils import EBVmixin

__all__ = ["GalSimGalaxies"]

def radiansToArcsec(value):
    return 3600.0*numpy.degrees(value)

class GalSimBase(InstanceCatalog, CameraCoords):

    camera = None

    cannot_be_null = ['sedFilepath']

    column_outputs = ['x_pupil', 'y_pupil', 'sedFilepath','magNorm']

    transformations = {'x_pupil':radiansToArcsec,
                       'y_pupil':radiansToArcsec}
    default_formats = {'S':'%s', 'f':'%.9g', 'i':'%i'}

    def get_sedFilepath(self):
        return numpy.array([self.specFileMap[k] if self.specFileMap.has_key(k) else None 
                         for k in self.column_by_name('sedFilename')])

class GalSimGalaxies(GalSimBase, AstrometryGalaxies, EBVmixin):
    catalog_type = 'galsim_galaxy'
    column_outputs = GalSimBase.column_outputs
    column_outputs += ['redshift','positionAngle','galacticAv','galacticRv',
                       'internalAv','internalRv','majorAxis','minorAxis',
                       'sindex']
    default_columns = [('galacticAv', 0.1, float)]
