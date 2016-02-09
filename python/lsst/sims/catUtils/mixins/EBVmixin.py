import pyfits
import math
import numpy
import os

from lsst.sims.catalogs.measures.instance import cached
from lsst.sims.photUtils import EBVmap, EBVbase

__all__ = ["EBVmixin"]


class EBVmixin(EBVbase):
    """
    This mixin class contains the getters which a catalog object will use to call
    calculateEbv in the EBVbase class
    """

    #and finally, here is the getter
    @cached
    def get_EBV(self):
        """
        Getter for the InstanceCatalog framework
        """

        galacticCoordinates=numpy.array([self.column_by_name('glon'),self.column_by_name('glat')])

        EBV_out=numpy.array(self.calculateEbv(galacticCoordinates=galacticCoordinates,interp=True))
        return EBV_out

    @cached
    def get_galacticAv(self):
        """
        Getter to return galactic(Milky Way) Av based on EBV values
        from getter (reading dustmaps) and the RV value from galacticRv
        """
        return self.column_by_name('EBV')*self.column_by_name('galacticRv')


