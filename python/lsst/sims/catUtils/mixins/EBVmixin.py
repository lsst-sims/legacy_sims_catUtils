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
    
    # Default value for RV in MW
    RVMW = 3.1
    
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

        EBV = self.column_by_name('EBV')
        av = EBV
        try:
            av = numpy.array(EBV/self.RVMW)
        except:
            pass
        return av

    # @cached    
    def get_galacticRv(self):
        """
        Returns galactic RV from self.RVMW which is set to 3.1 by
        default, but can be modified by the user.
        """
        num = 1
        EBV=self.column_by_name('EBV')
        try:
            num = len(EBV)
        except:
            pass
        return self.RVMW * numpy.ones(num)
    

