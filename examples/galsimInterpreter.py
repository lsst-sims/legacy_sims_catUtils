import os
import numpy
import galsim
from lsst.sims.photUtils import Sed, Bandpass


class GalSimInterpreter(object):

    def __init__(self, scale=0.2):
        self.imsimband = Bandpass()
        self.imsimband.imsimBandpass()
        self.scale = scale
        self.bigfft = galsim.GSParams(maximum_fft_size=10000)
        self.sedDir = os.getenv('SIMS_SED_LIBRARY_DIR')
        self.data = None
        self.image = None


    def readCatalog(self,catalogFile):
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
        firstline = lines[0].replace("#","").strip()
        ll = numpy.array(firstline.split(', '))

        dtype = numpy.dtype(
                [(ww, dataTypes[ww]) if ww in dataTypes else (ww, defaultType) for ww in ll]
                )

        self.data = numpy.genfromtxt(catalogFile, dtype=dtype, delimiter=', ')
   
        if self.image is None:
            self.xMin = self.data['x_pupil'].min()-2.0*self.data['halfLightRadius'].max()
            self.xMax = self.data['x_pupil'].max()+2.0*self.data['halfLightRadius'].max()
            self.yMin = self.data['y_pupil'].min()-2.0*self.data['halfLightRadius'].max()
            self.yMax = self.data['y_pupil'].max()+2.0*self.data['halfLightRadius'].max()

            nx = int((self.xMax-self.xMin)/self.scale)
            ny = int((self.yMax-self.yMin)/self.scale)
            
            self.xCenter = 0.5*(self.xMin + self.xMax)
            self.yCenter = 0.5*(self.yMin + self.yMax)

            self.image = galsim.Image(nx,ny)

    def drawCatalog(self, bandPass):
        self.bandPass = galsim.Bandpass(bandPass)
        if self.data is not None:
            for entry in self.data:
                print entry
                if entry['galSimType'] == 'galaxy':
                    print 'drawing'
                    self.drawGalaxy(entry)
                else:
                    print entry['galSimType']

    def drawGalaxy(self,entry):

        sedFile = os.path.join(self.sedDir,entry['sedFilepath'])
        
        obj = galsim.Sersic(n=entry['sindex'], half_light_radius=entry['halfLightRadius'],
                            gsparams=self.bigfft)
        
        obj = obj.shear(q=entry['minorAxis']/entry['majorAxis'], beta=entry['positionAngle']*galsim.radians)
    
        dx=entry['x_pupil']-self.xCenter
        dy=entry['y_pupil']-self.yCenter
    
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
        
        self.image = obj.drawImage(bandpass=self.bandPass, scale=self.scale, image=self.image,
                                  add_to_image=True, method='real_space')
 


bandPass = os.path.join(os.getenv('THROUGHPUTS_DIR'),'baseline','total_g.dat')

gs = GalSimInterpreter()
gs.readCatalog('galsim_example.txt')
gs.drawCatalog(bandPass)
   
gs.image.write('testImage.fits')
