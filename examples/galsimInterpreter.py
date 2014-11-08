import os
import numpy
import galsim
from lsst.sims.photUtils import Sed, Bandpass

imsimband = Bandpass()
imsimband.imsimBandpass()
scale = 0.2
bigfft = galsim.GSParams(maximum_fft_size=10000)
image = None

bandPassDir = os.path.join(os.getenv('THROUGHPUTS_DIR'),'baseline')
bandPassFile = os.path.join(bandPassDir,'total_g.dat')

sedDir = os.getenv('SIMS_SED_LIBRARY_DIR')
dataNeeded = ['x_pupil', 'y_pupil', 'magNorm', 'sedFilepath',
              'redshift', 'positionAngle',
              'galacticAv', 'galacticRv', 'internalAv', 'internalRv',
              'majorAxis', 'minorAxis', 'sindex', 'halfLightRadius']

dataTypes={
    'x_pupil':float, 'y_pupil':float, 'magNorm':float,
    'sedFilepath':(str,126), 'redshift':float, 'positionAngle': float,
    'galacticAv':float, 'galacticRv':float, 'internalAv':float,
    'internalRv':float, 'majorAxis':float, 'minorAxis':float,
    'sindex':float, 'halfLightRadius':float
    }

defaultType = (str,126)

fileName = 'galSim_example.txt'

data = open(fileName,'r')
lines = data.readlines()
data.close()
firstline = lines[0].replace("#","").strip()
ll = numpy.array(firstline.split(', '))


dtype = numpy.dtype(
        [(ww, dataTypes[ww]) if ww in dataTypes else (ww, defaultType) for ww in ll]
        )

data = numpy.genfromtxt(fileName, dtype=dtype, delimiter=', ')

bandPass = galsim.Bandpass(bandPassFile)

xCenter = data['x_pupil'].mean()
yCenter = data['y_pupil'].mean()

xMin = data['x_pupil'].min()-2.0*data['halfLightRadius'].max()
xMax = data['x_pupil'].max()+2.0*data['halfLightRadius'].max()
yMin = data['y_pupil'].min()-2.0*data['halfLightRadius'].max()
yMax = data['y_pupil'].max()+2.0*data['halfLightRadius'].max()

nx = int((xMax-xMin)/scale)
ny = int((yMax-yMin)/scale)

image = galsim.Image(nx,ny)

for (xPupil, yPupil, magNorm, redshift, sedpath, internalAv, internalRv,
     galacticAv, galacticRv, majorAxis, minorAxis, halfLightRadius,
     positionAngle, index) in \
     zip(data['x_pupil'], data['y_pupil'], data['magNorm'], data['redshift'],
         data['sedFilepath'], data['internalAv'], data['internalRv'],
         data['galacticAv'], data['galacticRv'], data['majorAxis'],
         data['minorAxis'], data['halfLightRadius'], data['positionAngle'],
         data['sindex']):


    sedFile = os.path.join(sedDir,sedpath)
        
    obj = galsim.Sersic(n=index, half_light_radius=halfLightRadius,
                        gsparams=bigfft)
        
    obj = obj.shear(q=minorAxis/majorAxis, beta=positionAngle*galsim.radians)
    
    dx=xPupil-xCenter
    dy=yPupil-yCenter
    
    print xPupil-xMin,yPupil-yMin,halfLightRadius
    
    obj = obj.shift(dx, dy)
                         
    sed = Sed()
    sed.readSED_flambda(sedFile)
    fNorm = sed.calcFluxNorm(magNorm, imsimband)
    sed.multiplyFluxNorm(fNorm)
    a_int, b_int = sed.setupCCMab()
    sed.addCCMDust(a_int, b_int, A_v=internalAv, R_v=internalRv)
    sed.redshiftSED(redshift, dimming=False)
    sed.addCCMDust(a_int, b_int, A_v=galacticAv, R_v=galacticRv)
        
    spectrum = galsim.SED(spec = lambda ll:
                                 numpy.interp(ll, sed.wavelen, sed.flambda),
                                 flux_type='flambda')
        
    obj = obj*spectrum
        
    if image is None:
        image = obj.drawImage(bandpass=bandPass, scale=scale, method='real_space')
    else:
        image = obj.drawImage(bandpass=bandPass, scale=scale, image=image, add_to_image=True, method='real_space')
        
image.write('testImage.fits')
