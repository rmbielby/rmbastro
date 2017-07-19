import mypython as mp
from mypython.ifu import muse
from mypython.ifu import muse_utils as utl
from mypython.ifu import muse_source as msc
import os,sys
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt; plt.ion()
from astropy.io import fits
import os
from PIL import Image
import scipy.misc

class cube_data(object):
    def __init__(self, cubename, helio=0.0):
        """
        chunk([z0, z0_1], z_lower, z_upper)
        Class for fitting chunk of data
        """
        cubdata,vardata,wcsc,wavec,regi = utl.readcube(fcube,helio=helio)
        self.cubdata=cubdata
        self.vardata=vardata
        self.wcsc=wcsc
        self.wavec=wavec
        self.regi=regi


parser = OptionParser()
parser.add_option('-c','--cubefile',dest='fcube',
                  help='Cube file path',type='string',
                  default='')
parser.add_option('-i','--image',dest='image',
                  help='Image file path',type='string',
                  default='')

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error ... check usage with -h ";
    sys.exit(1)



fcube   = options.fcube
print fcube
image   = options.image
path    = '/'.join(image.split('/')[:-1])
print path
# outspec = path+'/Spectra'
# if ((os.path.isdir(outspec)) == False): os.mkdir(outspec)
#find sources
objects = msc.findsources(image, fcube, check=True, output=path,
    nsig=3., minarea=5., regmask=None, clean=True, spectra=False)
nobj = len(objects)

# print objects.dtype.names
data    = cube_data(fcube)
srcmask = fits.getdata(path+'/source.fits')
#loop over detections

naxis3  = len(data.wavec) # Length of rou spectra
id      = np.arange(nobj)+1 # Integer array - identification number of object
x       = objects["x"] # Array of object image x-coords
y       = objects["y"] # Array of object image y-coords
ra,dec,dummy = data.wcsc.wcs_pix2world(x,y,y+0+1,1)
# type   =  String array - each science object is given a type 'P', filled below.
# Initialize flux, variance & sky arrays (sky not mandatory)
intensity,variance,sky = np.zeros((len(id),naxis3)),np.zeros((len(id),naxis3)),np.zeros((len(id),naxis3))

print data.wcsc.wcs.crval[2]
print data.wcsc.wcs.cd[2,2]
print data.wcsc.wcs.crpix[2]

# Fill the flux, variance and sky arrays here.
type = []
for i, obj in enumerate(objects):
# for i in np.arange(30,32):
#     obj = objects[i]
    # world = data.wcsc.wcs_pix2world([x[0],y[0]], 1)
    # print world
    tmpmask=np.zeros(srcmask.shape, dtype=np.bool)
    tmpmask[srcmask == i+1] = True
    # savename = "{}/id{}.fits".format(outspec, i+1)
    wavec, spec_flx, spec_err, spec_med = utl.cube2spec(fcube, obj['x'], obj['y'], None,
        shape='mask', mask=tmpmask, tovac=True, twod=False, write='temp.fits')
    type.append('P')
    # Remove !nans
    # print np.nan,np.median(spec_err),np.nanmedian(spec_err)
    print 'ID:',id[i]
    if np.isnan(np.median(spec_err)):
        print 'Removing np.nan'
        spec_err[np.isnan(spec_err)] = np.nanmedian(spec_err)*1.e6
    # print np.nan,np.median(spec_err),np.nanmedian(spec_err)
    if np.isnan(np.median(spec_err)):
        print np.min(wavec),np.nanmin(wavec),np.max(wavec),np.nanmax(wavec)
        print spec_err
        # f,ax = plt.subplots(1,1,figsize=(18,6))
        # ax.plot(wavec,spec_flx,label='Flux')
        # ax.plot(wavec,spec_err,label='Sigma')
        # # ax.plot(wavec,spec_gle,label='Global sigma')
        # ax.plot(wavec,wavec*0.,color='k',linestyle=':')
        # ax.legend()
        # ax.set_xlim(np.nanmin(wavec),np.nanmax(wavec))
        # plt.pause(5)
    intensity[i,:] = spec_flx
    variance[i,:]  = spec_err**2
    # sky[i,:]       = spec_gle
    print '{0:5.0f} {1:7.4f} {2:7.4f}'.format(i,np.median(spec_flx),np.median(spec_err))
    print("done {} spectra!".format(i+1))

# print np.shape(sky)
master_var = np.nanmedian(variance,axis=0)
master_sky = np.nanmedian(sky,axis=0)
skymeds = np.nanmedian(sky,axis=1)

varmeds = np.median(variance,axis=1)
if np.isnan(np.median(varmeds)):
    sel = np.where(np.isnan(varmeds))[0]
    print sel
    for line in sel:
        variance[line,:] = master_var
        # f,ax = plt.subplots(1,1)
        # ax.plot(wavec,intensity[line,:],label='Intensity')
        # ax.plot(wavec,variance[line,:]**0.5,label='Variance*0.5')
        # ax.plot(wavec,sky[line,:],label='Sky')
        # ax.legend()
        # ax.set_xlim(np.nanmin(wavec),np.nanmax(wavec))
        # plt.pause(5)

skymeds = np.nanmedian(sky,axis=1)
varmeds = np.median(variance,axis=1)
print skymeds,varmeds
print np.min(skymeds),np.min(varmeds)
# f,ax = plt.subplots(1,1)
# ax.plot(skymeds,label='SkyMeds')
# ax.plot(varmeds,label='VarMeds')
# ax.legend()
# plt.pause(5)

# temporary section to write out as txt file for Joe
import sys
from astropy.table import Table, Column, MaskedColumn
from astropy.io import ascii
t = Table(intensity, names=wavec)
ascii.write(t, '{0}/ExtractedSpectra_intensity.txt'.format(path))

# Choose a filename:
marzfile = '{0}/ExtractedSpectra_marz.fits'.format(path)

# Set-up the fits file
marz_hdu = fits.HDUList()
marz_hdu.append(fits.ImageHDU(intensity))
marz_hdu.append(fits.ImageHDU(variance))
marz_hdu.append(fits.ImageHDU(sky))
marz_hdu[0].header.set('crval1', data.wcsc.wcs.crval[2]*1.e10)
marz_hdu[0].header.set('crpix1', data.wcsc.wcs.crpix[2])
marz_hdu[0].header.set('cdelt1', data.wcsc.wcs.cd[2,2]*1.e10)
marz_hdu[0].header.set('ctype1', 'WAVE')
marz_hdu[0].header.set('cunit1', 'ANGSTROM')
marz_hdu[0].header.set('naxis1', naxis3)
marz_hdu[1].header.set('extname', 'VARIANCE')
marz_hdu[2].header.set('extname', 'SKY')

# Add in source parameters as fits table
c1 = fits.Column(name='source_id', format='80A', array=id)
c2 = fits.Column(name='RA', format='D', array=ra/180.*np.pi) # Note: these need to be in radians!!!!
c3 = fits.Column(name='DEC', format='D',array=dec/180.*np.pi)
c4 = fits.Column(name='X', format='J',array=x)
c5 = fits.Column(name='Y', format='J', array=y)
c6 = fits.Column(name='TYPE', format='1A', array=type)
coldefs = fits.ColDefs([c1, c2, c3, c4, c5,c6])
marz_hdu.append(fits.BinTableHDU.from_columns(coldefs))
marz_hdu[3].header.set('extname', 'FIBRES')

# And write out.
marz_hdu.writeto(marzfile,clobber=True)
