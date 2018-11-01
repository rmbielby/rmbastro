import mypython as mp
from mypython.ifu import muse
from mypython.ifu import muse_utils as utl
from mypython.ifu import muse_source as msc
import os,sys
from optparse import OptionParser
import numpy as np
import matplotlib.pyplot as plt; plt.ion()
from astropy.io import fits
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
parser.add_option('-s','--nsig',dest='nsig',
                  help='Sigma detection threshold',type='float',
                  default=5)
parser.add_option('-e','--catalogue',dest='catalogue',
                  help='External catalogue',type='string',
                  default=None)
parser.add_option('-g','--segmentation_map',dest='segmentation_map',
                  help='External segmentation_map',type='string',
                  default=None)
parser.add_option('-o','--outfile',dest='outfile',
                  help='External segmentation_map',type='string',
                  default=None)
parser.add_option('-a','--agriffiths',dest='agriffiths',
                  help='Set formatting to A. Griffiths MARZ',type='string',
                  default='N')
parser.add_option('-m','--maskgalacsi',dest='maskgalacsi',
                  help='Mask the gap in the spectrum due to GALACSI',type='string',
                  default='N')

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
# outspec = path+'/'
# if ((os.path.isdir(outspec)) == False): os.mkdir(outspec)
#find sources

print options.nsig

if options.catalogue == None:
    objects = msc.findsources(image, fcube, check=True, output=path,
        nsig=options.nsig, minarea=5., regmask=None, clean=True,spectra =False)
    id      = np.arange(len(objects))+1 # Integer array - identification number of object
else:
    print 'Reading catalogue from ',options.catalogue
    objects =  np.genfromtxt(options.catalogue,unpack=True,usecols=(0,15,16,17,18,20),names="number,x,y,ra,dec,flag")
    id      = objects['number'] # Integer array - identification number of object
    ra      = objects['ra']
    dec     = objects['dec']

nobj = len(objects)

print 'Found {0} objects'.format(nobj)

# print objects.dtype.names
data    = cube_data(fcube)

if options.segmentation_map == None:
    srcmask = fits.getdata(path+'/source.fits')
    # print np.shape(srcmask)
else:
    srcmask = np.array([fits.getdata(options.segmentation_map)])

#loop over detections

naxis3  = len(data.wavec) # Length of rou
x       = objects["x"] # Array of object image x-coords
y       = objects["y"] # Array of object image y-coords
if options.catalogue == None:
    ra,dec,dummy = data.wcsc.wcs_pix2world(x,y,y+0+1,1)
# type   =  String array - each science object is given a type 'P', filled below.
# Initialize flux, variance & sky arrays (sky not mandatory)
intensity,variance,sky = np.zeros((len(id),naxis3)),np.zeros((len(id),naxis3)),np.zeros((len(id),naxis3))

print data.wcsc.wcs.crval[2]
print data.wcsc.wcs.cd[2,2]
print data.wcsc.wcs.crpix[2]

# Fill the flux, variance and sky arrays here.
type = []
flaglim = 256
print 'max flag = ',np.max(objects['flag'])


for i, obj in enumerate(objects[:32]):
    sid = obj['number']
    type.append('P')
    if (obj['flag'] < flaglim):
        print 'ID:',sid,obj['x'],obj['y'],obj['ra'],obj['dec']
        tmpmask=np.zeros(srcmask.shape, dtype=np.bool)
        sel = np.where(srcmask.ravel(srcmask.shape[1]*srcmask.shape[2]) != sid)[0]
        # print len(sel),srcmask.shape[1]*srcmask.shape[2]
        # print np.shape(srcmask)
        # print srcmask.ravel(srcmask.shape[0]*srcmask.shape[1])[sel]
        if len(sel) == srcmask.shape[1]*srcmask.shape[2]:
            print 'Adding object to source mask!!!'
            for j in np.arange(11):
                for k in np.arange(11):
                    radius  = ((j-5)**2+(j-5)**2)**0.5
                    if (radius <= 5) & (obj['y']-5+k < srcmask.shape[1]) & (obj['x']-5+j < srcmask.shape[2]):
                        srcmask[0,np.int(obj['y'])-5+k,np.int(obj['x'])-5+j] = sid
            # print srcmask[0,np.int(obj['y'])-5:np.int(obj['y'])+5,np.int(obj['x'])-5:np.int(obj['y'])+5]
        tmpmask[srcmask == sid] = True
        # print 'tmpmask == True',np.where(tmpmask == True)
        # savename = "{}/id{}.fits".format(outspec, i+1)
        wavec, spec_flx, spec_err, spec_med = utl.cube2spec(fcube, obj['x'], obj['y'], None,
            shape='mask', mask=tmpmask, tovac=True, twod=False, write='temp.fits')
        # Remove !nans
        # print np.nan,np.median(spec_err),np.nanmedian(spec_err)
        if np.isnan(np.median(spec_err)):
            print 'Removing np.nan'
            spec_err[np.isnan(spec_err)] = np.nanmean(spec_err)*1.e6
        if options.maskgalacsi == 'Y':
            print 'Forcing high variance in wavelength gap (GALACSI)'
            galacsi_cent = (5748.9+6047.8)/2.
            galacsi_wdth = (6047.8-5748.9)/2.
            galacsi_sel = np.where(np.abs(wavec-galacsi_cent)<galacsi_wdth)
            spec_err[galacsi_sel] = np.nanmean(spec_err)*1.e6
        # print np.nan,np.median(spec_err),np.nanmedian(spec_err)
        # if np.isnan(np.median(spec_err)):
        # print np.min(wavec),np.nanmin(wavec),np.max(wavec),np.nanmax(wavec)
        # print spec_err
        # f,ax = plt.subplots(1,1,figsize=(12,4))
        # ax.plot(wavec,spec_flx,label='Flux')
        # ax.plot(wavec,spec_err,label='Sigma')
        # # ax.plot(wavec,spec_gle,label='Global sigma')
        # ax.plot(wavec,wavec*0.,color='k',linestyle=':')
        # ax.legend()
        # ax.set_xlim(np.nanmin(wavec),np.nanmax(wavec))
        # plt.pause(0.5)
        if np.isnan(np.nanmedian(spec_flx))==False:
            intensity[i,:] = spec_flx
            variance[i,:]  = spec_err**2
        print '{0:5.0f} {1:7.2f} {2:7.2f} {3:4.0f} {4:7.4f} {5:7.4f} {6:7.4f}'.format(sid,obj['x'], obj['y'],obj['flag'],np.nanmedian(intensity[i,:]),np.nanmedian(variance[i,:]),np.nanmax(variance[i,:]))
        print("done {} !".format(sid))

# print np.shape(sky)
fsel = np.where(objects['flag']<flaglim)[0]
master_var = np.nanmedian(variance[fsel,:],axis=0)
master_sky = np.nanmedian(sky[fsel,:],axis=0)

skymeds = np.nanmedian(sky,axis=1)
varmeds = np.median(variance,axis=1)

# print 'sky & var meds:',skymeds,varmeds
# print np.min(skymeds),np.min(varmeds)

if np.isnan(np.median(varmeds)):
    sel = np.where(np.isnan(varmeds))[0]
    # print sel
    for line in sel:
        variance[line,:] = master_var
        f,ax = plt.subplots(1,1)
        ax.plot(wavec,intensity[line,:],label='Intensity')
        ax.plot(wavec,variance[line,:]**0.5,label='Variance*0.5')
        ax.plot(wavec,sky[line,:],label='Sky')
        ax.legend()
        ax.set_xlim(np.nanmin(wavec),np.nanmax(wavec))
        plt.pause(5)

skymeds = np.nanmedian(sky,axis=1)
varmeds = np.median(variance,axis=1)
# print 'sky & var meds:',skymeds,varmeds
# print np.min(skymeds),np.min(varmeds)
# f,ax = plt.subplots(1,1)
# ax.plot(skymeds,label='SkyMeds')
# ax.plot(varmeds,label='VarMeds')
# ax.legend()
# plt.pause(5)

# temporary section to write out as txt file for Joe
# import sys
# from astropy.table import Table, Column, MaskedColumn
# from astropy.io import ascii
# t = Table(intensity, names=wavec)
# ascii.write(t, './{0}/Extracted_intensity.txt'.format(path))

# Choose a filename:
if options.outfile == None:
    marzfile = './{0}/Extracted_marz.fits'.format(path)
else:
    marzfile = options.outfile
# Set-up the fits file
specmeds = np.nanmean(intensity,axis=1)
specmax = np.nanmax(intensity,axis=1)
print np.shape(objects['flag'])
print np.shape(specmeds)
badobj = np.where((objects['flag']>=flaglim))[0]
# print 'flagged:',badobj,objects['flag']
badobj = np.where((np.isnan(specmeds))| (specmax == 0.))[0]
# badobj = np.where((np.isnan(specmeds)))[0]
# print 'zero-flux:',badobj,specmeds
# print objects['number'][badobj]
for line in badobj:
    objects['flag'][line] = flaglim

goodobj = np.where(objects['flag']<flaglim)[0]
fsel = np.where(objects['flag']<flaglim)[0]

print np.shape(ra)
print np.shape(dec)
print np.shape(fsel)
# print objects['flag'][fsel]
print np.shape(intensity[fsel,:])
print np.shape(intensity)

marz_hdu = fits.HDUList()
marz_hdu.append(fits.ImageHDU(intensity[fsel,:]))
marz_hdu.append(fits.ImageHDU(variance[fsel,:]))
marz_hdu.append(fits.ImageHDU(sky[fsel,:]))
marz_hdu[0].header.set('crval1', data.wcsc.wcs.crval[2]*1.e10)
marz_hdu[0].header.set('crpix1', data.wcsc.wcs.crpix[2])
marz_hdu[0].header.set('cdelt1', data.wcsc.wcs.cd[2,2]*1.e10)
marz_hdu[0].header.set('ctype1', 'WAVE')
marz_hdu[0].header.set('cunit1', 'ANGSTROM')
marz_hdu[0].header.set('naxis1', naxis3)
marz_hdu[1].header.set('extname', 'VARIANCE')
marz_hdu[2].header.set('extname', 'SKY')

# Add in source parameters as fits table

name = []
for i,num in enumerate(id):
    name.append('CatID-{0}'.format(num))
name = np.array(name)
type = np.array(type)

if options.agriffiths == 'Y':
    ra  = ra/180.*np.pi
    dec = dec/180.*np.pi

c1 = fits.Column(name='source_id', format='80A', array=id[fsel])
c2 = fits.Column(name='NAME', format='80A', array=name[fsel])#/180.*np.pi) # Note: these need to be in radians!!!! Maybe not any more?
c3 = fits.Column(name='RA',   format='D', array=ra[fsel])#/180.*np.pi) # Note: these need to be in radians!!!! Maybe not any more?
c4 = fits.Column(name='DEC',  format='D',array=dec[fsel])#/180.*np.pi)
c5 = fits.Column(name='X',    format='J',array=x[fsel])
c6 = fits.Column(name='Y',    format='J', array=y[fsel])
c7 = fits.Column(name='TYPE', format='1A', array=type[fsel])
coldefs = fits.ColDefs([c1, c2, c3, c4, c5,c6])
marz_hdu.append(fits.BinTableHDU.from_columns(coldefs))
marz_hdu[3].header.set('extname', 'FIBRES')

# And write out.
marz_hdu.writeto(marzfile,clobber=True)
