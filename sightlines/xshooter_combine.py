import numpy as np
import astropy.io.fits as fits
import os,sys
from optparse import OptionParser
import glob
import matplotlib.pyplot as plt

def readxsh(xshfits):
    hdulist   = fits.open(xshfits)
    try:
        spdata    = hdulist[1].data
        header    = hdulist[1].header
        wave      = spdata.field('wave')[0]*10.
        flux      = spdata.field('flux')[0]
        error     = spdata.field('err_flux')[0]
        try:
            cont      = spdata.field('continuum')[0]
            flux      = flux/cont
        except:
            print 'No continumm'
            cont = flux*0.+1.
    except:
        print hdulist[0].header
        print hdulist[1].header
    return wave,flux,error,cont,header

def write_pyigmfits(pyigmfile,wave,fnorm,error,cont):
    pyigm_hdu = fits.HDUList()
    # Add in source parameters as fits table
    c1 = fits.Column(name='WAVE', format='D',array=wave) # Note: these need to be in radians!!!!
    c2 = fits.Column(name='FLUX', format='D',array=fnorm)
    c3 = fits.Column(name='ERROR', format='D',array=error)
    c4 = fits.Column(name='CONT', format='D', array=cont)
    coldefs = fits.ColDefs([c1, c2, c3, c4])
    pyigm_hdu.append(fits.BinTableHDU.from_columns(coldefs))
    # pyigm_hdu[0].header.set('extname', 'FIBRES')
    # And write out.
    pyigm_hdu.writeto(pyigmfile,clobber=True)

def write_pyigmascii(pyigmfile,wave,fnorm,error,cont):
    x = open(pyigmfile,'w')
    for i in np.arange(len(wave)):
        x.write('{0:10.4f} {1:12.4e} {2:10.4e} {3:10.4e} \n'.format(wave[i],fnorm[i],error[i],cont[i]))
    x.close()

def dispelem(esofile):
    hdulist   = fits.open(esofile)
    return hdulist[0].header['DISPELEM']


parser = OptionParser()
parser.add_option('-q','--qso',dest='qsoname',
                  help='Name of the QSO',type='string',
                  default='')
parser.add_option('-d','--datadir',dest='datadir',
                  help='Diretory hosting all QSO data',type='string',
                  default='/home/rich/Dropbox/muselp/data/qsospectra/')
parser.add_option('-o','--outdir',dest='outdir',
                  help='Diretory for output QSO data-file',type='string',
                  default='/home/rich/Dropbox/MUSELLSLP/data/qsospectra/')

parser.add_option('-l','--llim',dest='llim',
                  help='Wavelength limit for output spectrum',type='string',
                  default=0)
parser.add_option('-u','--ulim',dest='ulim',
                  help='Wavelength limit for output spectrum',type='string',
                  default=1.e8)


try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error ... check usage with -h ";
    sys.exit(1)

print options.qsoname


specfiles = glob.glob('{0}/XSHOOTER/{1}/*fits'.format(options.datadir,options.qsoname))
# specfile = '{0}/XSHOOTER/{1}/{1}.fits'.format(options.datadir,options.qsoname)
# print specfiles

f,ax = plt.subplots(3)
armfiles = ['','','']
for specfile in specfiles:
    if dispelem(specfile) == 'UVB':
        armfiles[0] = specfile
    elif dispelem(specfile) == 'VIS':
        armfiles[1] = specfile
    elif dispelem(specfile) == 'NIR':
        armfiles[2] = specfile

wl = []
fl = []
er = []
ct = []

for i,specfile in enumerate(armfiles):
    wave,flux,err_flux,continuum,hd = readxsh(specfile)
    ax[i].plot(wave,flux)
    ax[i].plot(wave,err_flux)
    ax[i].plot(wave,continuum)
    wl = np.append(wl,wave)
    fl = np.append(fl,flux)
    er = np.append(er,err_flux)
    ct = np.append(ct,continuum)

print options.outdir,options.qsoname,np.float(options.llim),np.float(options.ulim)
sel = np.where((wl>=np.float(options.llim)) & (wl < np.float(options.ulim)))[0]
print sel
write_pyigmascii('{0}/XSHOOTER/{1}.dat'.format(options.outdir,options.qsoname),wl[sel],fl[sel]*ct[sel],er[sel],ct[sel])
plt.pause(5.6)
