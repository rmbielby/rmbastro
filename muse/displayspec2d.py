import sys
import numpy as np
import matplotlib.pyplot as plt
import mypython as mp
from rmbastro.muse import utilities as rmbutl
from mypython.ifu import muse_utils as utl
from mypython.ifu import muse_source as msc
from astropy.io import fits
import os
import matplotlib.gridspec as gridspec
from scipy.ndimage.filters import gaussian_filter1d
from scipy.ndimage.filters import gaussian_filter
import os,sys
from optparse import OptionParser

uscle = 3.6


#
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

def key_event(e):
    global curr_pos

    if e.key == "right":
        curr_pos = curr_pos + 1
    elif e.key == "left":
        curr_pos = curr_pos - 1
    elif e.key == "n":
        curr_pos = input("Choose a new ID: ")-1
    else:
        return
    curr_pos = curr_pos #% len(plots)
    wave,spec,spec_err,uwrap,simage = get_plots(curr_pos)
    print 'Updating plots with ID: ',curr_pos+1
    ax0.cla()
    ax1.cla()
    ax2.cla()
    ax0.plot(wavec,spec,label=curr_pos+1)
    ax0.plot(wavec,spec_err)
    ax0.plot(wavec,wavec*0.,linestyle=':',color='k')
    ax0.legend()
    ax0.set_xlim(np.min(wavec),np.max(wavec))
    ax0.set_ylabel(r'Flux')
    ax1.imshow(np.transpose(uwrap),aspect='auto',origin='lower',
        extent=(np.min(wave),np.max(wave),0,2*uscle*a[idobj]+1),cmap='RdYlBu_r')
    islice = 1
    while islice < 2*uscle*a[idobj]:
        ax1.plot(wavec,wavec*0.+islice,color='k',linestyle=':')
        islice += 1
    ax1.set_xlim(np.min(wavec),np.max(wavec))
    ax1.set_ylim(0,np.int(2*uscle*a[idobj])+1)
    ax2.imshow(simage,aspect='auto',origin='lower',
        extent=(y[idobj]-uscle*a[idobj],y[idobj]+uscle*a[idobj]+1,x[idobj]-uscle*a[idobj],x[idobj]+uscle*a[idobj]+1),
        cmap='RdYlBu_r',interpolation="bessel")
    print 'Plots updated!'
    fig.canvas.draw()

def get_plots(idobj):
    print 'Loading data for: ',curr_pos+1
    print y[idobj],x[idobj]
    wlu,uwrap,uvar,simage  = rmbutl.cube2uwrp(fcube,  np.int(y[idobj]),np.int(x[idobj]), np.int(uscle*a[idobj]), tovac=True)
    tmpmask=np.zeros(srcmask.shape, dtype=np.bool)
    tmpmask[srcmask == idobj+1] = True
    wavec, spec_flx, spec_err, spec_med = utl.cube2spec(fcube, x[idobj], y[idobj], None,shape='mask', mask=tmpmask, tovac=True)
    spec_flx = gaussian_filter1d(spec_flx,1.6)
    # scut = np.percentile(spec_flx,0.64)
    # spec_flx[spec_flx<scut] = scut
    uwrap = gaussian_filter(uwrap,1.2)
    lcut = np.percentile(uwrap,4.8)
    cut = np.percentile(uwrap,96)
    uwrap[uwrap>cut]  = cut
    uwrap[uwrap<lcut] = lcut
    return wavec,spec_flx,spec_err,uwrap,simage

parser = OptionParser()
parser.add_option('-c','--cubefile',dest='cfile',
                  help='Cube file path',type='string',
                  default='COMBINED_CUBE_FINAL.fits')
parser.add_option('-i','--imagefile',dest='ifile',
                  help='Image file path',type='string',
                  default='COMBINED_IMAGE_FINAL.fits')
parser.add_option('-p','--path',dest='path',
                  help='Image file path',type='string',
                  default='')

try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error ... check usage with -h ";
    sys.exit(1)
# field    = 'J094932+033531'
# print field
path     = options.path
fcube    = '{0}/{1}'.format(path,options.cfile)
image    = '{0}/{1}'.format(path,options.ifile)
ctlgfile = '{0}/catalogue.fits'.format(path)

cat=fits.open(ctlgfile)
nsrc=len(cat[1].data)


# data    = cube_data(fcube)
srcmask = fits.getdata(path+'/source.fits')

x= cat[1].data['x']+0.5
y= cat[1].data['y']+0.5
a= cat[1].data['a']
b= cat[1].data['b']
theta= cat[1].data['theta']

a[a<2] = 2


# naxis3  = len(data.wavec) # Length of rou spectra
# x       = objects["x"] # Array of object image x-coords
# y       = objects["y"] # Array of object image y-coords
# ra,dec,dummy = data.wcsc.wcs_pix2world(x,y,y+0+1,1)



fig = plt.figure(1,figsize=(16,6.4))
gs     = gridspec.GridSpec(2,8)


fig.canvas.mpl_connect('key_press_event', key_event)
print len(x)

curr_pos = 0
idobj    = 1*curr_pos
wavec,spec_flx,spec_err,uwrap,simage = get_plots(curr_pos)
ax0 = plt.subplot(gs[0,3:])
ax0.plot(wavec,spec_flx)
ax0.plot(wavec,spec_err)
ax0.plot(wavec,wavec*0.,linestyle=':',color='k')
ax0.set_xlim(np.min(wavec),np.max(wavec))
ax0.set_ylabel(r'Flux')

ax1 = plt.subplot(gs[1,3:])
ax1.imshow(np.transpose(uwrap),aspect='auto',origin='lower',interpolation='nearest',
    extent=(np.min(wavec),np.max(wavec),0,2*uscle*a[idobj]+1),cmap='RdYlBu_r',vmax=np.percentile(uwrap,96))
islice = 1
while islice < np.int(2*uscle*a[idobj])+1:
    ax1.plot(wavec,wavec*0.+islice,color='k',linestyle=':')
    islice += 1
ax1.set_xlim(np.min(wavec),np.max(wavec))
ax1.set_ylim(0,np.int(2*uscle*a[idobj])+1)
ax1.set_xlabel(r'Wavelength ($\AA$)')

ax2 = plt.subplot(gs[:,:3])
ax2.imshow(simage,aspect='auto',origin='lower',
    extent=(y[idobj]-uscle*a[idobj],y[idobj]+uscle*a[idobj]+1,x[idobj]-uscle*a[idobj],x[idobj]+uscle*a[idobj]+1),
    cmap='RdYlBu_r',interpolation="bessel")

plt.subplots_adjust(left=0.04,right=0.98,top=0.96,hspace=0.02,wspace=0.42)
plt.show()
