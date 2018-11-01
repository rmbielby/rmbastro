import numpy as np
import matplotlib.pyplot as plt
import pyfits
import glob
import subprocess
from optparse import OptionParser
import os,sys
from astropy.coordinates import SkyCoord
from astropy import units as u
parser = OptionParser()

parser.add_option('-b','--dobias',dest='dobias',
       help='Calculate bias 0/1',type='int',
       default=0)
parser.add_option('-f','--doflats',dest='doflats',
       help='Calculate flats 0/1',type='int',
       default=0)
parser.add_option('-t','--standards',dest='standards',
       help='Process standard frames 0/1',type='int',
       default=0)
parser.add_option('-p','--process',dest='process',
       help='Process science frames 0/1',type='int',
       default=0)
parser.add_option('-s','--runsex',dest='runsex',
       help='Run sextractor on science frames 0/1',type='int',
       default=0)
parser.add_option('-a','--posang',dest='posang',
       help='Position angle',type='float',
       default=0)
parser.add_option('-q','--quasar',dest='quasar',
       help='quasar/field',type='string',
       default=0)


try:
    options,args = parser.parse_args(sys.argv[1:])
except:
    print "Error ... check usage with -h ";
    sys.exit(1)


# Instrument params
xsz   = 2048
ysz   = 1034
pxszx = 0.2519
pxszy = 0.2519
PA    = options.posang #135.

# Load mask file
# hdumask = pyfits.open('../flag.fits')
# mask    = 1.-hdumask[0].data
# print 'Mask size: ',np.shape(mask)


# Header changes
rmkeys = ['PROJP1','PROJP3','PV1_1','PV1_2','PV2_1','PV2_2']

if options.dobias == 1:
    biases  = glob.glob('bias/FORS2*.fits')
    nbias   = len(biases)
    print nbias
    biasarr1= np.zeros((ysz,xsz,nbias/2))
    biasarr2= np.zeros((ysz,xsz,nbias/2))
    n1,n2 = 0,0
    for i,bias in enumerate(biases):
        hdulist = pyfits.open(bias)
        biasim  = hdulist[0].data
        print bias,np.shape(biasim)
        if hdulist[0].header['EXTNAME'] == 'CHIP1':
            biasarr1[:,:,n1] = biasim.astype(float)
            n1 += 1
        elif hdulist[0].header['EXTNAME'] == 'CHIP2':
            biasarr2[:,:,n2] = biasim.astype(float)
            n2 += 1

    masterbias1 = np.nanmedian(biasarr1,axis=2)
    pyfits.writeto('bias/master_bias_chip1.fits', masterbias1, hdulist[0].header,
        clobber=True, output_verify='ignore')
    masterbias2 = np.nanmedian(biasarr2,axis=2)
    pyfits.writeto('bias/master_bias_chip2.fits', masterbias2, hdulist[0].header,
        clobber=True, output_verify='ignore')
else:
    try:
        hdulist = pyfits.open('./bias/master_bias_chip1.fits')
        masterbias1 = hdulist[0].data
        hdulist = pyfits.open('./bias/master_bias_chip2.fits')
        masterbias2 = hdulist[0].data
    except:
        print "Couldn't find master bias"
# plt.figure(1)
# plt.imshow(masterbias)
# plt.show()

filters = ['z']

if options.doflats == 1:
    allflats = np.array(glob.glob('flat/FORS2*.fits'))
    band     = np.array(allflats)
    band[:]  = 'z'
    for filter in filters:
        # flats   = np.genfromtxt('flat_{0}.list'.format(filter),unpack=True,dtype=None)
        print np.shape(np.where(band == filter)[0])
        print np.where(band == filter)[0]
        print allflats[np.where(band == filter)[0]]
        flats = allflats[np.where(band == filter)[0]]
        print flats
        nflat   = len(flats)
        flatarr1 = np.zeros((ysz,xsz,nflat/2))
        flatarr2 = np.zeros((ysz,xsz,nflat/2))
        n1,n2 = 0,0
        for i,flat in enumerate(flats):
            print i,flat
            file = glob.glob(flat)[0]
            print file
            hdulist = pyfits.open(file)
            if hdulist[0].header['EXTNAME'] == 'CHIP1':
                flatim  = (hdulist[0].data-masterbias1) #*mask
                flatim  = flatim.astype(float)/np.nanmedian(flatim.astype(float))
                flatarr1[:,:,n1] = flatim
                n1 += 1
            if hdulist[0].header['EXTNAME'] == 'CHIP1':
                flatim  = (hdulist[0].data-masterbias2) #*mask
                flatim  = flatim.astype(float)/np.nanmedian(flatim.astype(float))
                flatarr2[:,:,n2] = flatim
                n2 += 1
        masterflat1 = np.nanmedian(flatarr1,axis=2)
        masterflat1 = masterflat1/np.nanmedian(masterflat1)
        pyfits.writeto('./flat/master_flat_chip1_{0}.fits'.format(filter), masterflat1, hdulist[0].header,clobber=True, output_verify='ignore')
        masterflat2 = np.nanmedian(flatarr2,axis=2)
        masterflat2 = masterflat2/np.nanmedian(masterflat2)
        pyfits.writeto('./flat/master_flat_chip2_{0}.fits'.format(filter), masterflat2, hdulist[0].header,clobber=True, output_verify='ignore')



keys   = ['EXPTIME','ACAMFILT','RA','DEC']

print options.standards
if options.standards == 1:
    # Process the individual images.
    field = 'std'
    print field
    files = glob.glob('{0}/FORS*.fits'.format(field))
    nfiles = len(files)
    print files
    for i,file in enumerate(files):
        print file
        file = file[len(field)+1:]
        print file
        hdulist    = pyfits.open('{0}/{1}'.format(field,file))
        image      = hdulist[0].data
        image      = image.astype(float)
        header0    = hdulist[0].header
        header     = hdulist[0].header

        filter     = header0['HIERARCH ESO INS FILT1 NAME'][0]
        print filter
        if header['EXTNAME'] == 'CHIP1':
            hdulist    = pyfits.open('./flat/master_flat_chip1_{0}.fits'.format(filter))
        if header['EXTNAME'] == 'CHIP2':
            hdulist    = pyfits.open('./flat/master_flat_chip2_{0}.fits'.format(filter))
        masterflat = hdulist[0].data
        mask       = np.array(masterflat)
        mask[masterflat>=0.9] = 1.
        mask[masterflat<0.9]  = 0.
        # Modify header
        coords = SkyCoord(header0['RA'],header0['DEC'],frame='icrs', unit=(u.hourangle, u.deg))
        print coords.ra.degree,coords.dec.degree
        # header.set("ACAMFILT",value=header0['ACAMFILT'],comment=header0.comments["ACAMFILT"])
        header.set("CRVAL1",value=coords.dec.degree,comment=header.comments["CRVAL1"])
        header["CRVAL2"] = coords.ra.degree
        header["CRPIX1"] = ysz/2+1
        header["CRPIX2"] = xsz/2+1
        header["CD1_1"]  = pxszy/3600.*np.cos(PA*np.pi/180.)
        header["CD1_2"]  = -pxszx/3600.*np.sin(PA*np.pi/180.)
        header["CD2_1"]  = pxszx/3600.*np.sin(PA*np.pi/180.)
        header["CD2_2"]  = pxszy/3600.*np.cos(PA*np.pi/180.)
        header["CUNIT1"] = 'deg'
        header["CUNIT2"] = 'deg'
        header.set("CTYPE1",'DEC--TAN',comment='WCS projection type for this axis')
        header.set("CTYPE2",'RA---TAN',comment='WCS projection type for this axis')

        for rmkey in rmkeys:
            if rmkey in header:
                del header[rmkey]
        plt.figure(3)
        plt.hist(image.ravel(),bins=256)

        weightim  = masterflat*1.*mask
        if header['EXTNAME'] == 'CHIP1':
            image     = (image-masterbias1)/(masterflat.astype(float))*mask
        elif header['EXTNAME'] == 'CHIP2':
            image     = (image-masterbias2)/(masterflat.astype(float))*mask
        outfile   = '{0}/p{1}'.format(field,file)
        outweight = '{0}/p{1}_weight.fits'.format(field,file[:-5])
        pyfits.writeto(outfile, image, header,clobber=True, output_verify='ignore')
        pyfits.writeto(outweight, weightim, header,clobber=True, output_verify='ignore')
        fitsheader = subprocess.check_output(['imhead',outfile])
        ahead = open('{0}/p{1}.ahead'.format(field,file[:-5]),'w')
        ahead.write(fitsheader)
        ahead.close()

if options.quasar != '':
    fields = [options.quasar]
else:
    fields = np.genfromtxt('field.list',unpack=True,dtype=None)

if options.process == 1:

    # Process the individual images.
    for field in fields:
        print field
        files = glob.glob('{0}/FORS*.fits'.format(field))
        nfiles = len(files)
        print files
        for i,file in enumerate(files):
            print file
            file = file[len(field)+1:]
            print file
            hdulist    = pyfits.open('{0}/{1}'.format(field,file))
            image      = hdulist[0].data
            image      = image.astype(float)
            header0    = hdulist[0].header
            header     = hdulist[0].header

            filter     = header0['HIERARCH ESO INS FILT1 NAME'][0]
            print filter
            if header['EXTNAME'] == 'CHIP1':
                hdulist    = pyfits.open('./flat/master_flat_chip1_{0}.fits'.format(filter))
            if header['EXTNAME'] == 'CHIP2':
                hdulist    = pyfits.open('./flat/master_flat_chip2_{0}.fits'.format(filter))
            masterflat = hdulist[0].data
            mask       = np.array(masterflat)
            mask[masterflat>=0.9] = 1.
            mask[masterflat<0.9]  = 0.
            mask[image<1000]  = 0.
            # Modify header
            coords = SkyCoord(header0['RA'],header0['DEC'],frame='icrs', unit=(u.hourangle, u.deg))
            print coords.ra.degree,coords.dec.degree
            # header.set("ACAMFILT",value=header0['ACAMFILT'],comment=header0.comments["ACAMFILT"])
            # header.set("CRVAL1",value=coords.dec.degree,comment=header.comments["CRVAL1"])
            # header["CRVAL2"] = coords.ra.degree
            # header["CRPIX1"] = ysz/2+1
            # header["CRPIX2"] = xsz/2+1
            # header["CD1_1"]  = pxszy/3600.*np.cos(PA*np.pi/180.)
            # header["CD1_2"]  = -pxszx/3600.*np.sin(PA*np.pi/180.)
            # header["CD2_1"]  = pxszx/3600.*np.sin(PA*np.pi/180.)
            # header["CD2_2"]  = pxszy/3600.*np.cos(PA*np.pi/180.)
            # header["CUNIT1"] = 'deg'
            # header["CUNIT2"] = 'deg'
            # header.set("CTYPE1",'DEC--TAN',comment='WCS projection type for this axis')
            # header.set("CTYPE2",'RA---TAN',comment='WCS projection type for this axis')

            for rmkey in rmkeys:
                if rmkey in header:
                    del header[rmkey]
            plt.figure(3)
            plt.hist(image.ravel(),bins=256)

            weightim  = masterflat*1.*mask
            if header['EXTNAME'] == 'CHIP1':
                image     = (image-masterbias1)/(masterflat.astype(float))*mask
            elif header['EXTNAME'] == 'CHIP2':
                image     = (image-masterbias2)/(masterflat.astype(float))*mask
            outfile   = '{0}/p{1}'.format(field,file)
            outweight = '{0}/p{1}_weight.fits'.format(field,file[:-5])
            pyfits.writeto(outfile, image, header,clobber=True, output_verify='ignore')
            pyfits.writeto(outweight, weightim, header,clobber=True, output_verify='ignore')
            fitsheader = subprocess.check_output(['imhead',outfile])
            ahead = open('{0}/p{1}.ahead'.format(field,file[:-5]),'w')
            ahead.write(fitsheader)
            ahead.close()

# Create SExtractor catalogues for SCAMP.
if options.runsex == 1:
    for field in fields:
        pfiles = glob.glob('{0}/pFORS2.????-??-??T??:??:??.???.fits'.format(field))
        print pfiles
        for pfile in pfiles:
            ldacfile = '{0}.ldac'.format(pfile[:-5])
            weightfile = '{0}_weight.fits'.format(pfile[:-5])
            sexcommand = ['sextractor', pfile,'-c','sex_ast.conf','-CATALOG_NAME',ldacfile,'-WEIGHT_TYPE','MAP_WEIGHT','-WEIGHT_IMAGE',weightfile]
            print sexcommand
            subprocess.call(sexcommand)
