def makeGaussian(size, fwhm = 3, center=None):
    import numpy as np
    """ Make a square gaussian kernel.
    size is the length of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """

    x = np.arange(0, size, 1, float)
    y = x[:,np.newaxis]

    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]

    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)


def completeness(filename,zeropoint,seeing,nsource=100.,mags=None):
    import sys
    import numpy as np
    import astropy.io.fits as fits
    import matplotlib.pyplot as plt
    import pysex
    import astropy.io.ascii as text
    from datetime import datetime
    # nsource = 32.
    thsz = 32

    print 'Running completeness sims on:\n{0}\nwith seeing {1}".'.format(filename,seeing)
    print datetime.now().time()
    hdulist  = fits.open(filename)
    image    = hdulist[0].data
    imsz     = np.shape(image)
    try:
        pixscale = np.abs(hdulist[0].header['CDELT1'])
    except:
        pixscale = np.abs(hdulist[0].header['CD1_1'])

    seeing_pix = seeing/3600./pixscale
    # print seeing_pix

    master_source = makeGaussian(2.*thsz,fwhm=seeing_pix)
    master_source = master_source/np.nansum(master_source)

    # f1,ax = plt.subplots(1,1)
    if len(mags) < 2:
        mags = np.arange(22.,32.,0.5)

    frac = mags*0
    ncalc = np.float(nsource*len(mags))
    print 'Running {0} x {1} sims'.format(nsource,len(mags))
    for i,mag in enumerate(mags):
        # print 'Running mag = ',mag
        ndet = 0
        flux  = 10.**(-0.4*(mag - zeropoint))
        for j in np.arange(nsource):
            sys.stdout.write("\rDone {0:8.2f}%".format((i*nsource+j)/(ncalc)*100.))
            sys.stdout.flush()
            thsig = 3.
            thstd = 0.
            thmed = 0.
            ntry = 0
            while (thstd == 0) | (thmed == 0) | (thsig >= 3.):
                xpix = np.int(np.random.rand()*(imsz[0]-thsz*2.)+thsz)
                ypix = np.int(np.random.rand()*(imsz[1]-thsz*2.)+thsz)
                thumb = image[xpix-thsz:xpix+thsz,ypix-thsz:ypix+thsz]
                thstd = np.nanstd(thumb)
                thmed = np.nanmedian(thumb)
                thsig = np.nanmean(thumb[thsz-2:thsz+3,thsz-2:thsz+3])/thstd
                ntry += 1
                if ntry > 10: thsig = thsig/2.
                if ntry > 20: thsig = thsig/4.
            # print xpix,ypix,np.nanmedian(thumb),np.std(thumb),thumb[thsz,thsz]/thstd
            thwsrc = thumb + flux*master_source
            cat = pysex.run(thwsrc, params=['X_IMAGE', 'Y_IMAGE', 'FLUX_APER','MAG_AUTO'],
                conf_args={'DETECT_THRESH':2.0,'PHOT_APERTURES':5,'MAG_ZEROPOINT':zeropoint})
            # print cat
            # print '------------------------'
            xpos,ypos,mag = [],[],[]
            try:
                for k,x in enumerate(cat['X_IMAGE']):
                    xpos.append(x)
                    ypos.append(cat['Y_IMAGE'][k])
                    mag.append(cat['MAG_AUTO'][k])
                xpos,ypos = np.array(xpos),np.array(ypos)
                offsets = ((xpos-(thsz+1))**2+(ypos-(thsz+1))**2)**0.5
                if np.nanmin(offsets) < 1.6*seeing_pix:
                    ndet += 1.
                # print xpos[np.argmin(offsets)],ypos[np.argmin(offsets)],mag[np.argmin(offsets)]
            except:
                ndet = ndet
                # print 'Source not detected!'

            # ax.imshow(thwsrc,origin='lower',vmin = np.nanpercentile(thwsrc,16),vmax = np.nanpercentile(thwsrc,98))
            # plt.pause(0.8)
        frac[i] = ndet/nsource
    print '\n============================================================\n'
    return mags,frac

# import numpy as np
# mags,frac = completeness('data/wht-acam/fields/HB890232-042/HB890232-042_r.fits',26.1719,0.8)
# for i in np.arange(len(mags)):
#     print mags[i],frac[i]
