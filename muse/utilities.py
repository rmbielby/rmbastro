def cube2uwrp(cube,x,y,s,write=None,shape='box',helio=0,mask=None,twod=True,tovac=False,idsource=None):

    """
    Extract a 2D spectrum from a cube at position x,y in box or circle of radius s

    If shape = 'mask', then mask is a boolean mask and pixels within it will be extracted form
    argument mask. Mask is a datacube [e.g. from cubex]

    idsource -> if > 0, then only pixels in mask with that ID will be extracted

    helio passes an heliocentric correction in km/s [should be 0 with pipeline v1.2.1]

    twod -> also reconstruct a 2D spec

    tovac -> if true, return wavelengths in vacuum

    write -> output file

    """
    import matplotlib.pyplot as plt
    import numpy as np
    from astropy.io import fits
    from mypython.ifu import muse_utils as utl

    #read the cube
    cubdata,vardata,wcsc,wavec,regi=utl.readcube(cube,helio=helio)
    cubdata=np.nan_to_num(cubdata)

    #if mask extract all True pixels
    if('mask' in shape):
        if(idsource):
            goodpix=np.nonzero(mask == idsource)
        else:
            goodpix=np.nonzero(mask)
        xpix=goodpix[1]
        ypix=goodpix[2]
    else:
        #If user defined region, grab inner pixels
        #cut region of interest according to shape
        xpix=[]
        ypix=[]
        xside=np.arange(x-s-1,x+s+1,1)
        yside=np.arange(y-s-1,y+s+1,1)
        for xx in xside:
            for yy in yside:
                if('box' in shape):
                    if((abs(xx-x) <= s) & (abs(yy-y) <= s)):
                        xpix.append(xx)
                        ypix.append(yy)
                if('circ' in shape):
                    dist=np.sqrt((xx-x)**2+(yy-y)**2)
                    if(dist <= s):
                        xpix.append(xx)
                        ypix.append(yy)

    #Some checks...
    #cbmed=np.median(cubdata,axis=0)
    #cbmed[xpix,ypix]=100000
    #imgplot=plt.imshow(cbmed,origin='lower')
    #imgplot.set_clim(-5,5)
    #plt.show()

    #now make space for the 1d spectrum
    spec_flx=np.zeros(len(wavec))
    spec_var=np.zeros(len(wavec))
    spec_med=np.zeros(len(wavec))

    #if want 2d, prepapre space for it
    #This simulates a slit in the x direction
    #adding up all the flux on the y
    uxpix=np.sort(list(set(xpix)))
    uypix=np.sort(list(set(ypix)))
    npix=len(uxpix)
    nwv=len(wavec)
    twodspec=np.zeros((nwv,npix**2))
    twoderr=np.zeros((nwv,npix**2))

    print np.shape(uxpix),np.shape(uypix)
    print np.shape(xpix),np.shape(ypix)

    #loop over all wavelength to fill in spectrum
    for ii in range(npix):
        for jj in range(npix):
            #add all the pixels in y
            if (uxpix[jj] > -1) & (uypix[jj] > -1):
                twodspec[:,jj+npix*ii] = cubdata[:,uxpix[jj],uypix[ii]]
                twoderr[: ,jj+npix*ii]  = vardata[:,uxpix[jj],uypix[ii]]
    print uxpix[0],uxpix[-1],uypix[0],uypix[-1]
    simage = np.nansum(cubdata[:,uxpix[0]:uxpix[-1],uypix[0]:uypix[-1]],axis=0)
    print np.shape(simage)
    #extract the 2D image with a small buffer around

    #mean in aperture
    totpix=len(xpix)


    #if set, convert to vacuum using airtovac.pro conversion
    if(tovac):
        #save current wave
        wavec=np.array(wavec,dtype=np.float64)
        wave_air=wavec

        sigma2 = (1e4/wavec)**2.
        fact = 1.+5.792105e-2/(238.0185-sigma2)+1.67917e-3/(57.362-sigma2)
        wavec = wavec*fact

    #tested and working
    #fl=open('test.txt','w')
    #for rr in range(len(wavec)):
    #    fl.write("{} {}\n".format(wave_air[rr],wavec[rr]))
    #fl.close()

    #if write, write
    if(write):
        hduflx  = fits.PrimaryHDU(spec_flx) #mean in region
        hduerr  = fits.ImageHDU(spec_err) #associated errors
        hduwav  = fits.ImageHDU(wavec)    #wave
        hdumed  = fits.ImageHDU(spec_med) #median spectrum
        if(twod): #twod
            hdu2flx  = fits.ImageHDU(twodspec)
            hdu2err  = fits.ImageHDU(twoderr)
            hduimg   = fits.ImageHDU(twodimg)
            hdulist = fits.HDUList([hduflx,hduerr,hduwav,hdumed,hdu2flx,hdu2err,hduimg])
        else:
            hdulist = fits.HDUList([hduflx,hduerr,hduwav,hdumed])
        hdulist.writeto(write,clobber=True)

    return wavec, twodspec, twoderr,simage
