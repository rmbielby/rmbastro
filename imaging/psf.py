def stack_psf(image,xpos,ypos,sample=1,width=48):
# Create PSF model from stars in an image.
# Takes an image and the star positions and x and y coords.
    import numpy as np
    import scipy.misc as misc
    from scipy import interpolate
    import scipy.optimize as opt
    import astropy.stats as ast
    from scipy.ndimage.filters import gaussian_filter
# Set input object list variables
    nstar    = len(xpos)
    xint     = xpos.astype(np.int)
    yint     = ypos.astype(np.int)

# Define thumbnail and output image parameters
    wdth_th  = np.int(width*1.2)
    wdth_mod = width * sample
    wdth_thr = wdth_th * sample
    starcube = np.zeros((2*wdth_mod+1,2*wdth_mod+1,nstar))
    xold = np.linspace(0,2*wdth_th,num=2*wdth_th+1)*sample
    yold = np.linspace(0,2*wdth_th,num=2*wdth_th+1)*sample
    x,y = np.meshgrid(xold, yold)

# Extract cube of stellar images
    for i in np.arange(nstar):
        thumb = image[xint[i]-wdth_th:xint[i]+wdth_th+1,yint[i]-wdth_th:yint[i]+wdth_th+1]
        thumb = thumb-np.nanmedian(thumb[thumb<np.nanpercentile(thumb,56)])
        # Resample if necessary
        if sample != 1:
            xnew = np.linspace(0,2*wdth_thr,num=2*wdth_thr+1)
            ynew = np.linspace(0,2*wdth_thr,num=2*wdth_thr+1)
            f = interpolate.interp2d(xold, yold, thumb, kind='cubic')
            thumb = f(xnew, ynew)
            x,y = np.meshgrid(xnew, ynew)
        # Perform fit to ascertain centre of profile
        initial_guess = (thumb[wdth_thr,wdth_thr],
            wdth_thr+(xpos[i]-xint[i])*sample,wdth_thr+(ypos[i]-yint[i])*sample,
            sample,sample)
        popt, pcov = opt.curve_fit(Gaus2D, (x,y), thumb.ravel(), p0 = initial_guess)
        # Add the thumbnail to the cube
        thumb_out = thumb[popt[1]-wdth_mod:popt[1]+wdth_mod+1,popt[2]-wdth_mod:popt[2]+wdth_mod+1]
        thumb_out[thumb_out>1.2*popt[0]] = np.nanmedian(thumb_out)
        starcube[:,:,i] = thumb_out/popt[0]
    # Median the cube to give the final PSF model image
    starcube = ast.sigma_clip(starcube,axis=2,sigma_lower=3.6,sigma_upper=1.6,iters=16)
    psfmod = np.nanmedian(starcube,axis=2)
    psfmod = psfmod - np.nanmedian(psfmod[psfmod<np.nanpercentile(psfmod,56)])
    psfmod = psfmod/np.max(psfmod)
    med = np.nanmedian(psfmod[psfmod<np.percentile(psfmod,56)])
    std = np.std(psfmod[psfmod<np.percentile(psfmod,56)])
    psfmod[gaussian_filter(psfmod,1.0)<med] = med
    return psfmod,sample


def Gaus2D((x, y), amplitude, x0, y0, sigma_x, sigma_y):
    import numpy as np
    gaus = amplitude * np.exp(-((x-x0)**2/(2*sigma_x**2) + (y-y0)**2/(2*sigma_y**2)))
    return gaus.ravel()

def fit_psf(image,psfimage,sample=1):
    import numpy as np
    from scipy import interpolate
    import scipy.optimize as opt
# Set up interpolation model as global variable, to make available in psffunc()
# below.
    global psfmod
# Set up array params
    image   = image-np.nanmedian(image[image<np.percentile(image,64)])
    imsz    = np.shape(image)
    imx     = np.arange(-(imsz[0]-1)/2,+(imsz[0]-1)/2+1)
    imy     = np.arange(-(imsz[1]-1)/2,+(imsz[1]-1)/2+1)
# Resample the data if asked to
    if sample != 1:
        f = interpolate.interp2d(imx, imy, image, kind='cubic')
        imx = np.linspace(-(imsz[0]-1)/2,+(imsz[0]-1)/2+1,num=imsz[0]*sample)
        imy = np.linspace(-(imsz[1]-1)/2,+(imsz[1]-1)/2+1,num=imsz[1]*sample)
        image = f(imx, imy)
# Set up the PSF model using interp2d
    psfsz   = np.shape(psfimage)
    psfx    = np.arange(-(psfsz[0]-1)/2,+(psfsz[0]-1)/2+1)
    psfy    = np.arange(-(psfsz[1]-1)/2,+(psfsz[1]-1)/2+1)
    psfmod  = interpolate.interp2d(psfx, psfy, psfimage, kind='cubic')
# Now fit the data
    psfpara = (image[(imsz[0]-1)/2,(imsz[1]-1)/2],0,0)
    psfopt,psfcov = opt.curve_fit(psffunc, (imx,imy), 1.08*image.ravel(), p0 = psfpara)
    print psfopt
    psffit  = psffunc((imx,imy),*psfopt)
    psffit  = psffit
    return psffit

def psffunc((x, y), amplitude, x0, y0):
    global psfmod
    psf = amplitude*psfmod(x-x0,y-y0)
    return psf.ravel()
