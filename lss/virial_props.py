def veldisp(redshift):
    import numpy as np
    zref = np.median(redshift)
    dvel = 3.e5*(redshift-zref)/(1.+0.5*(redshift+zref))
    return 3.**0.5*np.std(dvel)

def virial_radius(ra,dec,redshift):
    import numpy as np
    from astropy.cosmology import WMAP9 as cosmo
    zref   = np.median(redshift)
    zreferr  = np.std(redshift)
    thinvsum  = 0.
    nmemb  = len(ra)
    rafact = np.cos(np.median(dec)/180.*np.pi)
    for i in np.arange(nmemb-2):
        for j in np.arange(i+1,nmemb):
            thinvsum += 1./((((ra[i] - ra[j])*rafact)**2+(dec[i]-dec[j])**2)**0.5/180.*np.pi)
    thinvmean = 2.*thinvsum/((nmemb-2)*nmemb)
    rh       = np.pi*cosmo.comoving_distance(zref)/(1.+zref)/2./thinvmean
    rherr    = rh*((zreferr/zref)**2.+2./nmemb)**0.5
    print 'Virial radius = {0}+-{1}'.format(rh,rherr)
    return rh,rherr


def virrad_ram(ra,dec,redshift):
# Estimate Virial Radius using relations from Ramella+89
    import numpy as np
    from astropy.cosmology import WMAP9 as cosmo
    zref   = np.median(redshift)
    zreferr  = np.std(redshift)
    thsum  = 0.
    nmemb  = len(ra)
    rafact = np.cos(np.median(dec)/180.*np.pi)
    print rafact
    for i in np.arange(nmemb-2):
        for j in np.arange(i+1,nmemb):
            thsum += 1./((((ra[i] - ra[j])*rafact)**2+(dec[i]-dec[j])**2)**0.5/180.*np.pi)
    rh       = np.pi*cosmo.comoving_distance(zref).value*np.sin(0.5*(nmemb*(nmemb-1)/2./thsum))
    rherr    = rh*((zreferr/zref)**2.+2./nmemb)**0.5
    print 1./thsum/nmemb**2/np.pi*180.*60.,'arcmin'
    print 0.5*(nmemb*(nmemb-1)/2./thsum)/np.pi*180.*60.,'arcmin'
    print 'Virial radius = {0}+-{1}'.format(rh,rherr)
    return rh,rherr


def virmass_ram(sigv2,sigv_err,rh,rherr):
# Estimate Virial Mass using relations from Ramella+89
    import numpy as np
    mpc2m    = 206265.*1.e6*1.5e11 # Conversion factor from Mpc to metres
    mh       = 6.*sigv2*1.e6*rh*mpc2m/6.67e-11/2.e30
    mherr    = mh*((sigv_err**2/sigv2)**2+(rherr/rh)**2)**0.5
    print 'Virial mass = {0:8.2e}+-{1:8.2e}'.format(mh,mherr)
    return mh,mherr

def virmass(ra,dec,redshift):
# Estimate Virial Mass using relations from Ramella+89
    import numpy as np
    from astropy.cosmology import WMAP9 as cosmo
    mpc2m    = 206265.*1.e6*1.5e11 # Conversion factor from Mpc to metres
    zref     = np.median(redshift)
    sigvp    = np.std((3.e5*(redshift-zref)/(1.+(redshift+zref)/2.)))
    rsum     = 0.
    rden     = 0
    for i in np.arange(len(ra)-1):
        roff = ((ra[i]-ra[i+1:])**2+(dec[i]-dec[i+1:])**2)**0.5
        rsum += np.nansum(1./roff)
        rden += len(roff)
    rinv     = rsum/cosmo.comoving_distance(zref)
    rden     = len(ra) * (len(ra)-1.)/2.
    rpv      = 1./(rinv.value/rden) # Projected virial radius in Mpc
    print sigvp
    print '0.002*sig_p=',0.002*sigvp
    print 'R_VP = ',rpv
    rv       = 0.002*sigvp #np.pi/2.*rpv
    print 'R_V = ',rv
    mh       = 3.*(sigvp*1.e3)**2/(6.67e-11)*(rv*mpc2m)/2.e30
    # mherr    = mh*((sigv_err**2/sigv2)**2+(rherr/rh)**2)**0.5
    print 'Virial mass = {0:8.2e}+-{1:8.2e}'.format(mh,0.1)
    return mh,rv
