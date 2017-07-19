
def virial_radius(ra,dec,redshift):
# Estimate Virial Radius using relations from Ramella+89
    import numpy as np
    h0 = 73.
    zref   = np.median(redshift)
    zreferr  = np.std(redshift)
    thsum  = 0.
    nmemb  = len(ra)
    rafact = np.cos(np.mean(dec)/180.*np.pi)
    for i in np.arange(nmemb-2):
        for j in np.arange(i+1,nmemb):
            thsum += 1./((((ra[i] - ra[j])*rafact)**2+(dec[i] - dec[j])**2)**0.5/180.*np.pi)
    rh       = np.pi*3.e5*zref/h0*np.sin(0.5*(nmemb*(nmemb-1)/2./thsum))
    rherr    = rh*((zreferr/zref)**2.+2./nmemb)**0.5
    print 'Virial radius = {0}+-{1}'.format(rh,rherr)
    return rh,rherr


def virial_mass(sigv2,sigv_err,rh,rherr):
# Estimate Virial Mass using relations from Ramella+89
    import numpy as np
    mpc2m    = 206265.*1.e6*1.5e11 # Conversion factor from Mpc to metres
    mh       = 6.*sigv2*1.e6*rh*mpc2m/6.67e-11/2.e30
    mherr    = mh*((sigv_err**2/sigv2)**2+(rherr/rh)**2)**0.5
    print 'Virial mass = {0:8.2e}+-{1:8.2e}'.format(mh,mherr)
    return mh,mherr
