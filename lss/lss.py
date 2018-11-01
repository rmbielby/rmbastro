import numpy as np


def densmap(ra,dec,mag,maglim=19.5,cellsize=0.5,gridsize=0.064):
    import numpy as np
    meandec = (np.max(dec)+np.min(dec))/2.
    minra   = np.min(ra)
    maxra   = np.max(ra)
    mindec  = np.min(dec)
    maxdec  = np.max(dec)
    nra     = np.int((np.max(ra) - np.min(ra))*np.cos(meandec*np.pi/180.)/
(gridsize/60.))+1
    ndec    = np.int((np.max(dec) - np.min(dec))/(gridsize/60.))+1
    dengrid = np.zeros((ndec,nra))
    rabin   = np.zeros(nra)
    decbin  = np.zeros(ndec)
    limsel = np.where(mag<maglim)[0]
    for i in np.arange(nra):
        rabin[i] = minra+(i+0.5)*gridsize/60./np.cos(meandec*np.pi/180.)
        for j in np.arange(ndec):
            decbin[j]    = mindec+(j+0.5)*gridsize/60.
            ngal         = len(np.where((((ra[limsel]-rabin[i])*np.cos(meandec*np.pi/180.))**2 +
                (dec[limsel]-decbin[j])**2)**0.5<=cellsize/60.)[0])
            dengrid[j,i] = np.float(ngal)/np.pi/(cellsize)**2.
    return rabin,decbin,dengrid

def radecred_to_xyz(coords):
    import numpy as np
    from astropy.cosmology import Planck15 as cosmo
    # print np.shape(coords)
    ra,dec,redshift = coords[0,:],coords[1,:],coords[2,:]
    # print np.shape(ra)
    ngal = len(ra)
    x = np.zeros(ngal)
    y = np.zeros(ngal)
    z = np.zeros(ngal)
    for i in xrange(len(ra)):
        dist = cosmo.comoving_distance(redshift[i]).value
        x[i] = dist*np.cos(dec[i]*np.pi/180.)*np.cos(ra[i]*np.pi/180.)
        y[i] = dist*np.cos(dec[i]*np.pi/180.)*np.sin(ra[i]*np.pi/180.)
        z[i] = dist*np.sin(dec[i]*np.pi/180.)
    #        print x[i],y[i],z[i]
    #    print '***********************'
    xyz = np.array([x,y,z])
    print np.shape(xyz)
    return xyz

def paircount(xyz1,xyz2,bins=None):
    import numpy as np
    x1,y1,z1 = xyz1[0,:],xyz1[1,:],xyz1[2,:]
    x2,y2,z2 = xyz2[0,:],xyz2[1,:],xyz2[2,:]
    n1    = len(x1)
    n2    = len(x2)
    # print n1,n2
    s2    = np.zeros(n2)
    nbin = len(bins)
    pairs = np.zeros(nbin,dtype=np.uint32)
    sigarr = []
    piarr  = []
    pairs,sigpairs = 0.,0.
    s2 = (x2**2+y2**2+z2**2)**0.5
    s1 = (x1**2+y1**2+z1**2)**0.5
    for i in xrange(n1):
        ang    = (x1[i]*x2+y1[i]*y2+z1[i]*z2)/(s1[i]*s2)
        sel    = np.where(ang <1)[0]
        theta  = np.arccos(ang[sel])
        pi     = abs(s1[i]-s2[sel])
        sigma  = 0.5*theta*(s1[i]+s2[sel])
        sigarr,piarr = np.append(sigarr,sigma),np.append(piarr,pi)
        hdist  = (sigma**2+pi**2)**0.5
        pairs  += np.histogram(hdist,bins=bins)[0]
        sigpairs += np.histogram(sigma,bins=bins)[0]
    if len(sigarr) > 0:
        print 'min(sigma)',np.min(sigarr)
    if len(piarr) > 0:
        print 'min(pi)',np.min(piarr)
    return pairs,sigpairs

def landszal(coords1,coords2,rand1,rand2,type='radecz',bins=None):
    import numpy as np
    if type == 'radecz':
            xyz1  = radecred_to_xyz(coords1)
            xyz2  = radecred_to_xyz(coords2)
            rxyz1 = radecred_to_xyz(rand1)
            rxyz2 = radecred_to_xyz(rand2)
    rfac1 = len(rxyz1[0])/np.float(len(xyz1[0]))
    rfac2 = len(rxyz2[0])/np.float(len(xyz2[0]))
    # print rfac1,rfac2
    # print 'rand1: '
    # for i in np.arange(len(rand1[0,:])): print rand1[:,i]
    # print 'rand2: '
    # for i in np.arange(len(rand2[0,:])): print rand2[:,i]
    # print 'rxyz1: '
    # for i in np.arange(len(rxyz1[0,:])): print rxyz1[:,i]
    # print 'rxyz2: '
    # for i in np.arange(len(rxyz2[0,:])): print rxyz2[:,i]
    dd,ddsig = paircount(xyz1,xyz2,bins=bins)
    rd,rdsig = paircount(rxyz1,xyz2,bins=bins)
    dr,drsig = paircount(xyz1,rxyz2,bins=bins)
    rr,rrsig = paircount(rxyz1,rxyz2,bins=bins)
    wtheta  = (dd-rd/rfac1-dr/rfac2+rr/rfac1/rfac2)/(rr/rfac1/rfac2+1.e-5)
    wtherr  = np.abs((1.+wtheta))/((dd+0.5)**0.5)
    return dd,rd,dr,rr,wtheta,wtherr

def integral_constraint(r,rr,slope=1.8):
    import numpy as np
    return np.sum(rr*r**(-np.abs(slope)))/np.sum(rr)

def plaw_slope(r,*p):
    r0,slope = p
    return (r/r0)**(-np.abs(slope))

def plaw_fixg(r,*p):
    r0    = p
    slope = fixg
    return (r/r0)**(-np.abs(slope))

def fit_plaw(s,xi,xie,r0=3.,gamma=1.8):
    from scipy.optimize import curve_fit
    p0 = [r0,gamma]
    coeff, var_matrix = curve_fit(plaw_slope, s, xi,sigma=xie, p0=p0)
    return coeff[0],var_matrix[0,0]**0.5,coeff[1],var_matrix[1,1]**0.5
