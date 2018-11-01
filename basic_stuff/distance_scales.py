def pltdist(z):
    from astropy.cosmology import WMAP9 as cosmo
    import numpy as np
    import matplotlib.pyplot as plt
    # z      = (np.arange(36)+1.)/10.
    print cosmo.H0
    dczh0  = 3.e5*z/cosmo.H0
    dcom   = cosmo.comoving_distance(z)
    dlum   = dcom*(1.+z)
    dpro   = dcom/(1.+z)
    plt.figure(1)
    plt.plot(z,dczh0,label=r'$d_{czh0}$')
    plt.plot(z,dlum,label=r'$d_{lum}$')
    plt.plot(z,dcom,label=r'$d_{com}$')
    plt.plot(z,dpro,label=r'$d_{pro}$')
    plt.legend()
    plt.show()
