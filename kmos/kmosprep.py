def guides_refs(fcen):
# Retrieve external catalogue guide & reference stars.
# !!! Warning: make sure these line up with target astrometry !!!
    import subprocess
    import numpy as np
    try:
        x = subprocess.check_output(["aclient","axel.u-strasbg.fr","1660","ucac4","-c {0} {1}".format(fcen[0],fcen[1]),"-r 28","-m 1000000"])
    except:
        x = subprocess.check_output(["aclient_cgi","130.79.129.161","ucac4","-sr","-c {0} {1}".format(fcen[0],fcen[1]),"-r 28","-m 1000000"])
    ucacfile = open('ucac_temp.asc','w')
    ucacfile.write(x)
    ucacfile.close()
    ucaccat  = np.genfromtxt('ucac_temp.asc',dtype=None,delimiter='|')

    purpose = np.array(ucaccat["f0"])
    purpose[:] = 'N'
    band       = np.array(purpose)
    mag        = np.zeros(len(ucaccat))
    ra         = np.zeros(len(ucaccat))
    dec        = np.zeros(len(ucaccat))
    ucid = np.array(ucaccat["f0"])
    for i in np.arange(len(ucaccat)):
         ra[i]  = ucaccat["f1"][i][0:11]
         dec[i] = ucaccat["f1"][i][11:21]
         jmag = ucaccat["f7"][i][11:16]
         if jmag.replace(" ", "") != '---': jmag = np.float(jmag)
         rmag = ucaccat["f8"][i][33:38]
         if rmag.replace(" ", "") != '---': rmag = np.float(rmag)
         if (jmag <= 19) & (rmag > 13) & (rmag != '  ---'):
             purpose[i] = 'R'
             mag[i]     = jmag
             band[i]    = 'J'
         if (jmag <= 19) & (rmag <= 13) & (rmag > 12) & (rmag != '  ---'):
             purpose[i] = 'RG'
             mag[i]     = jmag
             band[i]    = 'J'
         elif (rmag <= 12) & (rmag != '  ---'):
             purpose[i] = 'G'
             mag[i]     = rmag
             band[i]    = 'R'
    gnr  = np.where(purpose != 'N')[0]
    ucid = ucid[gnr]
    return ucid,ra[gnr],dec[gnr],purpose[gnr],mag[gnr],band[gnr]

def write_kcat(kcat,name,ra,dec,target,mag,band,priority,comment):
    import numpy as np
    try:
        ntarg = len(ra)
        print '{0} objects found'.format(ntarg)
    except:
        ntarg = 1
    try:
        mag_float = mag.astype(np.float)
        for i in np.arange(len(mag)):
            mag[i] = '{0:5.2f}'.format(mag_float[i])
    except:
        print 'Note mag does not float.'
    try:
        priority = priority.astype(np.intc)
        pri_float = priority.astype(np.intc)
        for i in np.arange(len(priority)):
            priority[i] = '{0:1}'.format(pri_float[i])
    except:
        print 'Note pri does not float.'
    if ntarg > 1:
        for i in np.arange(ntarg):
            kcat.write("{0:16} {1:10.6f} {2:10.6f} {3} {4:5} {5:1} {6:1} {7} \n".format(name[i],ra[i],dec[i],target[i],mag[i],band[i],priority[i],comment[i]))
    else:
        kcat.write("{0:16} {1:10.6f} {2:10.6f} {3} {4:5} {5:1} {6:1} {7} \n".format(name,ra,dec,target,mag,band,priority,comment))

    return
