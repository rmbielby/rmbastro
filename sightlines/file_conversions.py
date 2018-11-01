def readqsodata(spf,source='XSHOOTER'):
    import numpy as np
    import astropy.io.fits as fits
    from scipy.ndimage.filters import gaussian_filter1d
    import os.path
    import glob
    if source == 'XSHOOTER':
        print spf
        if os.path.isfile(spf):
            if spf.split('.')[-1] == 'dat':
                wave,flux,flue,cont = np.genfromtxt(spf,unpack=True)
                return wave,flux,flue,cont
            elif spf.split('.')[-1] == 'fits':
                hdulist   = fits.open(spf)
                flux = hdulist[0].data
                flue = hdulist[1].data
                wave = hdulist[2].data
                cont = hdulist[3].data
                return wave,flux,flue,cont
        else:
            return 0
    elif source == 'MIKE':
        import astropy.io.fits as fits
        print spf
        if os.path.isfile(spf):
            hdu = fits.open(spf)
            # for i in np.arange(4): print hdu[i].data
            wave,flux,flue = hdu[2].data,hdu[0].data,hdu[1].data
            flux = gaussian_filter1d(flux,3.)
            flue = gaussian_filter1d(flue,3.)/3.**0.5
            cont = flux*0. + 1.
            return wave,flux,flue,cont
        else:
            return 0
    elif source == 'UVES':
        import astropy.io.fits as fits
        print spf
        if os.path.isfile(spf):
            hdu = fits.open(spf)
            print np.shape(hdu[0].data)
            # for i in np.arange(4): print hdu[i].data
            wl  =  10.**(hdu[0].header['CRVAL1']+(np.arange(1,hdu[0].header['NAXIS1']+1)-hdu[0].header['CRPIX1'])*hdu[0].header['CD1_1'])
            wave,flux,flue,cont = wl,hdu[0].data[0],hdu[0].data[1],hdu[0].data[3]
            flux = gaussian_filter1d(flux,3.)
            flue = gaussian_filter1d(flue,3.)/3.**0.5
            cont = flux*0. + 1.
            return wave,flux,flue,cont
        else:
            return 0
    elif source == 'HIRES':
        import astropy.io.fits as fits
        print spf
        if spf[-10:] == '_coadd.fits':
            hdu = fits.open(spf)
            wave,flux,flue = hdu[2].data,hdu[0].data,hdu[1].data
            flux = gaussian_filter1d(flux,3.)
            flue = gaussian_filter1d(flue,3.)/3.**0.5
            cont = flux*0. + 1.
            return wave,flux,flue,cont
        else:
            hdu = fits.open(spf)
            flux = hdu[0].data
            wave  =  10.**(hdu[0].header['CRVAL1']+(np.arange(1,hdu[0].header['NAXIS1']+1)-hdu[0].header['CRPIX1'])*hdu[0].header['CDELT1'])
            spfe = spf.split('f.fits')[0]+'e.fits'
            hdu = fits.open(spfe)
            flue = hdu[0].data
            cont = flux*0. + 1.
            return wave,flux,flue,cont
    elif source == 'ESI':
        import astropy.io.fits as fits
        print spf
        if os.path.isfile(spf):
            hdu = fits.open(spf)
            wave,flux,flue,cont = hdu[2].data,hdu[0].data,hdu[1].data,hdu[3].data
            print wave,flux,flue,cont
            return wave,flux,flue,cont
        else:
            return 0

def pyigm2alis(quasar,jfile,zqso,spectralist,zsys=3.0,vwin=2400.,resn = 40.,folder='./',dolyseries=6,zll=2.8):
# import rmbastro.sightlines.file_conversions as fconv
# x = fconv.pyigm2alis('IGM_model.json')
    import json
    import numpy as np
# Convert JSON file to python dictionary:


    outfile = '{0}{1}_{2:6.4f}_alis.mod'.format(folder,quasar,zsys)
    alisout = open(outfile,'w')
    # ALIS opening
    alis_intro = [['run ncpus -1'],
        ['run nsubpix 5'],
        ['run blind False'],
        ['run convergence False'],
        ['run convcriteria 0.2'],
        ['chisq atol 0.001'],
        ['chisq xtol 1.e-8'],
        ['chisq ftol 1.e-8'],
        ['chisq gtol 1.e-8'],
        ['chisq  miniter  10'],
        ['chisq  maxiter  3000'],
        ['out model True'],
        ['out fits True'],
        ['out verbose 1'],
        ['out overwrite True'],
        ['out covar datafile.mod.out.covar'],
        ['plot dims 3x3'],
        ['plot ticklabels True'],
        ['plot labels True'],
        ['plot fitregions True'],
        ['plot fits True']]

    for i in np.arange(len(alis_intro)):
        alisout.write('{0}\n'.format(alis_intro[i][0]))

    zwin  = vwin/3.e5*(1.+zsys)
    print 'zwin = ',zwin
    # l0 = (1.+za)*la
    # l0 = (1.+zb)*lb
    # (1.+za)*la = (1.+zb)*lb
    # zb = (1.+za)*la/lb-1.

    lyseries = np.array([1215.67,1025.7,972.5,949.7,937.8,912.])
    if dolyseries < len(lyseries):
        lyseries = lyseries[:dolyseries]
    nly = len(lyseries)

    zlys     = (1.+zsys)/(lyseries)*1215.67-1.
    zbeta = (1.+zsys)/(1025.7)*1215.67-1.
    zgama = (1.+zsys)/972.5*1215.67-1.
    zbalp = (1.+zsys)*1025.7/1215.67-1.
    # print zsys,zbeta,zgama,zbalp,zqso
    # zranges = np.zeros(3,2)
    zarng  = [zsys-zwin,zsys+zwin]
    zbrng  = [zbeta-zwin,zbeta+zwin]
    zgrng  = [zgama-zwin,zgama+zwin]


    # print zarng,zbrng,zgrng
    zranges = np.array([zarng])
    # print np.shape(zranges)
    # if zgrng[0] < zqso:
    #     zranges = np.array([zarng,zbrng,zgrng])

    for i in np.arange(1,len(lyseries)):
        print zranges,[zlys[i]-zwin,zlys[i]+zwin]
        if (zlys[i] < zqso):
            zranges = np.append(zranges,[[zlys[i]-zwin,zlys[i]+zwin]],axis=0)
    for i in np.arange(1,len(lyseries)):
        zbalp  = (1.+zsys)*lyseries[i]/1215.67-1.
        zbarng = [zbalp-zwin,zbalp+zwin]
        # if zbalp > 2.8:
        zranges = np.append(zranges,[zbarng],axis=0)
        print zbalp,zbarng
    # print np.shape(zranges)
    nzrng  = np.shape(zranges)[0]
    # print 'nzrng = ',  nzrng
    # print 'zranges = ',zranges

    # sfile = '{0}/{1}_z{2:6.4f}_alis.dat'.format(folder,quasar,zsys)

    wlims = [3600.,5000.]
    # resn = 40. # resolution in km/s, X-Shooter ~ 40km/s (1" slit)

    data_file = open(jfile)
    jdata      = json.load(data_file)
# Get the system IDs (zRedshift_Ion)
    systems   = jdata['cmps'].keys()

    civz,nciv   = [],0
    siivz,nsiiv = [],0
    siiiz,nsiii = [],0
    feiiz,nfeii = [],0
    mgiiz,nmgii = [],0
    for system in systems:
        ion = system.split('_')[-1]
        print ion,system
        if ion == 'CIV':
            if np.abs(jdata['cmps'][system]['zfit']-zsys)<zwin:
                civz.append(jdata['cmps'][system]['zfit'])
                nciv += 1
        if ion == 'SiIV':
            if np.abs(jdata['cmps'][system]['zfit']-zsys)<zwin:
                siivz.append(jdata['cmps'][system]['zfit'])
                nsiiv += 1
        if ion == 'SiII':
            if np.abs(jdata['cmps'][system]['zfit']-zsys)<zwin:
                siiiz.append(jdata['cmps'][system]['zfit'])
                nsiii += 1
        if ion == 'FeII':
            if np.abs(jdata['cmps'][system]['zfit']-zsys)<zwin:
                feiiz.append(jdata['cmps'][system]['zfit'])
                nfeii += 1
        if ion == 'MgII':
            if np.abs(jdata['cmps'][system]['zfit']-zsys)<zwin:
                mgiiz.append(jdata['cmps'][system]['zfit'])
                nmgii += 1

    print 'C IV redshifts = ',civz
    print 'SiIV redshifts = ',siivz
    print 'SiII redshifts = ',siiiz
    print 'MgII redshifts = ',mgiiz
    print 'FeII redshifts = ',feiiz
    civz = np.array(civz)
    siivz = np.array(siivz)
    siiiz = np.array(siiiz)
    mgiiz = np.array(mgiiz)
    feiiz = np.array(feiiz)


    nspec = nly*nzrng
    alisout.write('data read\n')

# Create ascii files of spectra suitable for ALIS.
    x = np.genfromtxt(spectralist,unpack=True,dtype=None)
    source_sfiles,instruments,resolution = x['f0'],x['f1'],x['f2']
    print source_sfiles,instruments,resolution
    print np.shape(source_sfiles),
    if np.shape(np.transpose([source_sfiles]))[0] == 1:
        source_sfiles = np.array([source_sfiles])
        instruments   = np.array([instruments])
        resolution    = np.array([resolution])
    sfiles = []
    for i in np.arange(len(source_sfiles)):
        data = readqsodata(source_sfiles[i],source=instruments[i])
        sfiles.append('{0}_{1}_alisspec.dat'.format(quasar,instruments[i]))
        good = np.where(data[2] > 0.)[0]
        alisdf = open(sfiles[i],'w')
        for i in np.arange(len(good)):
            alisdf.write('{0} {1} {2}\n'.format(data[0][good[i]],data[1][good[i]]/data[3][good[i]],data[2][good[i]]/data[3][good[i]]))
        alisdf.close()



    wlims       = np.array([[]])
    wlims_sfile = []
    wlims_resn  = []
    swlmax = []
    for specind,sfile in enumerate(sfiles):
        print '***********',sfile
        wl,flux,flue = np.genfromtxt(sfile,unpack=True,usecols=(0,1,2))
        swlmax.append(np.max(wl))
        wlcen       = []
        for i in np.arange(nzrng):
            for j in np.arange(nly):
                # print np.shape(zranges),i
                wlim  = lyseries[j]*(1.+zranges[i,:])
                flmean = np.nanmean(flux[(wl>wlim[0]) & (wl<wlim[1])])
                femed  = np.nanmedian(flue[(wl>wlim[0]) & (wl<wlim[1])])
                # print 'flmean=',flmean
                print lyseries[j],wlim,flmean,femed
                if wlcen != []:
                    diff = np.min(np.abs(np.array(wlcen)-((wlim[1]+wlim[0])/2.)))
                else:
                    diff = 10000.
                if (flmean>-0.1) & (femed<0.40) & (diff > (wlim[1]-wlim[0])/10.) & (wlim[1] > (zll+1.)*912.):
                    print 'appending',sfile
                    wlims = np.append(wlims,wlim)
                    wlims_sfile.append(sfile)
                    wlims_resn.append(resolution[specind])
                    wlcen.append((wlim[1]+wlim[0])/2.)
    wlims = np.reshape(wlims,(len(wlims)/2,2))
    wlcen = np.mean(wlims,axis=1)
    print np.shape(wlims_resn)
    # wlcen_int,unique = np.unique(wlcen.astype(np.int),return_index=True)
    # wlcen = wlcen[unique]
    # wlims = wlims[unique,:]
    wlwid = wlims[:,1]-wlims[:,0]
    # print wlims,wlcen
    nspec = len(wlcen)
    for i in np.arange(nspec):
        alisout.write('  {0} specid={1} fitrange=[{2:8.3f},{3:8.3f}] resolution=vfwhm({4}vh) columns=[wave:0,flux:1,error:2] plotone=True \n'.format(wlims_sfile[i], i, wlims[i,0], wlims[i,1], wlims_resn[i]))


    metlims = np.array([[]])
    specnum = nspec-1
    for specind,sfile in enumerate(sfiles):
        if (nciv > 0) & (1551.*(1.+zranges[0,1]) < swlmax[specind]):
            specnum += 1
            wlim  = [1548.*(1.+zranges[0,0]),1551.*(1.+zranges[0,1])]
            metlims = np.append(metlims,wlim)
            alisout.write('  {0} specid={1} fitrange=[{2:8.3f},{3:8.3f}] resolution=vfwhm({4}vh) columns=[wave:0,flux:1,error:2] plotone=True\n'.format(sfile,
                specnum,wlim[0],wlim[1],resolution[specind]))
        if (nsiiv > 0) & (1403.*(1.+zranges[0,1]) < swlmax[specind]):
            specnum += 1
            wlim  = [1393.*(1.+zranges[0,0]),1403.*(1.+zranges[0,1])]
            metlims = np.append(metlims,wlim)
            alisout.write('  {0} specid={1} fitrange=[{2:8.3f},{3:8.3f}] resolution=vfwhm({4}vh) columns=[wave:0,flux:1,error:2] plotone=True\n'.format(sfile,
                specnum,wlim[0],wlim[1],resolution[specind]))
        if (nmgii > 0)& (2804.*(1.+zranges[0,1]) < swlmax[specind]):
            specnum += 1
            wlim  = [2796.*(1.+zranges[0,0]),2804.*(1.+zranges[0,1])]
            metlims = np.append(metlims,wlim)
            alisout.write('  {0} specid={1} fitrange=[{2:8.3f},{3:8.3f}] resolution=vfwhm({4}vh) columns=[wave:0,flux:1,error:2] plotone=True\n'.format(sfile,
                specnum,wlim[0],wlim[1],resolution[specind]))
        if nfeii > 0:
            fewls    = [2344.,2382.,2600.]
            for fewl in fewls:
                specnum += 1
                wlim     = fewl*(1.+zranges[0,:])
                metlims  = np.append(metlims,wlim)
                alisout.write('  {0} specid={1} fitrange=[{2:8.3f},{3:8.3f}] resolution=vfwhm({4}vh) columns=[wave:0,flux:1,error:2] plotone=True\n'.format(sfile,
                    specnum,wlim[0],wlim[1],resolution[specind]))
    if specnum > nspec-1:
        metlims = np.reshape(metlims,(len(metlims)/2,2))
        # print np.shape(wlims),np.shape(metlims)
        wlims   = np.append(wlims,metlims,axis=0)
        # print np.shape(wlims)
    alisout.write('data end\n')
    nspec = specnum

    alisout.write('model read\n   fix vfwhm value True\n   lim voigt bturb [0.2,None]\n   lim voigt ColDens [8.0,22.0]\n')
    alisout.write(' emission\n')
    for i in np.arange(nspec):
        alisout.write('  constant   1.0000CONT  specid={0}\n'.format(i))
    alisout.write('  constant   1.0000CONT  specid={0}\n'.format(nspec))


    alisout.write(' absorption\n')
    for system in systems:
        ion = system.split('_')[-1]
        print ion,system
        if ion == 'HI':
            ion = '1H_I'
            coldens  = jdata['cmps'][system]['Nfit']
            bpara    = jdata['cmps'][system]['bfit']['value']
            redshift = jdata['cmps'][system]['zfit']
            # print redshift,zsys,np.abs(redshift-zsys),zwin
            # if (np.abs(redshift-zsys)<=zwin) | (np.abs(redshift-zbalp)<=zwin) | (np.abs(redshift-zbeta)<=zwin) | (np.abs(redshift-zgama)<=zwin):
            print '*** Accepted ***'
            sel = np.array([])
            for lylam in lyseries:
                slam = lylam*(1.+redshift)
                # print lylam,redshift
                lsel = np.where(np.abs(slam-wlcen)<=wlwid)[0]
                sel  = np.append(sel,lsel)
                # print slam,wlcen,wlwid
            sel = np.unique(sel)
            print redshift,sel,len(sel)
            # try:
            if len(sel) >= 1:
                specid_str = 'specid={0:1.0f}'.format(sel[0])
                print specid_str
                for i in np.arange(1,len(sel)):
                    specid_str = '{0},{1:1.0f}'.format(specid_str,sel[i])
                print ion,coldens,redshift,bpara,specid_str
                print '  voigt  ion={0}  {1:10.5f} {2:10.7f} {3:10.5f}  0.0000ZEROT  {4}\n'.format(ion,coldens,redshift,bpara,specid_str)
                alisout.write('  voigt  ion={0}  {1:10.5f} {2:10.7f} {3:10.5f}  0.0000ZEROT  {4}\n'.format(ion,coldens,redshift,bpara,specid_str))
            # except:
            #     print '!!!!! Skipped object !!!!!!!!!!'
        elif ion == 'CIV':
            if np.abs(jdata['cmps'][system]['zfit']-zsys)<zwin:
                ion = '12C_IV'
                coldens  = jdata['cmps'][system]['Nfit']
                bpara    = jdata['cmps'][system]['bfit']
                redshift = jdata['cmps'][system]['zfit']
                metwl = 1550.*(1.+redshift)
                print wlims[:,0],wlims[:,1]
                print metwl
                sel = np.where((metwl>wlims[:,0]) & (metwl<wlims[:,1]))[0]
                try:
                    print '*** Including {0} ***'.format(ion)
                    specid_str = 'specid={0:1.0f}'.format(sel[0])
                    alisout.write('  voigt  ion={0}  {1:10.5f} {2:10.7f} {3:10.5f}  0.0000ZEROT  {4}\n'.format(ion,coldens,redshift,bpara,specid_str))
                except:
                    print 'Failed to include ion {0}'.format(ion)
        elif ion == 'SiIV':
            if np.abs(jdata['cmps'][system]['zfit']-zsys)<zwin:
                ion = '28Si_IV'
                coldens  = jdata['cmps'][system]['Nfit']
                bpara    = jdata['cmps'][system]['bfit']
                redshift = jdata['cmps'][system]['zfit']
                metwl = 1397.*(1.+redshift)
                print wlims[:,0],wlims[:,1]
                print metwl
                sel = np.where((metwl>wlims[:,0]) & (metwl<wlims[:,1]))[0]
                try:
                    print '*** Including {0} ***'.format(ion)
                    specid_str = 'specid={0:1.0f}'.format(sel[0])
                    alisout.write('  voigt  ion={0}  {1:10.5f} {2:10.7f} {3:10.5f}  0.0000ZEROT  {4}\n'.format(ion,coldens,redshift,bpara,specid_str))
                except:
                    print 'Failed to include ion {0}'.format(ion)
        elif ion == 'MgII':
            if np.abs(jdata['cmps'][system]['zfit']-zsys)<zwin:
                ion = '24Mg_II'
                coldens  = jdata['cmps'][system]['Nfit']
                bpara    = jdata['cmps'][system]['bfit']
                redshift = jdata['cmps'][system]['zfit']
                metwl = 2798.*(1.+redshift)
                print wlims[:,0],wlims[:,1]
                print metwl
                sel = np.where((metwl>wlims[:,0]) & (metwl<wlims[:,1]))[0]
                try:
                    print '*** Including {0} ***'.format(ion)
                    specid_str = 'specid={0:1.0f}'.format(sel[0])
                    alisout.write('  voigt  ion={0}  {1:10.5f} {2:10.7f} {3:10.5f}  0.0000ZEROT  {4}\n'.format(ion,coldens,redshift,bpara,specid_str))
                except:
                    print 'Failed to include ion {0}'.format(ion)
        elif ion == 'FeII':
            if np.abs(jdata['cmps'][system]['zfit']-zsys)<zwin:
                ion = '56Fe_II'
                coldens  = jdata['cmps'][system]['Nfit']
                bpara    = jdata['cmps'][system]['bfit']
                redshift = jdata['cmps'][system]['zfit']
                for fewl in fewls:
                    metwl = fewl*(1.+redshift)
                    print metwl
                    # print lylam,redshift
                    msel = np.where((metwl>wlims[:,0]) & (metwl<wlims[:,1]))[0]
                    sel  = np.append(sel,msel)
                    # print slam,wlcen,wlwid
                sel = np.unique(sel)
                try:
                    print '*** Including {0} ***'.format(ion)
                    specid_str = 'specid={0:1.0f}'.format(sel[0])
                    if len(sel) > 1:
                        for i in np.arange(1,len(sel)):
                            specid_str = '{0},{1:1.0f}'.format(specid_str,sel[i])
                    alisout.write('  voigt  ion={0}  {1:10.5f} {2:10.7f} {3:10.5f}  0.0000ZEROT  {4}\n'.format(ion,coldens,redshift,bpara,specid_str))
                except:
                    print 'Failed to include ion {0}'.format(ion)

    alisout.write('model end\n')
    alisout.close()
    return jdata
