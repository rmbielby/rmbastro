def make_lephare_cat(filename,id,mag,magerr,context=None,zspec=None,string=None):
    import numpy as np
    lepharecat = open(filename,'w')
    nobj   = len(id)
    nfilt  = np.shape(mag)[0]
    for j in np.arange(nobj):
        lepharecat.write('{0:8.0f} '.format(id[j]))
        for i in np.arange(nfilt):
            if (np.abs(mag[i,j]) > 40) | (mag[i,j] < 10):
                mag[i,j] = -99.
                magerr[i,j] = -99.
            lepharecat.write('{0:8.4f} {1:8.4f} '.format(mag[i,j],magerr[i,j]))
        if context != None:
            sel = np.where(np.abs(mag[:,j]-20)<10)[0]
            lepharecat.write('{0:4.0f} {1:8.5f} '.format(context[j],zspec[j]))
            # if zspec[j] > -1:
            #     print zspec[j]
            if string[j] != None:
                lepharecat.write('{0}'.format(string[j]))
        lepharecat.write('\n')
    lepharecat.close()
    return

def make_gazpar_cat(filename,id,mag,magerr,zspec,flag,ra,dec,mask):
    import numpy as np
    lepharecat = open(filename,'w')
    nobj   = len(id)
    nfilt  = np.shape(mag)[0]
    for j in np.arange(nobj):
        lepharecat.write('{0} '.format(id[j]))
        for i in np.arange(nfilt):
            if (np.abs(mag[i,j]) > 40) | (mag[i,j] < 10):
                mag[i,j] = -99.
                magerr[i,j] = -99.
            lepharecat.write('{0:8.4f} {1:8.4f} '.format(mag[i,j],magerr[i,j]))
        lepharecat.write('{0:7.4f} {1:3.0f} {2:12.6f} {3:12.6f} {4:3.0f}'.format(zspec[j],flag[j],ra[j],dec[j],mask[j]))
        lepharecat.write('\n')
    lepharecat.close()
    return

def make_redshift_cat(filename,id,ra,dec,f140mag,frad,sexflag,mag,magerr,photz,grism,specz,bnames=['u','g','r','i','z','F140','F160']):
    import numpy as np
    lepharecat = open(filename,'w')
    nobj   = len(id)
    nfilt  = np.shape(mag)[0]
    lepharecat.write('#{0:8} {1:10} {2:10} {3:8} {4:8} '.format('ID','R.A.','Dec','F140-au','F140-rad'))
    for i in np.arange(nfilt):
        lepharecat.write('{0:8} {1:8} '.format(bnames[i]+'-ap',bnames[i]+'-er'))
    lepharecat.write(' {0:10} {1:10} {2:10} {3:4}'.format('Photz','Photz-min','Photz-max','PMod'))
    lepharecat.write(' {0:5} {1:10} {2:4}'.format('Gr-ID','Gr-z','Gr-q'))
    lepharecat.write(' {0:5} {1:10} {2:4}'.format('Sp-ID','Sp-z','Sp-q'))
    lepharecat.write(' {0:5}'.format('pFLAG'))
    lepharecat.write('\n')
    print np.shape(id),np.shape(photz),np.shape(grism),np.shape(specz)
    for j in np.arange(nobj):
        lepharecat.write('{0:8} {1:10.6f} {2:10.6f} {3:8.4f} {4:8.2f} '.format(id[j],ra[j],dec[j],f140mag[j],frad[j]))
        for i in np.arange(nfilt):
            lepharecat.write('{0:8.4f} {1:8.4f} '.format(mag[i,j],magerr[i,j]))
        lepharecat.write(' {0:10.5f} {1:10.5f} {2:10.5f} {3:4.0f}'.format(photz[0,j],photz[1,j],photz[2,j],photz[3,j]))
        lepharecat.write(' {0:5.0f} {1:10.5f} {2:4.0f}'.format(grism[0,j],grism[1,j],grism[2,j]))
        lepharecat.write(' {0:5.0f} {1:10.5f} {2:4.0f}'.format(specz[0,j],specz[1,j],specz[2,j]))
        lepharecat.write(' {0:5.0f}'.format(sexflag[j]))
        lepharecat.write('\n')
    lepharecat.close()
    return
