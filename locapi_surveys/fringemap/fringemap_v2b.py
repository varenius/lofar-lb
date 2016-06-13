# ---------------------------------------------------------------------
#        fringemap ( see instructions at end )  v.2 Neal Jackson 150831
#            version history: v1 Neal Jackson 150630
#                             v2 NJ cleaned/updated, bug in FQ fixed
#                             v2a NJ reversed data and process loops (faster)
#   Some instructions at the bottom of this script
#   Note the dependencies below: requires in particular numpy,scipy,
#   pyfits,astLib and all parseltongue libraries including Wizardry
#   If astLib not present, will still run but won't plot Wenss crosses
#   Extra files required: correlate.py, wenss2000.npy (if you want 
#   Wenss crosses), lofipi_aips.py
#
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList, AIPSMessageLog
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
import numpy as np; from numpy import fft; from scipy import ndimage
import re,sys,pickle,os,glob,time,warnings,pyfits,matplotlib,time
from matplotlib import pyplot as plt; from scipy.ndimage import measurements
from pyfits import getdata,getheader
from correlate import *
from lofipi_aips import *
try:
    import astLib; from astLib import astCoords
    have_astLib = True
except:
    print 'No astLib present - running stubbed version (without WENSS plots)'
    have_astLib = False
warnings.simplefilter('ignore')
twopi, logfiledir = 2.0*np.pi, './'
# --------------------------------------------------------------------

def uvw2dh (u,v,w,dec):
    L = np.sqrt(u**2+v**2+w**2)
    sindec,cosdec = np.sin(np.deg2rad(dec)),np.cos(np.deg2rad(dec))
    sind = (w*sindec+v*cosdec)/L
    cosd = np.sqrt(1-sind**2)
    sinHh = (u/L)/cosd
    cosHh = np.sqrt(1-sinHh**2)
    vtry = L*(sind*cosdec-cosd*sindec*cosHh)
    wtry = L*(sind*sindec+cosd*cosdec*cosHh)
    verr,werr = abs(vtry/v-1), abs(wtry/w-1)
    if verr>0.0001 or werr>0.0001:   # baseline direction ambiguity
        sind*=1              # yes, really!
        cosd = -np.sqrt(1-sind**2)
        sinHh = (u/L)/cosd
        cosHh = np.sqrt(1-sinHh**2)
        vtry = L*(sind*cosdec-cosd*sindec*cosHh)
        wtry = L*(sind*sindec+cosd*cosdec*cosHh)
        verr,werr = abs(vtry/v-1), abs(wtry/w-1)
    if verr>0.0001 or werr>0.0001:    # should never get here
        print 'FAILED IN uvw2dh: -->',verr,werr
    return sindec,cosdec,sind,cosd,sinHh,cosHh

def frd2xy (uvw,dec,incr,nchan,ttotal): # 2nd derivs of IPhEq wrt D,H,t,f
    L = np.sqrt((uvw**2).sum())
    sindec,cosdec,sind,cosd,sinHh,cosHh = uvw2dh(uvw[0],uvw[1],uvw[2],dec)
    l1 = twopi**2*L
    l2 = 4.485E-8*L                   # this is something/c  :)
    alpha = l1*ttotal*cosd*cosHh/twopi
    beta =  l1*ttotal*cosd*sindec*sinHh/twopi
    gamma = l2*incr*nchan*cosd*sinHh/twopi
    delta = l2*incr*nchan*(sind*cosdec-cosd*sindec*cosHh)/twopi
    return alpha,beta,gamma,delta     # matrix relating phase rates to H,D

def visstart (wdata,utstart,direction):
    tables = np.asarray(wdata.tables,dtype='S')
    nxtables = tables[tables[:,1]=='AIPS NX']
    startvis = 0
    if not len (nxtables):
        print 'No NX tables: if very slow, run INDXR with cparm 0,5,0'
    else:    # read the most recent NX table
        n_nx = int(max(np.asarray(nxtables[:,0],dtype='int')))
        nx = wdata.table('NX',n_nx)
        for i in range(len(nx)):
            if nx[i]['time']>utstart:
                break
            startvis = nx[i]['start_vis']+direction
    return startvis            

def fring_filter(a):
    ny,nx = a.shape
    iy,ix = measurements.maximum_position(a)
    for i in range(0,max(0,ix-3)) + range(min(nx-1,ix+3),nx-1):
        a[iy,i] = np.median(a[:,i])
    for i in range(0,max(0,iy-3)) + range(min(ny-1,iy+3),ny-1):
        a[i,ix] = np.median(a[i,:])
    return a

def writeit(fmap,ra,dec,fsiz,maxoff,outname):
    h = pyfits.PrimaryHDU(fmap)
    h.header.append(('CTYPE1','RA---SIN'))
    h.header.append(('CRVAL1',ra))
    h.header.append(('CRPIX1',0.5*fsiz))
    h.header.append(('CDELT1',-2.*maxoff/fsiz))
    h.header.append(('CTYPE2','DEC--SIN'))
    h.header.append(('CRVAL2',dec))
    h.header.append(('CRPIX2',0.5*fsiz))
    h.header.append(('CDELT2',2.*maxoff/fsiz))
    h.header.append(('EQUINOX',2000))
    pyfits.writeto(outname,h.data,h.header,clobber=True)

def addwenss (aipsname,indisk,ra,dec,maxoff,pmin,pmax,w_max=1000.,incl='IMAP'):
    if have_astLib:
        wenss = np.load('wenss2000.npy')
        wenss = wenss[wenss[:,2]<2]   # only single or centre of multiple
        wenss = wenss[wenss[:,3]>w_max]  # only sources >w_max mJy
        a = correlate (np.array([[ra,dec]]),0,1,wenss,0,1,maxoff)
        fo = open('starsfile','w')
        for i in a:
            idx = int(i[1])
            sra = astCoords.decimal2hms(wenss[idx,0],' ')
            sdec = astCoords.decimal2dms(wenss[idx,1],' ')
            fo.write('%s %s\n' % (sra,sdec))
        fo.close()
        stars (aipsname,incl,indisk,intext='./starsfile',logfiledir=logfiledir)
    greys (aipsname,incl,indisk,pmin,pmax,5,1,logfiledir=logfiledir)
    lwpla (aipsname,incl,indisk,'./'+aipsname+'_plot.ps',logfiledir=logfiledir)
    AIPSImage(aipsname,incl,indisk,1).zap()
    os.system('rm starsfile')

def shift_func(c,px,py,inshape):
    return (py[c[1],c[0]]+0.5*inshape[0],px[c[1],c[0]]+0.5*inshape[1])

def frget(d,h,i,uvw,times,maxoff,fsiz,dofilt):
    tdec = h['crval'][h['ctype'].index('DEC')]
    nc = h['naxis'][h['ctype'].index('FREQ')]*h['naxis'][h['ctype'].index('IF')]
    chwid = h['cdelt'][h['ctype'].index('FREQ')]
    fring = abs(fft.fftshift(fft.fft2(d)))
    fring = fring_filter(fring) if dofilt else fring
    alpha,beta,gamma,delta = frd2xy (uvw[i],tdec,chwid,nc,times[-1]-times[0])
    xoff,yoff = np.meshgrid(np.arange(fring.shape[1],dtype='float'),\
                            np.arange(fring.shape[0],dtype='float'))
    xoff -= 0.5*fring.shape[1]; yoff -= 0.5*fring.shape[0]
    det = 1.0/(alpha*delta-beta*gamma)
    xp,yp = det*(delta*xoff-beta*yoff), det*(-gamma*xoff+alpha*yoff)
    maxoff = np.deg2rad(maxoff)
    if maxoff==0.0:
        maxoff = max(xp.max(),-xp.min(),yp.max(),-yp.min())
    mat = np.asarray(np.meshgrid(np.arange(float(fsiz)),np.arange(float(fsiz))))
    mat = maxoff*(mat-0.5*fsiz)/(0.5*fsiz)
    px,py = alpha*mat[0]+beta*mat[1],gamma*mat[0]+delta*mat[1]
    fmap = ndimage.geometric_transform(fring,shift_func,\
      output_shape=(fsiz,fsiz),extra_arguments=(px,py,fring.shape)).T
    fmap = np.fliplr (fmap)    # HA -> RA, flip about centre
    return fmap,np.rad2deg(maxoff)

def get_outname (frifiledir,aipsname,pol,antenna):
    i=0
    ps = 'L' if pol else 'R'
    while True:
        test = frifiledir+'FR_'+aipsname+ps+'_'+antenna+'_'+str(i)+'.fits'
        if not os.path.isfile(test):
            return test
        i+=1

def procants (name, data, docomp=5):
    an = data.table('AN',0)
    nums = np.zeros(len(name),dtype='int')
    for anr in an:
        for iname in range(len(name)):
            ncomp = min (docomp, len(name[iname]))
            if name[iname][:ncomp] == anr.anname[:ncomp]:
                nums[iname] = int(anr.nosta)
    status = any(nums==0)
    if status:
        print 'Error: following antennas not in data:'
        print np.asarray(name,dtype='S')[nums==0.0]
    return nums, status

def fringemap(aipsno, aipsname, intl_name, refant, imaxoff=0.0, \
       aipsclass= 'FITS', utstart=0, utstop=0, utinc=0, fsiz=256, \
       dofilt = False, indisk=1, frifiledir='./', w_max=1000., zap=True):
    AIPS.userno = aipsno
    if zap:
        print 'Removing all files beginning FR_ from directory '+frifiledir
        os.system('rm '+frifiledir+'FR_*')
    wdata = WizAIPSUVData (aipsname,aipsclass, indisk, 1)
    intl_antennas,status = procants (intl_name, AIPSUVData(wdata))
    refant,status = procants ([refant], AIPSUVData(wdata))
    if status:
         return
    h = wdata.header
    ra = h['crval'][h['ctype'].index('RA')]
    dec = h['crval'][h['ctype'].index('DEC')]
    npol = min(h['naxis'][h['ctype'].index('STOKES')],2)  # only ll,rr
    nc = h['naxis'][h['ctype'].index('FREQ')]*h['naxis'][h['ctype'].index('IF')]
    f0 = h['crval'][h['ctype'].index('FREQ')]
    chw = h['cdelt'][h['ctype'].index('FREQ')]
    uvfac = 1.0 + 0.5*chw*(nc-1)/f0  # converts uvw to centre of band
    baseline = np.dstack((intl_antennas,list(refant)*len(intl_antennas)))[0]
    baseline = np.array([baseline]) if baseline.ndim==1 else baseline
    for i in range(len(baseline)):
        baseline[i] = np.sort(baseline[i])
    utstart = wdata[0].time-0.0001 if utstart==0.0 else utstart
    if utstop==0.0:
        try:          # quicker if NX table exists
            nx = AIPSUVData(wdata).table('NX',0)
            utstop = nx[len(nx)-1].time+0.5*nx[len(nx)-1].time_interval+0.0001
        except:       # find time of last visibility, can be slow
            print 'No NX table, so checking for end time manually....'
            utstop = wdata[len(wdata)-1].time
    if utinc==0.0:    # doing the whole file
        Utstart, Utstop = [utstart], [utstop]
    else:             # setup array of start and stop times for each chunk
        Utstart = np.arange(utstart,utstop,utinc)
        Utstop = Utstart + utinc
        if Utstop[-1]>utstop: # avoid truncated one at the end
            Utstart, Utstop = Utstart[:-1], Utstop[:-1]
    print 'Time range: %.4f %.4f, intervals %d'%(utstart,utstop,len(Utstart))
#    v0,v1 = visstart (wdata, u0, -1), visstart (wdata, u1, +1)
# ------------------
    ui=-1; amap=np.array([])
    for visibility in wdata:
        vtime = visibility.time
        if vtime < Utstart[0]:
            continue
        elif vtime >= Utstart[0] and vtime < Utstop[-1]:
            for j in range(len(Utstart)):
                if Utstart[j]<=vtime<Utstop[j]:
                    this_ui = j
            if this_ui!=ui:
                ui = this_ui
                if ui:
                    amap,maxoff = proc_chunk (d,u,times,baseline,Utstart[ui],\
                                Utstop[ui],uvfac,npol,intl_antennas,h,imaxoff,\
                                fsiz,dofilt,frifiledir,aipsname,intl_name,\
                                logfiledir,ra,dec,amap)
                d = [np.array([])]*len(baseline)
                u = [np.array([])]*len(baseline)
                times = [np.array([])]*len(baseline)
            else:
                continue_chunk (visibility,baseline,nc,npol,d,u,times)
#        elif vtime > Utstop[-1]:
#            amap,maxoff = proc_chunk (d,u,times,baseline,Utstart[-1],Utstop[-1],uvfac,\
#                          npol,intl_antennas,h,imaxoff,fsiz,dofilt,\
#                          frifiledir,aipsname,intl_name,logfiledir,ra,dec,amap)
#            break
        else:
            print 'WTF?'

    amap,maxoff = proc_chunk (d,u,times,baseline,Utstart[-1],Utstop[-1],uvfac,\
                  npol,intl_antennas,h,imaxoff,fsiz,dofilt,\
                  frifiledir,aipsname,intl_name,logfiledir,ra,dec,amap)
    np.save('amap',amap)
    if imaxoff != 0.0:
        amap = np.rollaxis (amap,2,0)
        acomb = imcomb (amap)
        writeit (acomb,ra,dec,fsiz,maxoff,'fring.fits')
        pload ('./fring.fits',aipsname,indisk,'IMAP',logfiledir=logfiledir,doindxr=False)
        b=np.ravel(acomb);b=np.sort(b[~np.isnan(b)])
        bmin,bmax = b[len(b)/10],b[99*len(b)/100]
        addwenss (aipsname,indisk,ra,dec,maxoff,bmin,bmax,w_max=w_max)

    for name in intl_name:
        nlist,ncube = np.sort(glob.glob('FR_*'+name+'_*.fits')), np.array([])
        for nfile in nlist:
            try:
                ncube = np.dstack((ncube, getdata(nfile)))
            except:
                ncube = np.copy(getdata(nfile))
        h = getheader(nfile)
        ncube = np.rollaxis (ncube,2,0)
        ncomb = imcomb (ncube)
        thisname = './FR_'+aipsname+'_'+name+'.fits'
        writeit (ncomb,h['CRVAL1'],h['CRVAL2'],h['NAXIS1'],\
                 h['NAXIS1']*h['CDELT2'],thisname)
        pload (thisname,aipsname[:6]+'_'+name[:5],indisk,'IMAP',\
               logfiledir=logfiledir,doindxr=False)
        b=np.ravel(ncomb);b=np.sort(b[~np.isnan(b)])
        bmin,bmax = b[len(b)/10],b[99*len(b)/100]
        addwenss (aipsname[:6]+'_'+name[:5],indisk,ra,dec,maxoff,bmin,bmax,w_max=w_max)

# shape is (nif,nc,npol,AphiW)

def continue_chunk (visibility,baseline,nc,npol,d,u,times):
    idx=-1      # ugly test, change for something more pythonic
    for j in range(len(baseline)):
        if all(baseline[j]==visibility.baseline):
            idx=j
    if idx==-1:
        return
    v = visibility.visibility
    dcol = np.zeros((nc,npol),dtype='complex')
    for spw in range(v.shape[0]):
        for i in range(v.shape[1]):
            fq = spw*v.shape[1]+i
            for j in range(npol):
                dcol[fq,j] = complex(v[spw,i,j,0],v[spw,i,j,1])
    try:
        d[idx] = np.dstack((d[idx],dcol))
        u[idx] = np.column_stack((u[idx],visibility.uvw))
    except:
        d[idx] = np.copy(dcol)
        u[idx] = np.copy(visibility.uvw)
    times[idx] = np.append(times[idx],visibility.time)
#    print 'leaving continue_chunk',len(d),d[0].shape

def proc_chunk (d,u,times,baseline,u0,u1,uvfac,npol,intl_antennas,h,imaxoff,fsiz,\
                dofilt,frifiledir,aipsname,intl_name,logfiledir,ra,dec,amap):
    for i in range(len(d)):
        try:
            d[i]=np.rollaxis(d[i],1,0)
        except:
            pass
    uvw = np.zeros((len(baseline),3))
    print '%s t=%.4f-%.4f:uvw*%.4f to band centre' % \
                 (time.ctime().split()[3],u0,u1,uvfac)
    for i in range (len(baseline)):     # apply correction to band centre
        uvw[i] = uvfac * np.array([np.mean(u[i][0]),\
                               np.mean(u[i][1]),np.mean(u[i][2])])
    for pol in range(min(2,npol)):
        for i in range(len(intl_antennas)):
            fmap,maxoff = frget(d[i][pol],h,i,uvw,times[i],imaxoff,fsiz,dofilt)
            outname = get_outname (frifiledir,aipsname,pol,intl_name[i])
            fo = open(logfiledir+aipsname+'.log','a')
            fo.write('**** Processing %s pol:%d ant:%s FOV:%.f\n' % \
                     (aipsname, pol, intl_name[i], maxoff))
            fo.close ()
            writeit (fmap,ra,dec,fsiz,maxoff,outname)
            if imaxoff != 0.0:
                try:
                    amap = np.dstack ((amap,fmap))
                except:
                    amap = np.copy(fmap)
    print '%s Processed chunk, reading data' % time.ctime().split()[3]
    return amap,maxoff

        
def fitsplot (aipsname, antenna, indisk=1, w_max=1000):
    a = np.sort(glob.glob(frifiledir+'FR_'+aipsname+'*'+antenna+'*.fits'))
    h = getheader(a[0])
    for i in a:
        data = getdata(i)
        try:
            acube = np.dstack((acube,data))
        except:
            acube = np.copy(data)
    acube = np.rollaxis(acube,2,0)
    amap = imcomb (acube)
    pyfits.writeto('./fring.fits',amap,h,clobber=True)
    pload('./fring.fits',aipsname,indisk,antenna,logfiledir=logfiledir,doindxr=False)
    b=np.ravel(amap);b=np.sort(b[~np.isnan(b)])
    bmin,bmax = b[len(b)/10],b[99*len(b)/100]
    addwenss (aipsname, indisk, h['crval1'], h['crval2'], \
              0.5*np.deg2rad(h['cdelt2']*h['naxis2']), bmin,bmax,incl=antenna)

def imcomb (acube):
    amean, astd, acen = np.array([]),np.array([]),np.array([])
    layers,ay,ax = acube.shape
    for a in acube:                     # for each layer in the cube
        b=np.sort(np.ravel(a[a!=0.0]))
        b=b[b<2*np.median(b)-b.min()]   # get histg of the main distrn
        a[a!=0.0]-=np.median(b)         # take off the background
        cen = a[ay/2-10:ay/2+10,ax/2-10:ay/2+10].sum()
        a/= cen
        astd=np.append(astd,np.std(b)/cen)
    aweight = np.copy(acube)
    wtot, atot = np.zeros_like(acube[0]),np.zeros_like(acube[0])
    for i in range(len(aweight)):
        np.putmask(aweight[i],aweight[i]!=0.0,1./astd[i]**2)
        atot += acube[i]*aweight[i]
        wtot += aweight[i]
    return atot/wtot        
 
# ----------------------------------------------------------------

#  Instructions (sample runs given below).
#  Arguments to fringemap:
#  1) Aips number
#  2) AIPS INNA of uv-data file
#  3) Array of names of the antennas for which maps are wanted
#  4) Name of the reference antenna
#      -------------  optional arguments ------------
#  5) imaxoff: Half-width of maps (deg) (default: 0=scale separately for each)
#  6) aipsclass: AIPS INCL of uv-data file (default FITS)
#  7) utstart: UT start (default 0 = beginning of data)
#  8) utend: UT end (default 0 = end of data)
#  9) utinc: UT chunk size for adding chunks (default 0 = only one chunk)
#  10) fsiz: Number of pixels of fringe-maps (default 256)
#  11) dofilt: Clean up fringemaps for diffraction spikes? (default no)
#  12) indisk: AIPS INDI of uv data file (default 1)
#  13) frifiledir: Disk to dump plots, nb trailing slash (default current)
#  14) w_max: Maximum WENSS flux to overplot on map (in mJy), default 1000
#  15) zap: remove all files in directory beginning FR_ first, default yes
#
# Outputs:
#  - Set of fits files with fringe maps for each polarization (max 2)
#             on baselines between requested antennas and ref ant
#  - (If argument 10 given) postscript file of combined fringe map, with
#             WENSS positions overplotted

#fringemap(341,'L255231',['DE601','DE602','DE603','DE605','FR606','SE607'],'ST001',aipsclass='FITS',fsiz=256,imaxoff=3.0,dofilt=True,indisk=1,frifiledir='./',w_max=500,utstart=0,utstop=0,utinc=0)

#fringemap(341,'3C266',['DE601','DE602','DE603','DE605','FR606','SE607'],'TS001',aipsclass='SPLAT',fsiz=256,imaxoff=1.0,dofilt=True,indisk=1,frifiledir='./',w_max=500,utstart=0,utstop=0,utinc=0.03,zap=True)
#
fringemap(341,'S356',['DE601','DE602','DE603','DE604','DE605','FR606','SE607','UK608'],'ST001',aipsclass='FITS',fsiz=256,imaxoff=10.0,dofilt=True,indisk=1,frifiledir='./',w_max=5000,utstart=0,utstop=0,utinc=0,zap=True)
#
#fringemap(340,'SBFRMAP',['DE602','DE603','DE604','DE605','FR606','SE607','UK608'],'DE601',aipsclass='FITS',fsiz=256,imaxoff=10.0,dofilt=True,indisk=1,frifiledir='./',w_max=5000,utstart=0,utstop=0,utinc=0,zap=True)
#
#fringemap(340,'J1153',['DE601','DE602','DE603','DE604','DE605','FR606','SE607','UK608'],'TS001',aipsclass='FITS',fsiz=256,imaxoff=3.0,dofilt=True,indisk=1,frifiledir='./',w_max=500,utstart=0.625,utstop=0.645,utinc=0.005,zap=True)
