#from AIPS import AIPS, AIPSDisk
#from AIPSTask import AIPSTask, AIPSList, AIPSMessageLog
#from AIPSData import AIPSUVData, AIPSImage, AIPSCat
#from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
#from lofipi_aips import *
import numpy as np; from numpy import fft; from scipy import ndimage
import re,sys,pickle,os,glob,time,warnings,pyfits,matplotlib,time,functools
from functools import partial
from matplotlib import pyplot as plt; from scipy.ndimage import measurements
from pyfits import getdata,getheader
from correlate import *
try:
    import astLib; from astLib import astCoords
except:
    print 'Operating without astLib'
try:
    import casac,inspect,string,sys
except:
    print 'Operating without casa libraries'

ncores = int(os.popen('nproc').read())


def parallel_function(f):   # Scott Sievert easy_parallelize code
    def easy_parallize(f, sequence):
        from multiprocessing import Pool
#        ncores=int(os.popen('nproc').read())  # work out how many cores
        ncores=4
        pool = Pool(processes=ncores) # depends on available cores
        print f
        print sequence
        result = pool.map(f, sequence) # for i in sequence: result[i] = f(i)
        cleaned = [x for x in result if not x is None] # getting results
        cleaned = np.asarray(cleaned)
        pool.close() # not optimal! but easy
        pool.join()
        return cleaned
    return partial(easy_parallize, f)

def lshift_thread (k):
    os.system('NDPPP lshift_%03d' % k)

def create_lobos_ms (infile='',ra='',dec='',nP=3,radius=2.5):
    infile='L401323_SB349_uv.dppp.MS'
    ra,dec,radius,nP='13:40:00','55:00:00',2.5,3
    try:    # ra, dec input can either be decimal or hh:mm:ss, dd:mm:ss
        ra, dec = 1.0*ra, 1.0*dec
        hexra = astCoords.decimal2hms(ra,':')
        hexdec = astCoords.decimal2dms(dec,':')
    except:
        hexra, hexdec = ra, dec
        ra = astCoords.hms2decimal(ra,':')
        dec = astCoords.dms2decimal(dec,':')
    if not os.path.isfile ('lobos_stats.sum'):
        os.system ('wget http://www.jb.man.ac.uk/~njj/lobos_stats.sum')
    lobos = np.loadtxt('lobos_stats.sum',dtype='S')
    for l in lobos:
        new = np.array([astCoords.hms2decimal(l[1],':'),\
                        astCoords.dms2decimal(l[2],':')])
        try:
            lobos_coord = np.vstack((lobos_coord,new))
        except:
            lobos_coord = np.copy(new)
    a = correlate(np.array([[ra,dec]]),0,1,lobos_coord,0,1,radius)
    nlobos = 0
    os.system('rm lshift_*')
    for i in np.asarray(a[:,1],dtype='int'):
        if lobos[i][5].count('P')>=nP:
            oroot = infile.replace('uv.dppp.MS','')
            namera = lobos[i,1].replace(':','').split('.')[0]
            namedec = lobos[i,2].replace(':','').split('.')[0]
            dpppra = lobos[i,1].replace(':','h',1).replace(':','m',1)+'s'
            dpppdec = lobos[i,2].replace(':','d',1).replace(':','m',1)+'s'
            hemisphere = '-' if '-' in hexdec else '+'
            f = open('lshift_%03d'%nlobos, 'w')
            f.write('msin = '+infile+'\n')
            f.write('msout = '+oroot+namera+hemisphere+namedec+'.ms\n')
            f.write('msin.datacolumn = DATA\n')
            f.write('steps = [shift,adder,filter,avg]\n')
            f.write('shift.type = \'phaseshift\'\n')
            f.write('shift.phasecenter = [\''+dpppra+'\',\''+dpppdec+'\']\n')
            f.write('adder.type = \'stationadder\'\n')
            f.write('adder.stations = {TS001:\'CS*\'}\n')
            f.write('filter.type = \'filter\'\n')
            f.write('filter.baseline = \'!CS*&*\'\n')
            f.write('avg.type = average\n')
            f.write('avg.freqstep = 16\n')
            f.write('avg.timestep = 20\n')
            f.close()
            nlobos+=1

    lshift_thread.parallel = parallel_function(lshift_thread)
    k = range(nlobos)
    parallel_result = lshift_thread.parallel(k)


def get_closure_phase(infile='L401323_SB349_uv.dppp.MS',\
                 triangle = ['TS001','DE601HBA','DE605HBA']):
    a=inspect.stack()
    stacklevel=0
    for k in range(len(a)):
        if (string.find(a[k][1],'ipython console')>0):
            stacklevel=k
    myf=sys._getframe(stacklevel).f_globals
    myf['__last_task']='mytask'
    myf['taskname']='mytask'
    tb=myf['tb']
    oroot = infile.split('uv')[0]
    for lfile in np.sort(glob.glob(oroot+'*ms')):
        os.system('ms2uvfits in='+lfile+' out='+lfile.replace('ms','fits')+' writesyscal=F')
        if lfile == infile:
            continue
        tb.open(lfile+'/ANTENNA')
        names = tb.getcol('NAME')
        trnum = []
        for itr in range(3):
            trnum.append(np.argwhere(names==triangle[itr])[0][0])
        tb.close()
        trnum.sort()
        tb.open(lfile)
        ant1 = tb.getcol('ANTENNA1')
        ant2 = tb.getcol('ANTENNA2')
        data = tb.getcol('DATA')
        ph12 = +np.angle(data[0,0,(ant1==trnum[0])&(ant2==trnum[1])])
        ph23 = +np.angle(data[0,0,(ant1==trnum[1])&(ant2==trnum[2])])
        ph31 = -np.angle(data[0,0,(ant1==trnum[0])&(ant2==trnum[2])])
        clph = ph12+ph23+ph31
        np.putmask(clph,clph>np.pi,clph-2.*np.pi)
        np.putmask(clph,clph<-np.pi,clph+2.*np.pi)
#        np.savetxt(lfile.replace('ms','txt'),np.unwrap(clph))
        np.savetxt(lfile.replace('ms','txt'),clph)

def plot_closure_phase (infile):
    nx,ny = 2,5
    a = np.sort(glob.glob(infile.split('uv')[0]+'*.txt'))
    plotnum = 0
    for aname in a:
        i=np.loadtxt(aname)
        np.putmask(i,i>np.pi,i-2.*np.pi)
        np.putmask(i,i<-np.pi,i+2.*np.pi)
        plt.subplot(ny,nx,1+plotnum%(nx*ny),xticks=[])
        matplotlib.rcParams.update({'font.size':6})
        plt.plot(i,'b+',markersize=0.4)
        plt.title(aname.split('_')[2].split('.')[0])
        plotnum+=1
        if not plotnum%(nx*ny):
            plt.savefig('lobos_closure_phase_%d.png'%(plotnum/(nx*ny)),\
                        bbox_inches='tight',clobber=True)
            plt.clf()
            
    plt.savefig('lobos_closure_phase_%d.png'%(1+plotnum/(nx*ny)),\
                        bbox_inches='tight',clobber=True)

#msin = L401323_SB349_uv.dppp.MS
#msin = TEMP.MS
#msout = STSHIFT.MS
#msin.datacolumn = DATA
#steps = [adder]
#steps = [shift,adder,filter]
#shift.type = 'phaseshift'
#shift.phasecenter = ['13h27m37.2s','55d04m06.0s']
#adder.type = 'stationadder'
#adder.stations={TS002:['CS002*','CS003*','CS004*','CS005*','CS006*','CS007*']}
#adder.stations = {TS001:'CS*'}
#filter.type = 'filter'
#filter.baseline = '!CS*&*'
