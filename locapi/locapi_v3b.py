#!/usr/bin/env python
#### locapi: connection script for long baseline calibration
#    version history:
#     v1 NJ 2013.04.05
#     v2 NJ 2013.04.17 from busy week 17, major changes by AD
#     v3 AD 2013.04.25 including checkgaps.py permanently and fixed naming.
#     v3b NJ 2015.06.30 more extensively parallelized using multiprocessing
#             library; option to phase up just superterp or allcore added
#assumptions:
#(1) Phase-up cal sources are in caldir with non-variable prefix calprf
#    and suffix calsuf; calibrator is in beam calbea; standard LOFAR
#    file-names; calibrator is at position calra, caldec
#(2) Target is in srcdir with filename prefix srcprf and suffix srcsuf
#(3) User wants to shift and average around position ra,dec
#
# Uses easy_parallize subroutine from Scott Sievert
# scottsievert.github.io/blog/2014/07/30/simple-python-parallelism
#
# subroutines:
#  0  mkpar(srcbea)  ---> makes parset files
#  1  avcal(caldir,calprf,calsuf,fstep,tstep)  ---> averages cal data
#  2  mkcal(calra,caldec) ---> runs BBS on cal data for CS stations
#  3  apcal(srcdir,srcprf,srcsuf,wrkdir) ---> runs BBS to apply CS solutions
#  4  shift(srcdir,srcprf,srcsuf,ra,dec,src) ---> 
#                      shifts and averages to target
#  ### if AD's checkgaps.py is present, will be run at this point ###
#  5  comb(nsub,src) ---> makes IFs consisting of nsub subbands, convert->FITS
#              nb! ra, dec, src must be arrays of same length. shift and
#              comb run once for each source.
# Output files: PP_Txxx.fits where xxx is the source name. These can be
#               further processed with lofipi
#
import os,glob,sys,re,numpy as np,multiprocessing#,pyrap.tables as pt
#### change the following things:
sbegin=5     # begin at subroutine sbegin (see above); or 100 for all sources in beam
ncores=int(os.popen('nproc').read())
caldir='./survey_cal/'
calprf='L229585_'
calsuf='_uv.dppp.MS'
calbea='BEAM_0'
calra = 123.400137917
caldec= 48.2172222222
srcdir='./survey_data/'
srcprf='L229587_'
srcbea='BEAM_0'
srcsuf='_uv.dppp.MS'
src   =['B1150+497','3C266']
ra    =['11h53m24.4s','11h45m43.4s']
dec   =['49d31m08.0s','49d46m08.0s']
mscorpol='./mscorpol.py'
nsub  = 10     # subbands per IF
fra   = 242.30541   # used for auto-search of field
fdec  = 26.69139
doflag = False
superterp = '[\'CS002HBA0\',\'CS002HBA1\',\'CS003HBA0\',\'CS003HBA1\',\'CS004HBA0\',\'CS004HBA1\',\'CS005HBA0\',\'CS005HBA1\',\'CS006HBA0\',\'CS006HBA1\',\'CS007HBA0\',\'CS007HBA1\']'
allcore = '\'CS*\''
tiedstation = superterp   # choose either superterp or allcore

#### probably leave the rest alone

wrkdir = os.getcwd()+'/'
if doflag:
    a = np.sort(glob.glob(caldir+calprf+'*'+calsuf))
    a = np.append(a,np.sort(glob.glob(srcdir+srcprf+'*'+srcsuf)))
    for i in a:
        f=open('NDPPP.parset','w')
        f.write('msin='+i+'\n')
        f.write('msout=\n')
        f.write('steps=[flag,count]\n')
        f.write('flag.type = aoflagger\n')
        f.close()
    os.system('NDPPP NDPPP.parset')

# page down to ENDMKPAR

def parallel_function(f):
    def easy_parallize(f, sequence):
        from multiprocessing import Pool
        ncores=int(os.popen('nproc').read())  # work out how many cores
        pool = Pool(processes=ncores) # depends on available cores
        print f
        print sequence
        result = pool.map(f, sequence) # for i in sequence: result[i] = f(i)
        cleaned = [x for x in result if not x is None] # getting results
        cleaned = np.asarray(cleaned)
        pool.close() # not optimal! but easy
        pool.join()
        return cleaned
    from functools import partial
    return partial(easy_parallize, f)

def mkpar (srcbea):
    f=open('mkcal.parset','w')
    f.write('Strategy.ChunkSize = 0\n')   # changed
    f.write('Strategy.Steps = [solve]\n')
    f.write('Step.solve.Baselines = CS*&RS*\n')
    f.write('Step.solve.Operation = SOLVE\n')
    f.write('Step.solve.Model.Sources = [*]\n')
    f.write('Step.solve.Model.Gain.Enable = T\n')
    f.write('Step.solve.Model.DirectionalGain.Enable = F\n')
#   we are using mscorpol later, which has the beam, so disable beam here
    f.write('Step.solve.Model.Beam.Enable = F\n')
    f.write('Step.solve.Model.Cache.Enable = T\n')
    f.write('Step.solve.Solve.Parms = ["Gain:0:0:*","Gain:1:1:*"]\n')
    f.write('Step.solve.Solve.CellSize.Freq = 0\n')
    f.write('Step.solve.Solve.CellSize.Time = 1\n')
    f.write('Step.solve.Solve.CellChunkSize = 0\n')   # changed
    f.write('Step.solve.Solve.Options.MaxIter = 200\n')
    f.write('Step.solve.Solve.Options.EpsValue = 1e-9\n')
    f.write('Step.solve.Solve.Options.EpsDerivative = 1e-9\n')
    f.write('Step.solve.Solve.Options.ColFactor = 1e-9\n')
    f.write('Step.solve.Solve.Options.LMFactor = 1.0\n')
    f.write('Step.solve.Solve.Options.BalancedEqs = F\n')
    f.write('Step.solve.Solve.Options.UseSVD = T\n')
    f.close()
    f=open('phcal.parset','w')
    f.write('Strategy.InputColumn = DATA\n')
    f.write('Strategy.TimeRange = []\n')
    f.write('Strategy.ChunkSize = 600\n')
    f.write('Strategy.UseSolver = F\n')
    f.write('Strategy.Steps = [solve]\n')
    f.write('Step.solve.Baselines = CS*&\n')
    f.write('Step.solve.Operation = SOLVE\n')
    f.write('Step.solve.Model.Sources = [*]\n')
    f.write('Step.solve.Model.Cache.Enable = T\n')
    f.write('Step.solve.Model.Gain.Enable = T\n')
    f.write('Step.solve.Model.Bandpass.Enable = F\n')
    f.write('Step.solve.Model.DirectionalGain.Enable = F\n')
    f.write('Step.solve.Model.Beam.Enable = F\n')
    f.write('Step.solve.Model.Ionosphere.Enable = F\n')
    f.write('Step.solve.Solve.Mode = COMPLEX\n')
    f.write('Step.solve.Solve.Parms = ["Gain:*"]\n')
    f.write('Step.solve.Solve.ExclParms = []\n')
    f.write('Step.solve.Solve.CalibrationGroups = []\n')
    f.write('Step.solve.Solve.CellSize.Freq = 0\n')
    f.write('Step.solve.Solve.CellSize.Time = 30\n')
    f.write('Step.solve.Solve.CellChunkSize = 10\n')
    f.write('Step.solve.Solve.PropagateSolutions = T\n')
    f.write('Step.solve.Solve.Options.MaxIter = 50\n')
    f.write('Step.solve.Solve.Options.EpsValue = 1e-9\n')
    f.write('Step.solve.Solve.Options.EpsDerivative = 1e-9\n')
    f.write('Step.solve.Solve.Options.ColFactor = 1e-9\n')
    f.write('Step.solve.Solve.Options.LMFactor = 1.0\n')
    f.write('Step.solve.Solve.Options.BalancedEqs = F\n')
    f.write('Step.solve.Solve.Options.UseSVD = T\n')
    f.close()
    f=open('apcal.parset','w')
    f.write('Strategy.InputColumn = DATA\n')
    f.write('Strategy.TimeRange = []\n')
    f.write('Strategy.ChunkSize = 600\n')
    f.write('Strategy.UseSolver = F\n')
    f.write('Strategy.Steps = [correct]\n')
    f.write('Step.correct.Baselines = CS*&CS*\n')
    f.write('Step.correct.Operation = CORRECT\n')
    f.write('Step.correct.Model.Sources = []\n')
    f.write('Step.correct.Model.Gain.Enable = T\n')
    f.write('Step.correct.Model.Beam.Enable = F\n')
    f.write('Step.correct.Model.Phasors.Enable = F\n')
    f.write('Step.correct.Output.Column = CORRECTED_DATA\n')
    f.close()

# ENDMKPAR

def avcal (caldir, calprf, calsuf,fstep=16,tstep=4):
    a = np.sort(glob.glob(caldir+calprf+'*'+calsuf))
    os.system ('rm -fr PP_M*')
    for i in a:
        f=open('NDPPP.parset','w')
        f.write('msin = '+i+'\n')
        inew = i.replace(caldir+calprf,'PP_M').replace(calsuf,'.ms')
        f.write('msout = '+wrkdir+inew+'\n')
        f.write('steps = [avg]\n')
        f.write('avg.type = average\n')
        f.write('avg.freqstep = %d\n'%int(fstep))
        f.write('avg.timestep = %d\n'%int(tstep))
        f.close()
        os.system('NDPPP NDPPP.parset')

def mkcal_thread (fname):
    os.system('calibrate-stand-alone -v -n -f '+fname+' mkcal.parset PP_M.mod')
    return True

def mkcal (calra,caldec):
    a = np.sort(glob.glob('PP_M*ms'))
    os.system('python vgsm.py PP_M.mod '+str(calra)+' '+str(caldec)+' 1.0')
    mkcal_thread.parallel = parallel_function(mkcal_thread)
    parallel_result = mkcal_thread.parallel(a)
    f = open('parmdbm_command','w')
    for k in a:
        f.write('open tablename=\''+k+'/instrument\'\n')
        f.write('export Gain* tablename=\''+k.replace('ms','table')+'\'\n')
    f.close()
    os.system('parmdbm < parmdbm_command')

def apcal_thread (k):
    os.system('calibrate-stand-alone -v -n -f --parmdb %s %s apcal.parset PP_M.mod'\
                 %(k[1],k[0]))
    return True

def apcal (srcdir,srcprf,srcsuf,wrkdir):
    a1 = np.sort(glob.glob(srcdir+srcprf+'*'+srcsuf))
    tab1 = np.sort(glob.glob('PP_M*.table'))
    anum = np.copy(a1)
    tabnum = np.copy(tab1)
    for i in range(len(anum)):
        anum[i] = anum[i].replace(srcdir,'Z')
        anum[i] = anum[i].replace(srcprf,'Z')
        anum[i] = anum[i].replace(srcsuf,'Z')
        anum[i] = re.search(r'\d+',anum[i]).group()
    for i in range(len(tabnum)):
        tabnum[i] = re.search(r'\d+',tabnum[i]).group()
    anum = np.intersect1d(anum,tabnum)
    a = np.array([]); tab = np.array([])
    for i in range(len(anum)):
        for j in a1:
            test = j.replace(srcdir,'Z')
            test = test.replace(srcprf,'Z')
            test = test.replace(srcsuf,'Z')
            test = re.search(r'\d+',test).group()
            if test==anum[i]:
                a=np.append(a,j)
        for j in tab1:
            test = j.replace(srcdir,'Z')
            test = test.replace(srcprf,'Z')
            test = test.replace(srcsuf,'Z')
            test = re.search(r'\d+',test).group()
            if test==anum[i]:
                tab=np.append(tab,j)

    apcal_thread.parallel = parallel_function(apcal_thread)
    k = np.dstack((a,tab))[0]
    parallel_result = apcal_thread.parallel(k)

def shift_thread (k):
    sbnum = k[0].split(srcprf)[1].split(srcsuf)[0]
    f=open('NDPPP'+sbnum+'.parset','w')
    f.write('msin = '+k[0]+'\n')
    outfile = k[0].replace(srcdir+srcprf,'PP_N').replace(srcsuf,k[1]+'.ms')
    f.write('msout = '+outfile+'\n')
    f.write('msin.datacolumn = CORRECTED_DATA\n')
    f.write('steps = [shift,avg, adder, filter]\n')
    f.write('shift.type=\'phaseshift\'\n')
    f.write('shift.phasecenter=[\''+k[2]+'\',\''+k[3]+'\']\n')
    f.write('avg.type = squash\n')
    f.write('avg.freqstep = %s\n'%k[4])
    f.write('avg.timestep = %s\n'%k[5])
    f.write('adder.type = \'stationadder\'\n')
    f.write('adder.stations = {TS001:%s}\n' % k[6])
    f.write('filter.type = \'filter\'\n')
    f.write('filter.baseline = \'!CS*&*\'\n')
    f.close()
    os.system('NDPPP NDPPP'+sbnum+'.parset')
    return True

def shift (srcdir,srcprf,srcsuf,ra,dec,src,tiedstation,fstep=1,tstep=1):
    a = np.sort(glob.glob(srcdir+srcprf+'*'+srcsuf))
    os.system('rm -fr PP_N*')
    for i_s in range (len(ra)):
        k = np.column_stack ((a, len(a)*[src[i_s]], len(a)*[ra[i_s]],\
             len(a)*[dec[i_s]], len(a)*[str(fstep)], len(a)*[str(tstep)],\
             len(a)*[tiedstation]))
        shift_thread.parallel = parallel_function(shift_thread)
        parallel_result = shift_thread.parallel (k)


def writeNDPPP(msin_list, src, fn):
     f = open('NDPPP.parset', 'w')
     f.write('msin = [')
     for j in range(len(msin_list)):
         f.write(msin_list[j])
         f.write(',' if j < len(msin_list) - 1 else ']\n')
         pass
     f.write('msin.datacolumn = DATA\n')
     f.write('msin.missingdata = True\n')
     f.write('msin.orderms = False\n')
     f.write('msout = PP_C'+str(fn).rjust(2,'0')+'_'+str(src)+'.ms\n')
     f.write('steps=[]\n')
     f.close()
     os.system('NDPPP')


def comb_thread (k):
    os.system(mscorpol+' -f '+k)
    return True

def comb (nsub,s):
    os.system('rm -fr PP_C*'+s+'.ms')
    filelist = np.sort(glob.glob('PP_N*'+s+'.ms'))
    fn       = 0
    chan1    = []
    lastchan = []
    if nsub <= 0:
        print 'ERROR: Number of subbands have to be positive.'
        sys.exit()
    for i in range(0,len(filelist)):
        chan1.append(float(os.popen('msoverview in='+filelist[i]+' | tail -2 | head -1 | awk \'{print $4}\'').readlines()[0]))
        print 'SB:', filelist[i][6:9],  '| 1st channel:', chan1[i]
        pass
    chan1_diff = np.diff(chan1)
    chan1_mult = (chan1_diff / chan1_diff.min()).astype(int) - 1
    check = all([chan1_mult[i] == 0 for i in range(len(chan1_mult))])
    if not check:
        print 'There have been missing subbands or they have not been evenly distributed. Will add "ghost" subbands.'
    offset = 0
    i = 0
    msin_list = []
    while i < len(filelist):
        if len(msin_list) == nsub:
            writeNDPPP(msin_list, s, fn); fn+=1
            msin_list = [] 
        if offset == 0:
            msin_list.append(filelist[i])
        else:
            msin_list.append('ghost')
            offset -= 1
            continue 
        if i == len(filelist) - 1:
            for ghost in range(nsub - len(msin_list)):
                msin_list.append('ghost')
            break 
        if (chan1_mult[i] > (nsub - len(msin_list))): 
            offset = chan1_mult[i] - (nsub - len(msin_list)) 
        for ghost in range(min(chan1_mult[i],(nsub - len(msin_list)))):
            msin_list.append('ghost') 
            if len(msin_list) == nsub:
                offset = chan1_mult[i] - (ghost + 1)
                break 
        i += 1
    writeNDPPP(msin_list, s, fn)
    os.system('rm -fr PP_T*'+s+'.ms')
    a=glob.glob('PP_C*'+s+'.ms')
    comb_thread.parallel = parallel_function(comb_thread)
    parallel_result = comb_thread.parallel(a)
    pt.msconcat (a, 'PP_T'+src+'.ms')
    os.system('ms2uvfits in=PP_T'+src+'.ms out=PP_T'+src+'.fits writesyscal=F')    

# -------------- main script -----------------
if sbegin <1:
    mkpar(srcbea)
if sbegin <2:
    avcal(caldir,calprf,calsuf)
if sbegin <3:
    mkcal(calra,caldec)
if sbegin <4:
    apcal(srcdir,srcprf,srcsuf,wrkdir)
if sbegin <5:
    shift(srcdir,srcprf,srcsuf,ra,dec,src,tiedstation)
if sbegin <6:
    for s in src:
        comb(nsub,s)

if sbegin == 100:         # do all sources in beam - say within 1 deg
    ra=[]; dec=[]; src=[]; field = 1.0
    nvss = np.load('NVSS_strong.npy')
    for i in nvss:
        decdiff = abs(i[1]-fdec)
        if decdiff > field:
            continue
        radiff = abs(i[0]-fra) * np.cos(np.deg2rad(fdec))
        if radiff > field:
            continue
        if np.hypot (radiff,decdiff) > field:
            continue
        rahr = int(i[0]/15)
        ramin = 60.*(i[0]/15-rahr)
        rasec = 60.*(ramin-int(ramin))
        ramin = int(ramin)
        decsgn = '+' if i[1]>0.0 else '-'
        deci = abs(i[1])
        decdeg = int(deci)
        decmin = 60.*(deci-decdeg)
        decsec = 60.*(decmin-int(decmin))
        decmin = int(decmin)
        ra.append('%02d'%rahr +'h'+'%02d'%ramin +'m'+ '%05.2f'%rasec +'s')
        dec.append(decsgn+'%02d'%decdeg+'d'+'%02d'%decmin+'m'+'%04.1f'%decsec+'s')
        src.append('S%02d'%rahr+'%02d'%ramin+'%02d'%int(rasec)+\
                   decsgn+'%02d'%decdeg+'%02d'%decmin+'%02d'%int(decsec))

    mkpar(srcbea)
    avcal(caldir,calprf,calsuf)
    mkcal(calra,caldec)
    apcal(srcdir,srcprf,srcsuf,wrkdir)
    shift(srcdir,srcprf,srcsuf,ra,dec,src,tiedstation)
    for s in src:
        comb(nsub,s)

