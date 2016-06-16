#!/usr/bin/env python
#### locapi: connection script for long baseline calibration
#    version history:
#     v1 NJ 2013.04.05
#     v2 NJ 2013.04.17 from busy week 17, major changes by AD
#     v3 AD 2013.04.25 including checkgaps.py permanently and fixed naming.
#     v3b NJ 2015.06.30 more extensively parallelized using multiprocessing
#             library; option to phase up just superterp or allcore added
#     v4 LM 2016.03.14 updated to use NDPPP; added argparser so editing 
#	      script is no longer necessary
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
import argparse

####### HELPER FUNCTIONS #######

def get_inttime(msname):

    os.system("taql 'select distinct EXPOSURE from %s' > inttime.log"%(msname))
    f = open('inttime.log','r')
    lines = f.readlines(3)
    inttime = np.int(np.float(lines[2].replace('\n','')))
    print 'Integration time: ',inttime,' seconds'
    os.system('rm inttime.log')
    return( inttime )

def get_nchans(msname):
    os.system("taql 'select distinct NUM_CHAN from %s/SPECTRAL_WINDOW' > intfreq.log"%(msname))
    f = open('intfreq.log','r')
    lines = f.readlines(3)
    nchans = np.int(np.float(lines[2].replace('\n','')))
    print 'Num channels: ',nchans
    os.system('rm intfreq.log')
    return( nchans )

## parallelization
def parallel_function(f):
    def easy_parallize(f, sequence):
        from multiprocessing import Pool
        #ncores=int(os.popen('nproc').read())  # work out how many cores
        #ncores=1
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

## write parsets
def mkpar (srcbea):

    f=open('mkcal.parset','w')
    f.write('msin.datacolumn=DATA\n')
    f.write('msout=.\n')
    f.write('steps=[filter,gaincal]\n')
    f.write('filter.baseline=CS*&RS*\n')
    f.write('gaincal.caltype=diagonal\n')
    f.write('gaincal.solint=1\n')
    f.write('gaincal.maxiter=200\n')
    f.write('gaincal.sources=[*]\n')
    f.write('gaincal.usebeammodel=False\n')
    f.close()

    f=open('phcal.parset','w')
    f.write('msin.datacolumn=CORRECTED_DATA\n')
    f.write('msout=.\n')
    f.write('steps=[filter,gaincal]\n')
    f.write('filter.baseline=CS*&\n')
    f.write('gaincal.caltype=phaseonly\n')
    f.write('gaincal.solint=15\n')
    f.write('gaincal.maxiter=50\n')
    f.write('gaincal.sources=[*]\n')
    f.write('gaincal.usebeammodel=False\n')
    f.close()

    f=open('apcal.parset','w')
    f.write('msin.datacolumn=DATA\n')
    f.write('msout=.\n')
    f.write('msout.datacolumn=CORRECTED_DATA\n')
    f.write('steps=[applycal]\n')
    f.write('applycal.correction=gain\n')
    f.close()

    f=open('apphcal.parset','w')
    f.write('msin.datacolumn=CORRECTED_DATA\n')
    f.write('msout=.\n')
    f.write('msout.datacolumn=CORRECTED_PHS_DATA\n')
    f.write('steps=[applycal]\n')
    f.write('applycal.correction=gain\n')
    f.close()


## run averaging of calibrator
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

## parallelization of calibrator gains
def mkcal_thread (fname):
    os.system('NDPPP mkcal.parset msin=%s gaincal.sourcedb=PP_M.mod gaincal.parmdb=%s'%(fname, fname.replace('.ms','.table')))
    return True

## find the calibrator gains
def mkcal (calra,caldec,calmod):
    ## this needs to be changed because I don't have vgsm.py on the system
    a = np.sort(glob.glob('PP_M*ms'))
    #os.system('python vgsm.py PP_M.mod '+str(calra)+' '+str(caldec)+' 1.0')
    #os.system('cp -r %s PP_M.skymod'%(calmod))
    ## NDPPP requires the output of makesourcedb instead of a text file
    os.system('makesourcedb in=%s out=PP_M.mod format="<"'%(calmod))
    mkcal_thread.parallel = parallel_function(mkcal_thread)
    parallel_result = mkcal_thread.parallel(a)
#    f = open('parmdbm_command','w')
#    for k in a:
#        f.write('open tablename=\''+k+'/instrument\'\n')
#        f.write('export Gain* tablename=\''+k.replace('ms','table')+'\'\n')
#    f.close()
#    os.system('parmdbm < parmdbm_command')

## parallelization of fixing gains
def fixcal_thread (tname):
    os.system('python %s %s'%(addIS,tname))
    os.system('python %s %s'%(updateIS,tname))
    return True

## fix the international station gains
def fixcal ():
    a = np.sort(glob.glob('PP_MS*table'))
    fixcal_thread.parallel = parallel_function(fixcal_thread)
    parallel_result = fixcal_thread.parallel(a)
    
## parallelization of applying gains
def apcal_thread (k):
    #os.system('calibrate-stand-alone -v -n -f --parmdb %s %s apcal.parset PP_M.mod'\
    os.system('NDPPP apcal.parset msin=%s applycal.parmdb=%s'%(k[0],k[1]))
    return True

## apply calibrator gains to target
def apcal (srcdir,srcprf,srcsuf,wrkdir):
    a1 = np.sort(glob.glob(srcdir+srcprf+'*'+srcsuf))
    tab1 = np.sort(glob.glob('def_M*.table'))
    anum = np.copy(a1)
    tabnum = np.copy(tab1)
    for i in range(len(anum)):
        anum[i] = anum[i].replace(srcdir,'Z')
        anum[i] = anum[i].replace(srcprf,'Z')
        anum[i] = anum[i].replace(srcsuf,'Z')
        anum[i] = re.search(r'\d+',anum[i]).group()
    for i in range(len(tabnum)):
        tabnum[i] = re.search(r'\d+',tabnum[i]).group()
    tmp = np.intersect1d(anum,tabnum)
    if len(tmp) == 0:
	## this is the case if the subband numbers are off by 244
	for t in tabnum:
	    ss = 'mv def_MSB%s.table def_MSB%s.table'%(t, str(int(t)-244))
	    os.system(ss)
	tab1 = np.sort(glob.glob('def_M*.table'))
	tabnum = np.copy(tab1)
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

## parallelization of fixing gains
def phcal_thread (fname):
    os.system('NDPPP phcal.parset msin=%s gaincal.sourcedb=PP_phs.mod gaincal.parmdb=%s'%(fname, fname.replace('.MS','.phtab')))
    os.system('NDPPP apphcal.parset msin=%s applycal.parmdb=%s'%(fname, fname.replace('.MS','.phtab')))
    return True

def phcal (srcdir,srcprf,srcsuf,srcmod):
    a = np.sort(glob.glob(srcdir+srcprf+'*'+srcsuf))
    phcal_thread.parallel = parallel_function( phcal_thread )
    os.system('makesourcedb in=%s out=PP_phs.mod format="<"'%(srcmod))
    parallel_result = phcal_thread.parallel(a)

## parallelization of phaseshift, average, stationadd, filter
def shift_thread (k):
    #sbnum = k[0].split(srcprf)[1].split(srcsuf)[0]
    sbnum = k[0][k[0].find('SB')+2:k[0].find('SB')+5]
    f=open('NDPPP'+sbnum+'.parset','w')
    f.write('msin = '+k[0]+'\n')
    #outfile = k[0].replace(srcdir+srcprf,'PP_N').replace(srcsuf,k[1]+'.ms')
    outfile = 'PP_NSB'+sbnum+k[1]+'.ms'
    f.write('msout = '+outfile+'\n')
    f.write('msin.datacolumn = CORRECTED_PHS_DATA\n')
    f.write('steps = [shift,avg, adder, filter]\n')
    f.write('shift.type=\'phaseshift\'\n')
    f.write('shift.phasecenter=[\''+k[2]+'\',\''+k[3]+'\']\n')
    f.write('avg.type = squash\n')
    f.write('avg.freqstep = %s\n'%k[4])
    f.write('avg.timestep = %s\n'%k[5])
    f.write('adder.type = \'stationadder\'\n')
    f.write('adder.stations = {TS001:%s}\n' % k[6])
    f.write('filter.type = \'filter\'\n')
    f.write('filter.baseline = \'!CS*&&*\'\n')
    f.write('filter.remove = True \n')
    f.close()
    os.system('NDPPP NDPPP'+sbnum+'.parset')
    return True

## phaseshift, average, stationadd, filter
def shift (srcdir,srcprf,srcsuf,ra,dec,src,tiedstation,fstep=1,tstep=1):
    a = np.sort(glob.glob(srcdir+srcprf+'*'+srcsuf))
    os.system('rm -fr PP_N*')
    for i_s in range (len(ra)):
        k = np.column_stack ((a, len(a)*[src[i_s]], len(a)*[ra[i_s]],\
             len(a)*[dec[i_s]], len(a)*[str(fstep)], len(a)*[str(tstep)],\
             len(a)*[tiedstation]))
        shift_thread.parallel = parallel_function(shift_thread)
        parallel_result = shift_thread.parallel (k)

## write combine parset
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

## parallelization of combine
def comb_thread (k):
    os.system(mscorpol+' -f '+k)
    return True

## combine subbands
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
        #chan1.append(float(os.popen('msoverview in='+filelist[i]+' | tail -2 | head -1 | awk \'{print $4}\'').readlines()[0]))
        chan1.append(float(os.popen('taql "select CHAN_FREQ[0] from '+filelist[i]+'/SPECTRAL_WINDOW" | tail -1 ').readlines()[0]))
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

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--doflag',dest='doflag',action="store_true",default=False,help="Set this flag to run aoflagger")
    parser.add_argument('--allcore',dest='allcore',action="store_true",default=False,help="Combine all core stations (default is to combine superterp only)")
    parser.add_argument('-n','--ncores',dest='ncores',type=int,default=1,help="Number of cores to use in parallel (default=1)")
    parser.add_argument('-s','--sbegin',dest='sbegin',type=int,default=1,help="begin at subroutine sbegin (see above); or 100 for all sources in beam (default 1)")
    parser.add_argument('--nsub',dest='nsub',type=int,default=10,help="Number of subbands per IF (default 10)")
    parser.add_argument('-m','--mscorpol',dest='mscorpol',type=str,default='./mscorpol.py',help="Path to mscorpol.py (default ./mscorpol.py")
    parser.add_argument('--intgains',dest='intgains',type=str,default='./',help="Path to IntlStationGainCal directory that contains addIS.py and updateISGains.py")
    parser.add_argument('-f','--file_desc',dest='file_desc',type=str,default='locapi_file_desc.txt',help="Name of file description (default locapi_file_desc.txt)")


    args = parser.parse_args()
    doflag = args.doflag
    allcore = args.allcore
    global ncores
    ncores = args.ncores
    sbegin = args.sbegin
    nsub = args.nsub
    global mscorpol
    mscorpol = args.mscorpol
    intgains = args.intgains
    file_desc = args.file_desc

    global addIS
    addIS = os.path.join( intgains, "addIS.py" )
    global updateIS
    updateIS = os.path.join( intgains, "updateISGains.py" )

    ## set superterp or all core
    if allcore:
	tiedstation = '\'CS*\''
    else:
	tiedstation = '[\'CS002HBA0\',\'CS002HBA1\',\'CS003HBA0\',\'CS003HBA1\',\'CS004HBA0\',\'CS004HBA1\',\'CS005HBA0\',\'CS005HBA1\',\'CS006HBA0\',\'CS006HBA1\',\'CS007HBA0\',\'CS007HBA1\']'

    ## read in file description
    ## should be in the following format:
    #caldir='./survey_cal/'
    #calprf='L427092_'
    #calsuf='_uv.dppp.MS'
    #calbea='BEAM_0'
    #calra = 123.400137917
    #caldec= 48.2172222222
    #calmod='/net/para33/data2/lofar/models/3C147_SH.skymodel'
    #srcdir='./survey_data/'
    #srcprf='L427102_'
    #srcbea='BEAM_0'
    #srcsuf='_uv.dppp.MS'
    #src   =['4C4315']
    #ra    =['07h35m21.887s']
    #dec   =['43d44m20.32s']
    #srcmod = '/net/para13/data2/morabito/4C4315/LBA/models/4C4315_10x10_25.MOSAIC.bbs_7deg_0.486Jy_A1.5deg'
    with open( file_desc ) as f:
	lines = f.readlines()
    caldir = lines[0].split("=")[1].rstrip('\n').strip("'")
    calprf = lines[1].split("=")[1].rstrip('\n').strip("'")
    calsuf = lines[2].split("=")[1].rstrip('\n').strip("'")
    calbea = lines[3].split("=")[1].rstrip('\n').strip("'")
    calra = lines[4].split("=")[1].rstrip('\n').strip("'")
    caldec = lines[5].split("=")[1].rstrip('\n').strip("'")
    calmod = lines[6].split("=")[1].rstrip('\n').strip("'")
    srcdir = lines[7].split("=")[1].rstrip('\n').strip("'")
    srcprf = lines[8].split("=")[1].rstrip('\n').strip("'")
    srcbea = lines[9].split("=")[1].rstrip('\n').strip("'")
    srcsuf = lines[10].split("=")[1].rstrip('\n').strip("'")
    src = lines[11].split("=")[1].rstrip('\n').strip("[").strip("]").strip("'").split(',')
    ra = lines[12].split("=")[1].rstrip('\n').strip("[").strip("]").strip("'").split(',')
    dec = lines[13].split("=")[1].rstrip('\n').strip("[").strip("]").strip("'").split(',')
    srcmod = lines[14].split("=")[1].rstrip('\n').strip("'")

    ## if flagging is requested, run AOFlagger
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

    ## run the steps requested
    if sbegin == 1:
	## write parsets
        mkpar(srcbea)
    if sbegin == 2:
	## average calibrator
	a = np.sort(glob.glob(caldir+calprf+'*'+calsuf))
	freqstep = int( get_nchans(a[0]) / 4 )
	timestep = int( get_inttime(a[0]) / 4 )
	print( "Averaging to 4 ch/SB, 4 sec" )
        avcal(caldir,calprf,calsuf,fstep=freqstep,tstep=timestep)
    if sbegin == 3:
	## calibrator gains
        mkcal(calra,caldec,calmod)
	## fix the calibrator gains
	fixcal()
    if sbegin == 4:
	## apply calibrator gains
        apcal(srcdir,srcprf,srcsuf,wrkdir)
    if sbegin == 5:
	## phase only calibration
	phcal (srcdir,srcprf,srcsuf,srcmod)
    if sbegin == 6:
	## phaseshift, average, combine stations, filter
        shift(srcdir,srcprf,srcsuf,ra,dec,src,tiedstation)
    if sbegin == 7:
	## combine sources
        for s in src:
            comb(nsub,s)


    ## I HAVEN'T TOUCHED THIS PART

    fra   = 242.30541   # used for auto-search of field
    fdec  = 26.69139

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
        mkcal(calra,caldec,calmod)
        apcal(srcdir,srcprf,srcsuf,wrkdir)
        shift(srcdir,srcprf,srcsuf,ra,dec,src,tiedstation)
        for s in src:
            comb(nsub,s)

if __name__ == "__main__":
    main()

