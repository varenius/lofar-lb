# Script to write PARSET files to combined multiple LOFAR SB files into one MS file per source per scan.
# After running the parsets generate by this file, run the writeCASASCRIPT-file to generate python scripts
# for CASA, which can be used to concatenate all files to one MS per source (instead of per source per scan).
# The CASA file will also given instructions for conversion to circular, as well as send
# This file has a lot of variations commented out below, to sum stations in different ways, after calibrating on 3C-scans.
# If not combing stations, these blocks of code (and subsequent DPPP parset files) can be skipped.
# Run this script as "python script.py msdir" where msdir contains the gazillion files downloaded from the LTA
import sys, os, re
import glob

msdir = sys.argv[1] # Directory containing all the gazillion ms files from the LOFAR LTA, i.e. one file per source per scan per subband.

ntot = 162 # Number of subbands per beam (Check your data, this varies. Total is 488, so if 3 beams then you probably have skipped two and have 486/3=162 per beam)
# Function to calculate which beam, i.e. subband range, we are in given a specific SB number
def get_sblist(sb):
    sb1 = range(0, ntot)
    sb2 = range(ntot,2*ntot)
    sb3 = range(2*ntot, 3*ntot)
    if sb in sb1:
        return sb1
    elif sb in sb2:
        return sb2
    elif sb in sb3:
        return sb3
  
scriptpath = os.path.dirname(os.path.realpath(sys.argv[0]))
data = open(scriptpath+'/allfiles.dat') # File with info from LTA average pipelinesummary page
RFIchans = []  # Optional, can remove specific subbands now, but will also be removed by AOFLAGGER later.
calsource = '3C295'

lnums = []
runids = []
sources = []

for line in data:
 if not line.startswith('#'):
  # Hard coded indices reflect coumns in the LTA-infopage file. Edit as needed.
  items = line.split()
  sources.append(items[3].split('/')[0]) # Source name
  runids.append(items[3].split('/')[1].split('.')[0]) # Run ID
  lnums.append(items[5]) # L-number, used to cross-ref some calibration IDs

# Get number of closest 3C-scan for another L-number, i.e. a 3C scan within the same "runID"
def get_calnum(Lnum):
 lid = lnums.index(Lnum)
 runid = runids[lid]
 answer = None
 for i, rid in enumerate(runids):
  if (rid==runid and sources[i] == calsource):
   answer = lnums[i]
   break
 if answer == None:
  print 'Found no calibrator scan for Lnum ' + str(Lnum)
 else:
  return answer

for i, Lnum in enumerate(lnums):
 files = glob.glob(msdir + '*' + str(Lnum) + '*')
 nfiles = len(files)
 if nfiles==0:
  print 'WARNING: Expected '+str(ntot)+' files for source '+sources[i] + ' Lnum '  + str(Lnum)  + ' but found '+str(nfiles) + '. Skipping these parsets...'
 else:
  if nfiles!=ntot:
   print 'WARNING: Expected '+str(ntot)+' files for source '+sources[i] + ' Lnum '  + str(Lnum)  + ' but found '+str(nfiles) + '. Making parsets anyway...'
  else:
   print 'Found '+str(nfiles)+' as expected for source '+sources[i] + ' Lnum '  + str(Lnum)  + '. Making parsets...'
  calnum = get_calnum(Lnum)   
  sbs = []
  for f in files:
   fn = os.path.basename(f)
   SB = int(fn.split('_')[1][2:])
   sbs.append(SB)
  sblist = get_sblist(sbs[0])
  if sblist is None:
   print 'sblist None! Check list of subbands to include ' + str(sbs[0])

  # Write parset to concatenate all SB-files for a specific source, as well as apply beam and average.
  outfile = Lnum + '_' + sources[i] + '_PRECAL.parset'
  with open(outfile, 'w') as f:
   f.write('msin = [\n')
   for sb in sblist:
    # Check for extracted filename
    msfilename =  'L' + Lnum + '_SB' + str(sb).zfill(3) + '_uv.dppp.MS'
    fullpath1 = os.path.join(msdir,msfilename)
    if (sb-sblist[0]) in RFIchans:
     f.write('\'' + Lnum + '_SB'+str(sb).zfill(3) + '_REMOVED_DUE_TO_RFI\'')
    elif os.path.exists(fullpath1):
     f.write('\'' + fullpath1+ '\'')
    # Found nothing, mark as missing, will be flagged by NDPPP
    else:
     f.write('\'' + Lnum + '_SB'+str(sb).zfill(3) + '_MISSING\'')
    # Finish line. Important, if comma at last row this will make fake channel at the end
    # Since we use missing MS flag true.
    if sb == sblist[-1]:
     f.write(']\n')
    else:
     f.write(',\n')
   f.write('msin.missingdata=True\n')
   f.write('msin.datacolumn = DATA\n')
   f.write('msin.orderms=False\n')
   f.write('msout = ' + Lnum + '_'+sources[i]+'_COMB_AVGD.MS\n')
   #f.write('preflagger.baseline = \'CS501HBA&&*\'\n')
   #f.write('aoflagger.type = aoflagger\n')
   f.write('applybeam.type=applybeam\n')
   f.write('applybeam.usechannelfreq=true\n')
   f.write('applybeam.invert=true\n')
   f.write('applybeam.beammode = \"default\"\n')
   f.write('averager.freqstep = 4\n') # to 1ch/SB assuming 4ch/SB input
   f.write('averager.timestep=2\n') # To 8s assuming 4sec input
   #f.write('steps = [preflagger, aoflagger, applybeam]\n')
   #f.write('steps = [aoflagger, applybeam, averager]\n')
   f.write('steps = [applybeam, averager]\n')
  
  # For the 3C-calibrator, also write GAINCAL-parset used to SUM stations (core or ears) below
  #if sources[i]==calsource:
  # outfile = Lnum + '_' + sources[i] + '_GAINCAL.parset'
  # with open(outfile, 'w') as f:
  #  f.write('msin =' + Lnum + '_'+sources[i]+'_COMB_AVGD.MS\n')
  #  f.write('msin.datacolumn = DATA\n')
  #  f.write('msout =\n')
  #  f.write('gaincal.baseline=CS*&CS*\n') # From Tammo Jan, only use CS-CS baselines
  #  f.write('gaincal.parmdb ='+Lnum + '_'+calsource+'.parmdb\n')
  #  f.write('gaincal.sourcedb = '+scriptpath+'/'+calsource+'.sourcedb\n') # This assumes a source-db file, which can be made from as skymodel-file using the LOFAR makesourcedb-tool
  #  f.write('gaincal.caltype =phaseonly\n')
  #  f.write('gaincal.solint = 0\n')
  #  f.write('steps = [gaincal]\n')

  # For all sources, write APPLYCAL-parsets to apply ear/core phasing solutions to all stations
  #outfile = Lnum + '_' + sources[i] + '_APPLYCAL.parset'
  #with open(outfile, 'w') as f:
  # f.write('msin ='+Lnum + '_'+sources[i]+'_COMB_AVGD.MS\n')
  # f.write('msin.datacolumn = DATA\n')
  # f.write('msout =.\n')
  # f.write('msout.datacolumn=CORRECTED_DATA\n')
  # f.write('applycal.parmdb ='+calnum + '_'+calsource+'.parmdb.timeext\n') # The "timeext-parmdbs" are made from the .parmdb-files by running e.g. "parmexportcal" to make them time-independent. I found parmexportcal to crash, so I use the "cal_expand" script instead to produce these files
  # #f.write('applycal.timeslotsperparmupdate = 45\n') # To fix bug with many SBs.
  # f.write('steps = [applycal]\n')
  
  # Write parset to sum all core stations
  #outfile = Lnum + '_' + sources[i] + '_SUMCORE.parset'
  #with open(outfile, 'w') as f:
  # f.write('msin ='+Lnum + '_'+sources[i]+'_COMB_AVGD.MS\n')
  # f.write('msin.datacolumn = CORRECTED_DATA\n')
  # f.write('msout ='+Lnum + '_'+sources[i]+'_SUMCORE.MS\n')
  # f.write('adder.type=stationadder\n')
  # f.write('adder.average=True\n')
  # f.write('adder.useweights=True\n')
  # f.write('adder.stations={CORE:[\'CS*\']}\n')
  # f.write('filter.baseline=\'!CS*&&*\'\n')
  # f.write('filter.remove=True\n')
  # f.write('steps = [adder, filter]\n')
  
  # Write file to sum partial stations, e.g. join HBA ears and sum a few core stations.
  # f.write('msin ='+Lnum + '_'+sources[i]+'_COMB_AVGD.MS\n')
  # f.write('msin.datacolumn = CORRECTED_DATA\n')
  # f.write('msout ='+Lnum + '_'+sources[i]+'_SUMPART.MS\n')
  # f.write('adder.type=stationadder\n')
  # f.write('adder.average=True\n')
  # f.write('adder.useweights=True\n')
  # # Station list from 
  # # From https://www.astron.nl/radio-observatory/astronomers/users/technical-information/lofar-array-configuration/lofar-array-conf
  # # And measurement sets
  # ears2join = [1, 11, 13, 17, 21, 24, 26, 30, 31, 32, 101, 201, 301, 302, 401, 501]
  # newsterp = 'STERP:[CS002HBA0,CS002HBA1,CS003HBA0,CS003HBA1,CS004HBA0,CS004HBA1,CS005HBA0,CS005HBA1,CS006HBA0,CS006HBA1,CS007HBA0,CS007HBA1]'
  # newcs = [', CORE'+str(n).rjust(3,'0')+':[CS'+str(n).rjust(3,'0')+'HBA0,CS'+str(n).rjust(3,'0')+'HBA1]' for n in ears2join]
  # f.write('adder.stations={'+newsterp + ''.join(newcs)+'}\n')
  # # Remove also Swedish station from the data because of bad amps
  # f.write('filter.baseline=!CS*&&*;!SE607*&&*\n')
  # f.write('filter.remove=True\n')
  # f.write('steps = [adder, filter]\n')
  
  # Need an extra PARSET file to filter out the summed stations, because of bug in DPPP. (Should be removed in SUMPART step, but were not)
  #outfile = Lnum + '_' + sources[i] + '_FILTERPARTSUMMED.parset'
  #with open(outfile, 'w') as f:
  # f.write('msin ='+Lnum + '_'+sources[i]+'_SUMPART.MS\n')
  # f.write('msin.datacolumn = DATA\n')
  # f.write('msout ='+Lnum + '_'+sources[i]+'_SUMFILT.MS\n')
  # # Remove also Swedish station from the data because of bad amps
  # f.write('filter.baseline=\'!CS*&&*;!SE607*&&*\'\n')
  # f.write('filter.remove=True\n')
  # f.write('steps = [filter]\n')
