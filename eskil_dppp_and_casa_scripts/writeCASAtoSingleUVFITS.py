import sys, os
import glob

msdir = sys.argv[1] # Directory containing all the gazillion ms files from the LOFAR LTA, i.e. one file per source per scan per subband.

objects = {'ARP299':[], '1127+5841':[], '1128+5925':[], '3C295':[]}

scriptpath = os.path.dirname(os.path.realpath(sys.argv[0]))
data = open(scriptpath+'/allfiles.dat') # File with info from LTA average pipelinesummary page

lnums = []
ids = []
sources = []

for line in data:
 if not line.startswith('#'):
  # Hard coded indices reflect coumns in the LTA-infopage file. Edit as needed.
  items = line.split()
  ids.append(items[0])
  sources.append(items[3].split('/')[0])
  lnums.append(items[5])

for obj,val in objects.iteritems():
 for i, lnum in enumerate(lnums):
  if obj in sources[i]:
   objects[obj].append(lnum)

files = glob.glob(msdir + '*.MS')
for obj,val in objects.iteritems():
 outfile = obj + '_SINGLE_CONCAT.py'
 with open(outfile, 'w') as f:
  f.write('import numpy as np\n')
  f.write('msin = [\n')
  for lnum in val:
   fname = msdir + lnum + '_' + obj+'_SUMFILT.MS' # Define name of input MS, i.e. the names of the "one MS per source per scan"-files
   if not fname in files:
    print 'FILE ' +fname + ' NOT FOUND!'
   f.write('\'' + fname + '\'')
   if lnum == val[-1]:
    f.write(']\n')
   else:
    f.write(',\n')
  f.write('outms = \'' + obj + '_SINGLE_CONCAT.MS\'\n')
  f.write('outfits = \'' + obj + '_SINGLE_CONCAT.UVFITS\'\n\n')
  f.write('#Concatenate all files to one single MS.\n')
  f.write('#os.system(\'rm -rf \' + outms)\n')
  f.write('concat(vis=msin, concatvis = outms)\n\n')
  f.write('# Change source name to sensible instead of BEAM_X\n')
  f.write("tb.open(os.path.join(outms,'FIELD'),nomodify=False)\n")
  f.write("oldname = tb.getcol('NAME')\n")
  f.write('# Change source name, but truncated to 8 chars for some AIPS tasks\n')
  f.write('newname= np.array([\'' + obj[0:8] + '\'])\n')
  f.write('print \'Changing fieldname from \' + oldname[0] + \' to \' + newname[0] + \'.\'\n')
  f.write("tb.putcol('NAME',newname)\n\n")
  f.write('print \'Changing from linear to circular polarisation...\'\n')
  f.write('# Change polarisation from linear to circular.\n')
  f.write('raw_input(\"PLEASE RUN THE COMMAND: taql \\"update \" + outms + \" set DATA=mscal.stokes(DATA,\'circ\')\\"\")\n')
  f.write('raw_input(\"PLEASE RUN THE COMMAND: taql \\"update \" + outms + \"/POLARIZATION set CORR_TYPE=[5,6,7,8]\\"\")\n')
  f.write('print \'...done!\'\n\n')
  f.write('print \'Exporting to UVFITS for AIPS import...\'\n')
  f.write('exportuvfits(vis=outms, fitsfile = outfits)\n')
  f.write('print \'...done! Ready with CASA pre-processing, run AIPS FRING etc.\'\n')
