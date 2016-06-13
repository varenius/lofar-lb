import os,sys,glob,numpy as np
flist = np.sort(glob.glob('SRM*.tar'))
for f in flist:
    newname = 'L'+f.split('L')[1].split('MS')[0]+'MS.tar'
    os.system('mv '+f+' '+newname)
    os.system('tar xvf '+newname)
    os.system('rm -f '+newname)

