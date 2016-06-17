#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from lofarpipe.support.data_map import DataMap
from lofarpipe.support.data_map import DataProduct
import numpy as np
import pyrap.tables, math
from correlate import *
from astropy import units as u
from astropy.coordinates import SkyCoord


# Colm Coughlan, March 2016
# used in locapi generic pipeline implementation


def main(infile='', mapfile_dir = '', filename='', wd = '', nP = 3, radius = 2.5):
    infile_map = infile
    outdir = wd
    
    nP = int(nP)
    radius = float(radius)
    
    print('Opened')
		
    tickermap = mapfile_dir+'/'+'TargetListToMapfileTicker'
    try:
        tickerdata = DataMap.load(tickermap)
        tick = int(tickerdata[0].host) + 1 # increment tick from last saved
    except:
        tick = 0	# if a ticker map exists, use it -> otherwise, we are on the first iteration

    ## if tick = 0, need to do the work to make directions files etc., otherwise just update the ticker
    
#    infile = ((DataMap.load(infile_map))[0]).file	# get the actual filename from the map provided

    bigfileid = os.path.join(mapfile_dir, filename)	# this big file holds all the directions
    if tick == 0:
        table = pyrap.tables.table(infile + '/FIELD', readonly = True)
        ra    = math.degrees(float(table.getcol('PHASE_DIR')[0][0][0] ) % (2 * math.pi))
        dec   = math.degrees(float(table.getcol('PHASE_DIR')[0][0][-1]))
        table.close()
        hexcoords = SkyCoord(ra, dec, unit = 'degree', frame='fk5')
        hexcoords = hexcoords.to_string('hmsdms', sep=':')
        hexra = [hexcoords.split(' ')[0] for i in hexcoords]
        hexdec = [hexcoords.split(' ')[1] for i in hexcoords]

        if not os.path.isfile (outdir+'/lobos_stats.sum'):
            os.system ('wget http://www.jb.man.ac.uk/~njj/lobos_stats.sum -P '+outdir)
        lobos = np.loadtxt(outdir+'/lobos_stats.sum',dtype='S')
        for l in lobos:
            newcoords = SkyCoord(l[1],l[2], unit=(u.hourangle, u.deg),frame='fk5')
            new = np.array([newcoords.ra.degree,newcoords.dec.degree])
            try:
                lobos_coord = np.vstack((lobos_coord,new))
            except:
                lobos_coord = np.copy(new)
        a = correlate(np.array([[ra,dec]]),0,1,lobos_coord,0,1,radius)
        map_out = DataMap([])
        name_out = DataMap([])
        for i in np.asarray(a[:,1],dtype='int'):
            if lobos[i][5].count('P')>=nP:
                namera = lobos[i,1].replace(':','').split('.')[0]
                namedec = lobos[i,2].replace(':','').split('.')[0]
                dpppra = lobos[i,1].replace(':','h',1).replace(':','m',1)+'s' ## phase-shift RA
                dpppdec = lobos[i,2].replace(':','d',1).replace(':','m',1)+'s'## phase-shift DEC
                hemisphere = '-' if '-' in hexdec else '+'
                outfile = namera+hemisphere+namedec             ## outfilename
                map_out.data.append(DataProduct('[\"'+dpppra+'\",\"'+dpppdec+'\"]', outfile, False ))
			
        map_out.save(bigfileid)			# save all directions
        current_coords = map_out[0].host	# save current direction
        current_name = map_out[0].file	# save current filename
        n = len(map_out)
        if(n == 0):
            print('Error! No suitable sources found!')
            return({'break' : True, 'error':1})
    else:
        data = DataMap.load(bigfileid)	# load all directions
        current_coords = data[tick].host	# save current direction
        current_name = data[tick].file
        n = len(data)
		
    fileid = os.path.join(mapfile_dir, tickermap)
    tickerdata = DataMap([])
    tickerdata.data.append(DataProduct(tick, current_coords, 0 ))	# make sure to increment tick in tickermap (saving correct value for next run)
    
    if (tick + 1) == n:
        do_break = True
    else:
        do_break = False
    tickerdata.save(fileid)
    result = {'tickermap': fileid, 'targetlist':bigfileid,'cords':current_coords,'cdir':current_name,'ndir':int(n),'tick':tick,'break':do_break}
    return result
