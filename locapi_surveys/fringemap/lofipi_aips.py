from math import *
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList, AIPSMessageLog
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
from scipy import ndimage; from scipy.ndimage import measurements
import matplotlib; from matplotlib import pyplot as plt
import pyfits; from pyfits import getdata,getheader
import re,sys,pickle,numpy as np,os,glob,time,warnings; from numpy import fft

def pfring (aipsname,refant,antennas,source,indisk,delaywin=600,ratewin=20,\
            solint=1,snr=2,logfiledir='./'):
    uvdata = AIPSUVData (aipsname,'FITS',indisk,1)
    fring = AIPSTask ('FRING')
    fring.refant = refant
    fring.indata = uvdata
    fring.calsour[1:] = [source]
    fring.antennas[1:] = antennas
    fring.solint = solint
    fring.aparm[1:] = [0,0,0,0,0,2,snr,0,0]
    fring.dparm[1:] = [0,delaywin,ratewin,0,0,0,0,0,0]
    fring.weightit = 1
    stdout = sys.stdout; sys.stdout = open(logfiledir+aipsname+'.log','a')
    fring.go()
    sys.stdout.close(); sys.stdout = stdout

def pclcal (aipsname,indisk,inver,logfiledir='./'):
    uvdata = AIPSUVData (aipsname,'FITS',indisk,1)
    clcal = AIPSTask ('clcal')
    clcal.indata = uvdata
    clcal.inver = inver
    clcal.snver = inver
    stdout = sys.stdout; sys.stdout = open(logfiledir+aipsname+'.log','a')
    clcal.go()
    sys.stdout.close(); sys.stdout = stdout

def psplit (aipsname,source,indisk,logfiledir='./'):
    uvdata = AIPSUVData (aipsname,'FITS',indisk,1)
    split = AIPSTask ('split')
    split.indata = uvdata
    split.outclass = 'SPLIT'
    split.docalib = 1
    stdout = sys.stdout; sys.stdout = open(logfiledir+aipsname+'.log','a')
    split.go()
    sys.stdout.close(); sys.stdout = stdout
    uvdata = AIPSUVData(source,'SPLIT',indisk,1)
    uvdata.rename(aipsname,'SPLIT',1)

def pload (filename,aipsname,indisk,outcl,logfiledir='./',doindxr=True):
    fitld = AIPSTask ('FITLD')
    fitld.datain = str(filename)
    fitld.outna = aipsname
    fitld.outcl = outcl
    fitld.outdisk = indisk
    stdout = sys.stdout; sys.stdout = open(logfiledir+aipsname+'.log','a')
    fitld.go ()
    if doindxr:
        uvdata = AIPSUVData (aipsname,'FITS',1,1)
        indxr = AIPSTask ('INDXR')
        indxr.cparm[1:] = [0,0,0.5,0,0,0,0,0,0,0]
        indxr.indata = uvdata
        indxr.go()
    sys.stdout.close(); sys.stdout = stdout

def stars (aipsname, incl, indisk, intext='./starsfile',logfiledir='./'):
    stars = AIPSTask('stars')
    stars.inname = aipsname
    stars.inclass = incl
    stars.indisk = indisk
    try:
        stars.stvers = 0    # does not exist in some AIPS versions
    except:
        pass
    stars.intext = './starsfile'
    stdout = sys.stdout; sys.stdout = open(logfiledir+aipsname+'.log','a')
    stars.go()
    sys.stdout.close(); sys.stdout = stdout

def greys (aipsname, incl, indisk, pmin, pmax, stfac, stvers, logfiledir='./'):
    greys = AIPSTask('greys')
    greys.inname = aipsname
    greys.inclass = incl
    greys.indisk = indisk
    greys.pixrange[1:] = [float(pmin),float(pmax)]
    greys.dotv = -1
    greys.stfac = stfac
    try:
        greys.stvers = stvers  # does not exist in some aips versions
    except:
        pass
    stdout = sys.stdout; sys.stdout = open(logfiledir+aipsname+'.log','a')
    greys.go()
    sys.stdout.close(); sys.stdout = stdout

def lwpla (aipsname,incl,indisk,outfile,logfiledir='./'):
    lwpla = AIPSTask('lwpla')
    lwpla.inname = aipsname
    lwpla.inclass = incl
    lwpla.indisk = indisk
    lwpla.outfile = outfile
    stdout = sys.stdout; sys.stdout = open(logfiledir+aipsname+'.log','a')
    lwpla.go()
    sys.stdout.close(); sys.stdout = stdout
