# Script to calculate smearing losses from Bridle and Schwabb (1999) formulas.
# Will by default create a LATEX table used in the Long Baseline chapter of the LOFAR cookbook.
# Author: Eskil Varenius (eskil.varenius@chalmers.se) 2014-10-14
import numpy as np
from scipy.constants import c
from scipy.optimize import fsolve
import scipy.special as spec

# Define the loss, i.e. the fractional amplitude that we allow. 
# For 5% amplitude decrease, we set loss=0.85
loss = 0.95

# The time averaging loss formula by B&S 1999, their formula 18-43
# Here it is defined as a function to minimize with Scipy, i.e.
# it will return 0 if the loss factor equals "loss" defined above.
def bs_time(offset_rad):
    tl = 1- 1.22e-9*(offset_rad/beam_rad)**2*tavg**2
    return loss - tl

# The time averaging loss formula by B&S 1999, their formula 18-24
# Here it is defined as a function to minimize with Scipy, i.e.
# it will return 0 if the loss factor equals "loss" defined above.
def bs_freq (offset_rad):
    fl = np.sqrt(np.pi)/(2*np.sqrt(np.log(2)))*beam_rad*nu/(offset_rad*bw)*spec.erf(np.sqrt(np.log(2))*(offset_rad*bw)/(beam_rad*nu))
    return loss-fl

# Set averaging time, baseline length (defining resolution) 
# and bandwidth (per channel)
tavg = 1 # seconds
b = 1000.0e3 # m, longest baseline
bw = 0.1953125*1e6/64.0 # Bandwidth, Hz, assuming 195kHz bandwidth per subband

# List all frequencies used in App B by Van Haarlem et al 2013.
nus = np.array([15,30,45,60,75,120,150,180,200,210,240])*1e6 # Hz

# Initialize lists to hold results
tloss = [] # loss factors from time averaging
floss = [] # loss factors from frequency averaging
res = [] # Beam sizes, FWHM
lams = [] # Wavelengths, m
for nu in nus:
    lam = c/nu # m
    lams.append(str(round(lam,2)))
    # Calculate beam size using lambda/b
    beam_rad = 0.8*lam/b # from van Haarlem 2013, last page.
    res.append(str(round(beam_rad*3600.0*180/np.pi,2)))
    # Set initual guess for smearing radius
    tguess = 1.0 * np.pi/180.0 # Initial guess 
    fguess = 1.0 * np.pi/180.0 # Initial guess 

    # Use numerical optimizer to find solutions when time and freq loss equals "loss".
    tlossdist_rad = fsolve(bs_time, tguess)[0]
    tlossdist_deg = tlossdist_rad*180.0/np.pi
    flossdist_rad = fsolve(bs_freq, fguess)[0]
    flossdist_deg = flossdist_rad*180.0/np.pi
    
    # Save values to array, multiplied by two to get FoV
    tloss.append(str(round(2*tlossdist_deg,2)))
    floss.append(str(round(2*flossdist_deg,2)))

# Station FoV to print table in the end
stat_fov = [19.39, 9.70, 6.46, 4.85, 3.88, 2.59, 2.07, 1.73, 1.55, 1.48, 1.29]

# Print LATEX table
for i, nu in enumerate(nus):
    print str(int(1e-6*nu)) + ' & ' + lams[i].ljust(4,'0') + ' & ' + res[i].ljust(4,'0') + ' & ' + str(stat_fov[i]).ljust(4,'0') + ' & ' + tloss[i].ljust(4,'0') + ' & ' + floss[i].ljust(4,'0') + '\\\\'
