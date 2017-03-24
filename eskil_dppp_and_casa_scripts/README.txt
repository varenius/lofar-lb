README FILE
These files were provided on request by Eskil. They come with very limited
inline comments, and no warranty. But, if you have questions try an email to
eskil.varenius@gmail.com. I have used similar scripts to bring LOFAR data from
the LTA-format (one MS per subband per scan per source) to AIPS UVFITS format
(one UVFITS per source).  Then, in AIPS, I usually run UVDEC (to skip some
channels to reach an even number), MOREIF (to split in IFs for FRING), UVMOD
(to multiply visibilites with 1e-6 for AIPS to not discard them as "infinity"
in some tasks) and WTMOD (to scale weights similarly, usually by around 1e6).
Then, one can run FRING. So, a quick summary is assuming the LTA files are
in the folder "msdir", and the code is in a folder called "code" is:

cd code
Edit the allfiles.dat to reflect the LTA info page for your osbervation
Edit the two .py files to make sure source names, subband numbers and allfiltes.dat column numbers are correct
mkdir ../averaged
cd ../averaged
python writeDPPPparsets.py ../msdir/
DPPP *PRECAL.parset
aoflagger *.MS
(optional DPPP *GAINCAL*, DPPP *APPLYCAL*, DPPP *SUMSTAT* etc.).
mkdir ../UVFITS
cd ../UVFITS
python ../code/writeCASAtoSingleUVFITS.py ../averaged/
start CASA, then run "execfile('BLA_SINGLE_SUMFILT_CONCAT.py')". This will
concatenate files in the ../averaged folder, and give you instructions to
convert to circular. Then it will convert to UVFITS.

Then in AIPS:
FITLD
MSORT
INDXR
UVDEC (should be done in DPPP, is faster )
UVMOD (should be done in DPPP, is faster)
WTMOD (should be done in DPPP, is faster)
USING ONLY LONG (>60klambda) baselines:
FRING
BPASS (assuming some source flux and spectral index, from e.g. NED)
CALIB (A&P)+IMAGR self-cal loop assuming some source flux and spectral index, from e.g. NED
Transfer cumulative to 3C-source (after loading that with FITLD-WTMOD above)
USING ONLY SHORT (<20klambda)baslines:
CALIB (phase only) assumign point
IMAGR for each IF, measuring the flux density. Compare with 3C-values given by AIPS task SETJY
If not flux cal OK, re-iterate above assuming different source flux (spectral index)
If flux cal OK, then transfer all to phase cal source
Run CALIB (phase only) for phase cal source (after loading with FITLD-WTMOD above)
Apply to target (after loading that with FITLD-WTMOD above)
Possibly IMAGR of target, and optionally selfcal. For MS-MFS imaging, export to UVFITS, then import again to MS with CASA, and image with e.g. CLEAN.


AGAIN: No warranty, but ask if something is very strange or if I have forgotten
something. I have scripts for the AIPS-part too, will perhaps add them later.
