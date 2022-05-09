from srtreduce import *
from srtobservation import *

basedir = '/Volumes/Toshiba5TB2/SRT/19-20/'
outdir = basedir+'n6946_19/'
plotdir = outdir + 'plots/'
toddir = outdir + 'tods/'
mapdir = outdir + 'maps/'

# Define observation parameters
srt = srtreduce(basedir,outdir,band='K')

# First go through the calibrators to analyse those
obs = srtobservation(srt, basedir+'20210205/20210204-231747-19-20-3C84_OPTIMIZED/')

# Then go through the main data files

print(obs.filelist)
exit()
