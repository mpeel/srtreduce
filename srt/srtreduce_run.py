from srtreduce import *
from srtdata import *
from srtobservation import *

basedir = '/Users/mpeel/Desktop/20210204-232112-19-20-3C84_OPTIMIZED/'
outdir = basedir+'../test_calibrator/'
plotdir = outdir + 'plots/'
toddir = outdir + 'tods/'
mapdir = outdir + 'maps/'

srt = srtreduce(basedir,outdir)
print(srt.rotang)
obs = srtobservation(basedir)
print(obs.filelist)
exit()
datafile = srtdata(srt, basedir+'20210204-232112-19-20-3C84_OPTIMIZED_001_001.fits0')
