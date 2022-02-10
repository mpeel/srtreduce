#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from srt_functions import *


basedir = '/Users/mpeel/Documents/maps/srt_perseus/'
basedir = '/Volumes/Toshiba5TB2/SRT/28-19/calibrated/'
outdir = basedir+'tod/'
plotdir = outdir + 'plots/'
toddir = outdir + 'tods/'
mapdir = outdir + 'maps/'
os.makedirs(outdir, exist_ok=True)
os.makedirs(plotdir, exist_ok=True)
os.makedirs(toddir, exist_ok=True)
os.makedirs(mapdir, exist_ok=True)

doplots = False

folderlist = os.listdir(basedir)
todolist = []
for folder in folderlist:
	if '.' not in folder and 'tod' not in folder:
		todolist.append(basedir+folder)
print(todolist)
for inputdir in todolist:
	prefix=inputdir.replace(basedir,'').replace('/','_')
	ext = '.dat'
	numext = 1
	inputlist = os.listdir(inputdir)
	print(inputlist)
	filelist = [f for f in inputlist if ext in f]
	trip = 0
	for filename in filelist:#[0:10]:
		# if 'AZ2' not in filename:
			# continue
		print(inputdir+'/'+filename)
		data = np.loadtxt(inputdir+'/'+filename).T
		ra = data[2]
		dec = data[3]
		data = data[4:]
		# data = data[4+7:-21]
		birdies = [0,1,2,3,4,5,6,7,64,65,66,96,97,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128]
		for birdie in birdies:
			data[birdie-1:birdie][:] = np.nan
		# print(data[0])
		# print(data)
		timestream = np.nanmedian(data,axis=0)
		if doplots:
			plt.plot(data)
			plt.savefig(plotdir+filename+'_bandpass.png')
			plt.clf()
			plt.plot(data.T)
			plt.savefig(plotdir+filename+'_tod.png')
			plt.clf()
			plt.pcolormesh(data.T)
			plt.savefig(plotdir+filename+'_data_waterfall.png')
			plt.clf()
			plt.plot(timestream)
			plt.savefig(plotdir+filename+'_timestream.png')
			plt.clf()
		writefits(toddir+prefix+'_'+filename.replace('.dat','')+'_tod.fits', ra, dec, timestream,filename)
		# exit()
