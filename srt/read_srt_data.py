#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
from scipy.optimize import curve_fit
from srt_functions import *


# basedir = '/Volumes/Toshiba5TB2/SRT/28-19/modified/K-Band/'
# basedir = '/Volumes/Maxtor4TB/SRT/33-18/'
# basedir = '/Volumes/Toshiba5TB2/SRT/19-20/20201124/'
basedir = '/Volumes/Toshiba5TB2/SRT/19-20/20201125/'
basedir = '/Volumes/Toshiba5TB2/SRT/19-20/20201127/'
# outdir = '/Volumes/Toshiba5TB2/SRT/28-19/reduce/'
# outdir = basedir+'../ngc6946_19/'
outdir = basedir+'../ngc6946_24/'
outdir = basedir+'../m51_19/'
plotdir = outdir + 'plots/'
toddir = outdir + 'tods/'
mapdir = outdir + 'maps/'
os.makedirs(outdir, exist_ok=True)
os.makedirs(plotdir, exist_ok=True)
os.makedirs(toddir, exist_ok=True)
os.makedirs(mapdir, exist_ok=True)

doplots = False

birdies = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 107, 108, 109, 110, 111, 128, 256, 384, 512, 513, 640, 768, 896, 972, 1023, 1024]
numchans = 1024
folderlist = os.listdir(basedir)
todolist = []
# This is for where we have subfolders
# for folder in folderlist:
# 	if '.' not in folder:
# 		subfolderlist = os.listdir(basedir+folder)
# 		search = 'PERSEUS'
# 		# search = 'PERSEUS_AZ'
# 		# search = '3C84_EL'
# 		# search = '3C84_OPTIMIZED'
# 		# search = 'W3OH'
# 		subfolderlist = [f for f in subfolderlist if search in f]
# 		for subfolder in subfolderlist:
# 			todolist.append(basedir+folder+'/'+subfolder)
# 		print(subfolderlist)
# This is for where we only have one folder
for folder in folderlist:
	# if 'N6946' in folder:
	if 'M51' in folder:
		todolist.append(basedir+folder)
print(folderlist)
print(todolist)
for inputdir in todolist:
	prefix=inputdir.replace(basedir,'').replace('/','_')
	ext = 'fits'
	numext = 1
	inputlist = os.listdir(inputdir)
	print(inputlist)
	filelist = [f for f in inputlist if ext in f]
	trip = 0
	for filename in filelist:#[0:10]:
		if 'summary' in filename:
			continue
		for i in range(numext):
			az_set, el_set, ra_set, dec_set, data = read_fits_file(inputdir+'/'+filename[:])#+str(i))
			horn_number = int(filename[-1])
			nanarr = np.zeros(len(az_set))
			nanarr[:] = np.nan
			# Just look at the first 1k channels for now - removing the first 10 channels as bad.
			# Ordering is LL, RR, U, Q
			# Only using LL, RR for now
			data = data[:,0:2048]
			print(np.shape(data))
			# data = np.asarray(data[:,10:1023])#2048*2-1]
			# data = data[:,8000:8300]
			# Blank the first N channels
			# data[:,1024:1034] = 1.0
			# print(np.nanargmax(np.nanmean(data,axis=0)))
			for birdie in birdies:
				data[:,birdie-1:birdie][:] = np.nan
				data[:,numchans+birdie-1:numchans+birdie] = np.nan
				# data[:,2*numchans+birdie-1:2*numchans+birdie] = 0.0#np.nan
				# data[:,3*numchans+birdie-1:3*numchans+birdie] = 0.0#np.nan
			# if horn_number == 2:
			data[:,0:60] = np.nan
			data[:,1024:1024+60] = np.nan
			# Horn 5 LL is missing
			if horn_number == 5:
				data[:,0:1024] = np.nan
			# Horns 3 and 6 RR is weird
			if horn_number == 3 or horn_number == 6:
				data[:,1024:] = np.nan
			# data[:,190:192] = 1.0#np.nan
			# data[:,226:308] = 1.0#np.nan
			# data[:,410:483] = 1.0#np.nan
			# data[:,505:534] = 1.0#np.nan
			# data[:,580:580] = 1.0#np.nan
			# data[:,715:715] = 1.0#np.nan
			# data[:,761:778] = 1.0#np.nan
			# data[:,812:814] = 1.0#np.nan
			if trip == 0:
				print(np.shape(data))
			# data[~np.isfinite(data)] = 0.0
			# if trip == 0:
			# 	plt.ylim(np.nanmin(data),np.nanmax(data))
			# 	plt.plot(np.nanmean(data,axis=0),'b', label='Mean in given freq bin')
			# 	plt.savefig(plotdir+filename+'_bandpass.png')
			# 	plt.clf()
			# 	plt.pcolormesh(data)
			# 	plt.savefig(plotdir+prefix+'_data_waterfall.png')
			# 	plt.clf()
			# plotdata = data.copy()
			# plotdata[~np.isfinite(plotdata)] = 0.0
			# print(np.nanmin(plotdata))
			if doplots:
				plt.plot(np.nanmean(data,axis=0),'b.', label='Mean in given freq bin')
				plt.savefig(plotdir+filename+'_bandpass.png')
				plt.clf()
				plt.pcolormesh(data)
				plt.savefig(plotdir+filename+'_data_waterfall.png')
				plt.clf()

			# # Remove bad values in the bandpass
			# threshold=3.0
			# bandpass = np.nanmean(data,axis=0)
			# # for test in range(0,len(bandpass)):
			# # 	print(str(test) + ' ' + str(bandpass[test]))
			# std = np.nanstd(bandpass)
			# median = np.nanmedian(bandpass)
			# # print(std)
			# bandpass[np.abs(bandpass-median) > std*threshold] = np.nan
			# # Do this twice.
			# std = np.nanstd(bandpass)
			# median = np.nanmedian(bandpass)
			# # print(std)
			# bandpass[np.abs(bandpass-median) > std*threshold] = np.nan
			# # timestream[timestream-median < std/threshold] = 1.0
			# if trip == 0:
			# 	plt.plot(bandpass)
			# 	plt.savefig(plotdir+prefix+'_bandpass_deglitch.png')
			# 	plt.clf()
			# # Flag based on differencing as well
			# diff = np.abs(np.diff(bandpass,append=1))
			# # diff = np.concatenate(np.asarray(diff), np.asarray([1000]))
			# bandpass[diff>5] = np.nan
			# for j in range(0,len(bandpass)-1):
			# 	if bandpass[j+1] == np.nan:
			# 		bandpass[j] = np.nan

			# if trip == 0:
			# 	plt.plot(bandpass)
			# 	plt.savefig(plotdir+prefix+'_bandpass_deglitch2.png')
			# 	plt.clf()
			# badvals = np.isnan(bandpass)
			# print(badvals==True)
			# print(np.sum(badvals==True))
			# data[:,badvals==True] = data[:,badvals==True]*np.nan

			data = data/np.nanmean(data,axis=0)
			maxval = np.nanargmax(np.nanstd(data,axis=0))
			if maxval > 1024:
				maxval -= 1024
			print(maxval)
			if trip == 0 and doplots:
				plt.pcolormesh(data)
				plt.savefig(plotdir+filename+'_data_waterfall_bp.png')
				plt.clf()
				print(np.nanmin(data))
				print(np.nanmax(data))
			std = np.nanmedian(data)
			if trip == 0:
				print(std)
			data[np.abs(data) > 5.0*std] = 0.0
			if doplots:
				plt.pcolormesh(data)
				plt.savefig(plotdir+filename+'_data_waterfall_bp_std.png')
				plt.clf()
				plt.plot(np.nanstd(data,axis=0),'b.', label='Mean in given freq bin')
				plt.savefig(plotdir+filename+'_bandpass_flat.png')
				plt.clf()
			# plt.plot(data)
			# plt.savefig('test_data.png')
			# plt.clf()
			# Now flatten the dataset
			timestream_set = np.nanmean(data,axis=1)
			# And remove a linear trend
			timestream_set = subtractbaseline(timestream_set,option=1,navg=len(timestream_set))
			# if trip == 0:
			if doplots:
				plt.plot(timestream_set)
				plt.savefig(plotdir+filename+'_timestream.png')
				plt.clf()
				print(np.shape(az_set))

			# Write out the TOD to disk
			writefits(toddir+filename.replace('.','_')+'_tod.fits', ra_set, dec_set, timestream_set,prefix)

			if trip == 0:
				# az = list(az_set.copy())
				# el = list(el_set.copy())
				ra = list(ra_set.copy())
				dec = list(dec_set.copy())
				timestream = list(timestream_set.copy())
				trip = 1
			else:
				# az = az + list(az_set)
				# el = el + list(el_set)
				ra = ra + list(ra_set)
				dec = dec + list(dec_set)
				timestream = timestream + list(timestream_set)
		# 	break
		# break
	print(np.shape(ra))
	# plt.xlim([50,60])
	# plt.ylim([27,36])
	plt.title(prefix)
	plt.plot(ra, dec,'b,')
	plt.xlabel('RA [deg]')
	plt.ylabel('Dec [deg]')
	plt.savefig(plotdir+prefix+'_ra_dec.png')
	plt.clf()
	# plt.plot(az)
	# plt.ylabel('Azimuth [Deg]')
	# plt.xlabel('Sample')
	# plt.savefig(prefix+'_az.png')
	# plt.clf()
	# plt.plot(az, el,'g,')
	# plt.xlabel('Azimuth')
	# plt.ylabel('Elevation')
	# plt.savefig(plotdir+prefix+'_az_el.png')
	# plt.clf()
	# plt.plot(timestream)
	# plt.savefig(plotdir+prefix+'_timestream.png')
	# plt.clf()
	# np.savetxt(plotdir+prefix+'_radec.txt',np.transpose([ra, dec]))
	# print('Hi')

	# # Create a mask of the different azimuth scans
	# az_mask = np.zeros(len(az)).astype(int)
	# # Find out whether we're going up or down in elevation
	# for i in range(1,len(az_mask)):
	# 	if az[i-1] < az[i] or az[i-1] > az[i]+1.0:
	# 		az_mask[i] = -1
	# 	else:
	# 		az_mask[i] = 1

	# # To catch noisy bits in transitions, require that sets of testlen all have to have the same value
	# testlen = 10
	# for i in range(0,int(len(az_mask)/testlen)-testlen):
	# 	if np.sum(az_mask[i*testlen:i*testlen+testlen]) != testlen and np.sum(az_mask[i*testlen:i*testlen+testlen]) != -testlen:
	# 		az_mask[i*testlen:i*testlen+testlen] = 0
	# # Count how many different sections we have, and update the mask so we can extract them.
	# num_azscans = 1
	# notzero = 0
	# for i in range(0,len(az_mask)):
	# 	if az_mask[i] != 0:
	# 		az_mask[i] = az_mask[i] * num_azscans
	# 		notzero = 1
	# 	else:
	# 		if notzero == 1:
	# 			num_azscans = num_azscans+1
	# 			notzero = 0
	# # print(num_azscans)
	# plt.plot(az_mask)
	# plt.savefig(prefix+'_azmask.png')
	# plt.clf()
	# az_mask2 = np.abs(az_mask.copy())
	# plt.plot(az_mask2)
	# plt.savefig(prefix+'_azmask2.png')
	# plt.clf()

	# timestream = np.asarray(timestream)
	# for i in range(1,np.max(az_mask2)+1):
	# 	# Destripe per azimuth strip
	# 	xline=np.arange(0,sum(az_mask2==i))
	# 	# print(timestream[az_mask==i])
	# 	ret=curve_fit(linfit,xline,timestream[az_mask2==i])
	# 	# print(ret)
	# 	A=ret[0][0]
	# 	B=ret[0][1]
	# 	# print(A,B)
	# 	timestream[az_mask2==i] = timestream[az_mask2==i] - (A*xline+B)
	# # Flag the change-over points.
	# timestream[az_mask==0]=0.0


	# # timestream = np.asarray(subtractbaseline(timestream))
	# plt.plot(timestream)
	# plt.savefig(prefix+'_timestream_sub.png')
	# plt.clf()

	# Remove glitches
	threshold=5.0
	timestream = np.asarray(timestream, dtype=np.float64)
	std = np.nanstd(timestream)
	median = np.nanmedian(timestream)
	# print(std)
	timestream[np.abs(timestream-median) > std*threshold] = np.nan
	# Do this twice.
	std = np.nanstd(timestream)
	median = np.nanmedian(timestream)
	# print(std)
	try:
		timestream[np.abs(timestream-median) > std*threshold] = np.nan
	except:
		continue
	# timestream[timestream-median < std/threshold] = 1.0
	plt.plot(timestream)
	plt.savefig(plotdir+prefix+'_timestream_sub_deglitch.png')
	plt.clf()


	numpix = 60
	min_ra = np.min(ra)
	max_ra = np.max(ra)
	pixsize_ra = (max_ra-min_ra)/(numpix-1)
	ras = np.arange(min_ra,max_ra+pixsize_ra,pixsize_ra)
	# print(len(ras))
	min_dec = np.min(dec)
	max_dec = np.max(dec)
	pixsize_dec = (max_dec-min_dec)/(numpix-1)
	decs = np.arange(min_dec,max_dec+pixsize_dec,pixsize_dec)
	maparr = np.zeros([numpix,numpix])
	weightarr = np.zeros([numpix,numpix])
	for i in range(len(timestream)):
		pix_ra = int((ra[i]-min_ra)/pixsize_ra)
		pix_dec = int((dec[i]-min_dec)/pixsize_dec)
		# print(pix_ra,pix_dec)
		maparr[pix_ra,pix_dec] += timestream[i]
		weightarr[pix_ra,pix_dec] += 1

	maparr[weightarr != 0] /= weightarr[weightarr != 0]
	maparr[weightarr==0] = np.nan
	# print(maparr)
	plt.pcolormesh(ras,decs,maparr.transpose())
	plt.savefig(mapdir+prefix+'_map.png')
	plt.clf()

	hdr = fits.header.Header()
	fits.writeto(mapdir+prefix+'_map.fits', maparr, hdr, overwrite=True) 
	fits.writeto(mapdir+prefix+'_weights_map.fits', weightarr, hdr, overwrite=True) 
