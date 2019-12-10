#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
from scipy.optimize import curve_fit

def linfit(x, A, B):
	return A*x+B

def subtractbaseline(data, option=0, navg=100):
	# Crude baseline removal
	npoints = len(data)
	# print(npoints)
	for i in range(0,npoints//navg):
		start = i*navg
		end = ((i+1)*navg)
		# print('start:' + str(start) + ', end: ' + str(end))
		if option == 0:
			data[start:end] = data[start:end] - np.median(data[start:end])
		else:
			xline = range(0,len(data[start:end]))
			if len(xline) != 0:
				A,B=curve_fit(linfit,xline,data[start:end])[0]
				# print(A,B)
				data[start:end] = data[start:end] - (A*xline+B)

	if option == 0:
		data[(npoints//navg)*navg-1:] = data[(npoints//navg)*navg-1:] - np.median(data[(npoints//navg)*navg-1:])
	else:
		xline = range(0,len(data[(npoints//navg)*navg-1:]))
		if len(xline) != 0:
			A,B=curve_fit(linfit,xline,data[(npoints//navg)*navg-1:])[0]
			# print(A,B)
			data[(npoints//navg)*navg-1:] = data[(npoints//navg)*navg-1:] - (A*xline+B)

	return data


def read_fits_file(filename):

	inputfits = fits.open(filename)
	print(inputfits.info())
	# print(inputfits[0].columns.names)
	print(inputfits[1].columns.names)
	print(inputfits[2].columns.names)
	print(inputfits[3].columns.names)
	print(inputfits[4].columns.names)
	print(inputfits[5].columns.names)
	ra = np.degrees(inputfits[4].data.field(3))
	dec = np.degrees(inputfits[4].data.field(4))
	az = np.degrees(inputfits[4].data.field(1))
	el = np.degrees(inputfits[4].data.field(2))
	# ra = np.degrees(inputfits[4].data.field(1))
	# dec = np.degrees(inputfits[4].data.field(2))
	# az = np.degrees(inputfits[4].data.field(3))
	# el = np.degrees(inputfits[4].data.field(4))
	data = inputfits[4].data.field(-1)
	# exit()
	# cols = inputfits[0].columns
	# col_names = cols.names
	# print(col_names)
	# nmaps = len(cols)
	# maps = []
	# if usehealpixfits:
	# 	maps = hp.read_map(indir+inputfile,field=None)
	# else:
	# 	for i in range(0,nmaps):
	# 		maps.append(inputfits[1].data.field(i))
	# 		print(len(maps[i]))
	# print(len(maps[0]))
	# # Check to see whether we have nested data, and switch to ring if that is the case.
	# if (inputfits[1].header['ORDERING'] == 'NESTED'):
	# 	maps = hp.reorder(maps,n2r=True)
	# newheader = inputfits[1].header.copy(strip=False)
	# inputfits.close()
	return ra, dec, az, el, data

basedir = '/Volumes/Maxtor4TB/SRT/28-19/'
basedir = '/Volumes/Maxtor4TB/SRT/33-18/'
folderlist = os.listdir(basedir)
todolist = []
for folder in folderlist:
	if '.' not in folder:
		subfolderlist = os.listdir(basedir+folder)
		search = 'PERSEUS_AZ'
		# search = '3C84_EL'
		# search = '3C84_OPTIMIZED'
		search = 'W3OH'
		subfolderlist = [f for f in subfolderlist if search in f]
		for subfolder in subfolderlist:
			todolist.append(basedir+folder+'/'+subfolder)
		print(subfolderlist)
print(folderlist)
print(todolist)
# inputdir='/Volumes/Maxtor4TB/SRT/28-19/20191205/20191205-224502-28-19-PERSEUS_AZ'
# prefix='perseus_20191205-224502'
# inputdir='/Volumes/Maxtor4TB/SRT/28-19/20191205/20191205-195834-28-19-PERSEUS_AZ'
# prefix='perseus_20191205-195834'

for inputdir in todolist:
	prefix=inputdir.replace(basedir,'').replace('/','_')
	ext = 'fits0'
	numext = 1
	inputlist = os.listdir(inputdir)
	print(inputlist)
	filelist = [f for f in inputlist if ext in f]
	trip = 0
	for filename in filelist:#[0:10]:
		for i in range(numext):
			az_set, el_set, ra_set, dec_set, data = read_fits_file(inputdir+'/'+filename[:-1]+str(i))
			# Just look at the first 1k channels for now - removing the first 10 channels as bad.
			# data = data[:,10:1023]#2048*2-1]
			data = data[:,8000:8300]
			# Blank the first N channels
			# data[:,1024:1034] = 1.0
			if trip == 0:
				print(np.shape(data))
			data[~np.isfinite(data)] = 0.0
			if trip == 0:
				plt.plot(data.mean(axis=0),'b', label='Mean in given freq bin')
				plt.savefig(prefix+'_bandpass.png')
				plt.clf()
				plt.pcolormesh(data)
				plt.savefig(prefix+'_data_waterfall.png')
			data = data/data.mean(axis=0)
			if trip == 0:
				plt.pcolormesh(data)
				plt.savefig(prefix+'_data_waterfall_bp.png')
				plt.clf()
				print(np.nanmin(data))
				print(np.nanmax(data))
			std = np.nanmedian(data)
			if trip == 0:
				print(std)
			data[np.abs(data) > 5.0*std] = 0.0
			if trip == 0:
				plt.pcolormesh(data)
				plt.savefig(prefix+'_data_waterfall_bp_std.png')
				plt.clf()
			# plt.plot(data)
			# plt.savefig('test_data.png')
			# plt.clf()
			# Now flatten the dataset
			timestream_set = data.mean(axis=1)
			if trip == 0:
				plt.plot(timestream_set)
				plt.savefig(prefix+'_timestream.png')
				plt.clf()
			# exit()
			print(np.shape(az_set))
			if trip == 0:
				az = list(az_set.copy())
				el = list(el_set.copy())
				ra = list(ra_set.copy())
				dec = list(dec_set.copy())
				timestream = list(timestream_set.copy())
				trip = 1
			else:
				az = az + list(az_set)
				el = el + list(el_set)
				ra = ra + list(ra_set)
				dec = dec + list(dec_set)
				timestream = timestream + list(timestream_set)
	print(np.shape(ra))
	# plt.xlim([50,60])
	# plt.ylim([27,36])
	plt.title(prefix)
	plt.plot(ra, dec,'b,')
	plt.xlabel('RA [deg]')
	plt.ylabel('Dec [deg]')
	plt.savefig(prefix+'_ra_dec.png')
	plt.clf()
	# plt.plot(az)
	# plt.ylabel('Azimuth [Deg]')
	# plt.xlabel('Sample')
	# plt.savefig(prefix+'_az.png')
	# plt.clf()
	plt.plot(az, el,'g,')
	plt.xlabel('Azimuth')
	plt.ylabel('Elevation')
	plt.savefig(prefix+'_az_el.png')
	plt.clf()
	plt.plot(timestream)
	plt.savefig(prefix+'_timestream.png')
	plt.clf()
	np.savetxt(prefix+'_radec.txt',np.transpose([ra, dec]))
	print('Hi')


	# Write out the TOD to disk before destriping
	ra_col = fits.Column(name='ra',format='D',array=np.asarray(ra))
	dec_col = fits.Column(name='dec',format='D',array=np.asarray(dec))
	col1 = fits.Column(name='ch1',format='D',array=np.asarray(timestream))
	cols = fits.ColDefs([ra_col, dec_col, col1])
	hdu = fits.BinTableHDU.from_columns(cols)

	hdr = fits.Header()
	hdr['INSTRUM'] = 'SRT-PERS'
	hdr['VERSION'] = '0.0'
	hdr['OBSNAME'] = prefix
	primary_hdu = fits.PrimaryHDU(header=hdr)
	hdul = fits.HDUList([primary_hdu, hdu])
	hdul.writeto(prefix+'_tod.fits',overwrite=True)

# 	# Create a mask of the different azimuth scans
# 	az_mask = np.zeros(len(az)).astype(int)
# 	# Find out whether we're going up or down in elevation
# 	for i in range(1,len(az_mask)):
# 		if az[i-1] < az[i] or az[i-1] > az[i]+1.0:
# 			az_mask[i] = -1
# 		else:
# 			az_mask[i] = 1

# 	# To catch noisy bits in transitions, require that sets of testlen all have to have the same value
# 	testlen = 10
# 	for i in range(0,int(len(az_mask)/testlen)-testlen):
# 		if np.sum(az_mask[i*testlen:i*testlen+testlen]) != testlen and np.sum(az_mask[i*testlen:i*testlen+testlen]) != -testlen:
# 			az_mask[i*testlen:i*testlen+testlen] = 0
# 	# Count how many different sections we have, and update the mask so we can extract them.
# 	num_azscans = 1
# 	notzero = 0
# 	for i in range(0,len(az_mask)):
# 		if az_mask[i] != 0:
# 			az_mask[i] = az_mask[i] * num_azscans
# 			notzero = 1
# 		else:
# 			if notzero == 1:
# 				num_azscans = num_azscans+1
# 				notzero = 0
# 	# print(num_azscans)
# 	plt.plot(az_mask)
# 	plt.savefig(prefix+'_azmask.png')
# 	plt.clf()
# 	az_mask2 = np.abs(az_mask.copy())
# 	plt.plot(az_mask2)
# 	plt.savefig(prefix+'_azmask2.png')
# 	plt.clf()

	timestream = np.asarray(timestream)
# 	for i in range(1,np.max(az_mask2)+1):
# 		# Destripe per azimuth strip
# 		xline=np.arange(0,sum(az_mask2==i))
# 		# print(timestream[az_mask==i])
# 		ret=curve_fit(linfit,xline,timestream[az_mask2==i])
# 		# print(ret)
# 		A=ret[0][0]
# 		B=ret[0][1]
# 		# print(A,B)
# 		timestream[az_mask2==i] = timestream[az_mask2==i] - (A*xline+B)
# 	# Flag the change-over points.
# 	timestream[az_mask==0]=0.0


# # 	timestream = np.asarray(subtractbaseline(timestream))
# 	plt.plot(timestream)
# 	plt.savefig(prefix+'_timestream_sub.png')
# 	plt.clf()

	# Remove glitches
	threshold=5.0
	std = np.std(timestream)
	median = np.median(timestream)
	# print(std)
	timestream[np.abs(timestream-median) > std*threshold] = 0.0
	# Do this twice.
	std = np.std(timestream)
	median = np.median(timestream)
	# print(std)
	timestream[np.abs(timestream-median) > std*threshold] = 0.0
	# timestream[timestream-median < std/threshold] = 1.0
	plt.plot(timestream)
	plt.savefig(prefix+'_timestream_sub_deglitch.png')
	plt.clf()


	numpix = 75
	min_ra = np.min(ra)
	max_ra = np.max(ra)
	pixsize_ra = (max_ra-min_ra)/(numpix-1)
	ras = np.arange(min_ra,max_ra+pixsize_ra,pixsize_ra)
	# print(len(ras))
	# exit()
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
	maparr[weightarr==0] = 1.0
	# print(maparr)
	plt.pcolormesh(ras,decs,maparr.transpose())
	plt.savefig(prefix+'_map.png')
	plt.clf()

	hdr = fits.header.Header()
	fits.writeto(prefix+'_map.fits', maparr, hdr, overwrite=True) 
	fits.writeto(prefix+'_weights_map.fits', weightarr, hdr, overwrite=True) 

	# exit()

