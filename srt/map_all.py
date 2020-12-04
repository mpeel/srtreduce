#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os

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


# Set up the map
# min_ra = 49.5
# max_ra = 60.0
# min_dec = 27.5
# max_dec = 36.0
min_ra = 308.2
max_ra = 309.2
min_dec = 59.8
max_dec = 60.5
min_ra = 202.1
max_ra = 203.8
min_dec = 47.0
max_dec = 47.5
# res = 0.4 / 60.0
res = 1.0 / 60.0
# res = 3.0 / 60.0
npix_ra = int((max_ra-min_ra)/res)
npix_dec = int((max_dec-min_dec)/res)
print('Num in RA:' + str(npix_ra))
print('Num in Dec:' + str(npix_dec))
pixsize_ra = res#(max_ra-min_ra)/(numpix-1)
ras = np.arange(min_ra,max_ra+pixsize_ra,pixsize_ra)
print(len(ras))
# exit()
# min_dec = np.min(dec)
# max_dec = np.max(dec)
pixsize_dec = res#(max_dec-min_dec)/(numpix-1)
decs = np.arange(min_dec,max_dec+pixsize_dec,pixsize_dec)
print(len(decs))
# Set 8 to do a combined one at the maxend-end
numhorns = 8
# numhorns = 1
for h in range(0,numhorns):

	maparr = np.zeros([npix_ra,npix_dec])
	weightarr = np.zeros([npix_ra,npix_dec])
	hitarr = np.zeros([npix_ra,npix_dec])
	# hdr = fits.header.Header()
	# fits.writeto('test_image.fits', maparr, hdr, overwrite=True) 
	# exit()

	# if h == 7:
	ext = '.fits'
	# else:
		# ext = 'fits'+str(h)
	numext = 1
	basedir = '/Volumes/Toshiba5TB2/SRT/28-19/reduce/'
	basedir = '/Volumes/Toshiba5TB2/SRT/28-19/calibrated/tod/'
	basedir = '/Volumes/Toshiba5TB2/SRT/19-20/ngc6946_19/'
	basedir = '/Volumes/Toshiba5TB2/SRT/19-20/ngc6946_24/'
	basedir = '/Volumes/Toshiba5TB2/SRT/19-20/m51_19/'
	inputlist = os.listdir(basedir+'tods/')
	filelist = [f for f in inputlist if ext in f]
	# print(filelist)
	filelist.sort()
	print(filelist)
	i = 0



	for file in filelist:
		print(file)
		# if '018' in file:
			# break
		if 'map.fits' not in file:
			inputfits = fits.open(basedir+'tods/'+file)
			cols = inputfits[1].columns
			# print(cols.info())

			ra = inputfits[1].data.field(0)
			dec = inputfits[1].data.field(1)
			timestream = inputfits[1].data.field(2)
			# print(np.median(ra)*60.0)
			# print(dec)
			# print(timestream)
			# print(len(ra))

			# timestream = np.asarray(subtractbaseline(timestream))
			# plt.plot(timestream)
			# plt.savefig(prefix+'_timestream_sub.png')
			# plt.clf()

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
			# plt.plot(timestream)
			# plt.savefig(prefix+'_timestream_sub_deglitch.png')
			# plt.clf()

			minstart=0
			maxend=len(timestream-1)
			print(len(timestream))
			for i in range(len(timestream)):
				# start=i-15
				# end=i+15
				# if start < minstart:
				# 	end += minstart-start
				# 	start=minstart
				# if end > maxend:
				# 	start -= maxend-end
				# 	end=maxend
				# print(start)
				# print(end)
				# start = 0
				# end = len(timestream)
				# variance = np.var(timestream[start:end])
				variance = np.var(timestream)
				if variance == 0:
					variance = 1.0
				if ~np.isfinite(variance):
					variance=1.0
				pix_ra = int((ra[i]-min_ra)/pixsize_ra)
				pix_dec = int((dec[i]-min_dec)/pixsize_dec)
				# print(pix_ra,pix_dec)
				maparr[pix_ra,pix_dec] += timestream[i]/variance
				weightarr[pix_ra,pix_dec] += 1.0/variance
				hitarr[pix_ra,pix_dec] += 1

	maparr[weightarr != 0] /= weightarr[weightarr != 0]
	maparr[hitarr == 0] = np.nan
	weightarr[hitarr == 0] = np.nan
	hitarr[hitarr == 0] = np.nan
	# Write out maps
	hdr = fits.header.Header()
	fits.writeto(basedir+'combined_map'+str(h)+'.fits', maparr, hdr, overwrite=True) 
	fits.writeto(basedir+'combined_weights_map'+str(h)+'.fits', weightarr, hdr, overwrite=True) 
	fits.writeto(basedir+'hitmap'+str(h)+'.fits', hitarr, hdr, overwrite=True) 


	# maparr[weightarr==0] = 1.0
	# print(maparr)
	try:
		plt.pcolormesh(ras,decs,maparr.transpose())
		plt.savefig(basedir+'combined_map_img'+str(h)+'.png')
		plt.clf()
		plt.pcolormesh(ras,decs,weightarr.transpose())
		plt.savefig(basedir+'combined_weights_img'+str(h)+'.png')
		plt.clf()
		plt.pcolormesh(ras,decs,hitarr.transpose())
		plt.savefig(basedir+'hitmap_img'+str(h)+'.png')
		plt.clf()
	except:
		continue

