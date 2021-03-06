#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
from astropy.wcs import WCS

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
# NGC6946
# ra_cen = 308.718015
# dec_cen = 60.153915
# min_ra = 308.2
# max_ra = 309.2
# min_dec = 59.8
# max_dec = 60.5
# M51
ra_cen = 202.484167
dec_cen = 47.230556
# min_ra = 202.2
# max_ra = 202.8
# min_dec = 47.0
# max_dec = 47.4
# NGC891
# ra_cen = 35.639224
# dec_cen = 42.349146
# res = 0.4 / 60.0
res = 0.5 / 60.0
res_ra = res/np.cos(dec_cen*np.pi/180.0)
# res = 3.0 / 60.0
npix_ra = 60
# npix_ra = int((max_ra-min_ra)/res)
# res_dec = res/np.cos(0.5*(min_dec+max_dec)*np.pi/180.0)
# min_dec = min_dec/np.cos(min_dec*np.pi/180.0)
# max_dec = max_dec/np.cos(max_dec*np.pi/180.0)
# npix_dec = int((max_dec-min_dec)/(res_dec))
# npix_dec = int((max_dec-min_dec)/(res))
npix_dec = npix_ra
min_ra = ra_cen - 0.5*npix_ra*res_ra
min_dec = dec_cen - 0.5*npix_dec*res
max_ra = ra_cen + 0.5*npix_ra*res_ra
max_dec = dec_cen + 0.5*npix_dec*res
print('Num in RA:' + str(npix_ra))
print('Num in Dec:' + str(npix_dec))
print('Min RA:' + str(min_ra))
print('Min Dec:' + str(min_dec))
pixsize_ra = res_ra#(max_ra-min_ra)/(numpix-1)
ras = np.arange(min_ra,max_ra+pixsize_ra,pixsize_ra)
print(len(ras))
# exit()
# min_dec = np.min(dec)
# max_dec = np.max(dec)
pixsize_dec = res#(max_dec-min_dec)/(numpix-1)
decs = np.arange(min_dec,max_dec+pixsize_dec,pixsize_dec)
print(len(decs))
wcs_list = WCS(naxis=2)
wcs_list.wcs.crpix = [npix_ra/2, npix_ra/2]
wcs_list.wcs.crval = [ra_cen, dec_cen]
wcs_list.wcs.cunit = ["deg", "deg"]
wcs_list.wcs.ctype = ["RA---TAN", "DEC--TAN"]
wcs_list.wcs.cdelt = [res_ra, res]
wcs_list.array_shape = [npix_ra, npix_dec]
fits_header = wcs_list.to_header()
print(fits_header)

# Set 8 to do a combined one at the maxend-end
# numhorns = 8
numhorns = 1
for h in range(0,numhorns):

	maparr = np.zeros([npix_ra,npix_dec])
	weightarr = np.zeros([npix_ra,npix_dec])
	hitarr = np.zeros([npix_ra,npix_dec])
	# hdr = fits.header.Header()
	# fits.writeto('test_image.fits', maparr, hdr, overwrite=True) 
	# exit()

	# if h == 7:
	ext = 'fits'
	# else:
	# ext = 'fits'+str(h)
	numext = 1
	# basedir = '/Volumes/Toshiba5TB2/SRT/28-19/reduce/'
	# basedir = '/Volumes/Toshiba5TB2/SRT/28-19/calibrated/tod/'
	# basedir = '/Volumes/Toshiba5TB2/SRT/19-20/ngc6946_19/'
	# basedir = '/Volumes/Toshiba5TB2/SRT/19-20/ngc6946_24/'
	# basedir = '/Volumes/Toshiba5TB2/SRT/19-20/n891_24/'
	# basedir = '/Volumes/Toshiba5TB2/SRT/19-20/n891_19/'
	# basedir = '/Volumes/Toshiba5TB2/SRT/19-20/m51_19/'
	basedir = '/Volumes/Toshiba5TB2/SRT/19-20/m51_24/'
	inputlist = os.listdir(basedir+'tods/')
	filelist = [f for f in inputlist if ext in f]
	# print(filelist)
	filelist.sort()
	# print(filelist)
	i = 0


	prevfile = ''
	for file in filelist:
		print(file)
		# if prevfile != '':
		# 	break
		# prevfile = file
		# if 'DEC' in file:
		# 	break
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

			timestream = np.asarray(subtractbaseline(timestream))
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
			# print(len(timestream))
			# print(ra)
			# print(dec)
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
				if pix_ra < npix_ra and pix_dec < npix_dec and timestream[i] != 0:
					maparr[pix_ra,pix_dec] += timestream[i]/variance
					weightarr[pix_ra,pix_dec] += 1.0/variance
					hitarr[pix_ra,pix_dec] += 1

	maparr[weightarr != 0] /= weightarr[weightarr != 0]
	maparr[hitarr == 0] = np.nan
	weightarr[hitarr == 0] = np.nan
	hitarr[hitarr == 0] = np.nan
	# Write out maps
	# hdr = fits.PrimaryHDU(header=fits_header)
	# hdr = fits.header.Header()
	hdr = fits_header
	print(hdr)
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

