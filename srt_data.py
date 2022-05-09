#!/usr/bin/env python
# This class handles reading in SRT data, and making basic modifications to it, such as correcting the pointing
# Note that fixed telescope information is recorded in srtreduce itself - this only records things that change each measurement
# Mike Peel
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
from scipy.optimize import curve_fit

class srtdata():
	def __init__(self,srtreduce, filename=''):
		# Directories
		self.filename = filename
		self.RightAsc = 0.0
		self.Declination = 0.0
		self.ra = []
		self.dec = []
		self.az = []
		self.el = []
		self.mjd = []
		self.pa = []
		self.rotang = []
		self.parang = []
		self.datacube = []
		self.humidity = 0
		self.temperature = 0
		self.pressure = 0
		self.srtreduce = srtreduce
		self.read_srt_fits_file()

	def read_srt_fits_file(self):
		print(self.filename)
		inputfits = fits.open(self.filename)
		# This is what the fits file looks like
		# print(inputfits.info())
		# 0  PRIMARY       1 PrimaryHDU      71   ()
		# 1  SECTION TABLE    1 BinTableHDU     29   7R x 7C   [J, 8A, D, J, D, D, D]
		# 2  RF INPUTS     1 BinTableHDU     34   14R x 9C   [J, J, 8A, D, D, D, D, D, J]
		# 3  FEED TABLE    1 BinTableHDU     21   7R x 4C   [J, D, D, D]
		# 4  DATA TABLE    1 BinTableHDU     40   84R x 11C   [1D, 1D, 1D, 1D, 1D, 1D, 1D, 1D, 1J, 3D, 4096J]
		# 5  ANTENNA TEMP TABLE    1 BinTableHDU     32   83R x 7C   [2D, 2D, 2D, 2D, 2D, 2D, 2D]
		# 6  SERVO TABLE    1 BinTableHDU     38   83R x 9C   [D, D, D, D, D, D, D, D, D]
		# print(inputfits[0].header)

		# The zero'th header contains various values, read those in - mostly through srtreduce
		if not self.srtreduce.is_set:
			self.srtreduce.set_params(inputfits)
		self.RightAsc = inputfits[0].header['RightAscension']
		self.Declination = inputfits[0].header['Declination']

		# The first header contains info about the bins, those are also read in by srtreduce so ignore them here
		# print(inputfits[1].header)
		# print(inputfits[1].columns.names)
		# print(inputfits[1].data)
		# ['id', 'type', 'sampleRate', 'bins', 'frequency', 'bandWidth', 'flux']
		#(7,)
		#[(0, 'stokes', 3000., 1024, 0., 1500., 0.)
	 	# (1, 'stokes', 3000., 1024, 0., 1500., 0.)
	 	# (2, 'stokes', 3000., 1024, 0., 1500., 0.)
	 	# (3, 'stokes', 3000., 1024, 0., 1500., 0.)
	 	# (4, 'stokes', 3000., 1024, 0., 1500., 0.)
	 	# (5, 'stokes', 3000., 1024, 0., 1500., 0.)
	 	# (6, 'stokes', 3000., 1024, 0., 1500., 0.)]

		# The second header gives info about the RF properties of each horn - these are also read in by srtreduce so ignore them here
		# print(inputfits[2].header)
		# print(inputfits[2].columns.names)
		# print(np.shape(inputfits[2].data))
		# print(inputfits[2].data)
		# [(0, 0, 'LCP', 24000., 1400., 23900.,  7., 14.8679312 , 0)
		#  (0, 1, 'RCP', 24000., 1400., 23900.,  7., 12.8185766 , 0)
		#  (1, 0, 'LCP', 24000., 1400., 23900.,  9., 23.78312133, 1)
		#  (1, 1, 'RCP', 24000., 1400., 23900.,  9., 20.31523133, 1)
		#  (2, 0, 'LCP', 24000., 1400., 23900., 11., 23.46308047, 2)
		#  (2, 1, 'RCP', 24000., 1400., 23900., 11., 22.4005672 , 2)
		#  (3, 0, 'LCP', 24000., 1400., 23900.,  7., 22.4089376 , 3)
		#  (3, 1, 'RCP', 24000., 1400., 23900.,  7., 19.67772647, 3)
		#  (4, 0, 'LCP', 24000., 1400., 23900.,  0., 31.35327867, 4)
		#  (4, 1, 'RCP', 24000., 1400., 23900.,  0., 25.12231087, 4)
		#  (5, 0, 'LCP', 24000., 1400., 23900.,  0., 18.44680767, 5)
		#  (5, 1, 'RCP', 24000., 1400., 23900.,  0., 17.39153967, 5)
		#  (6, 0, 'LCP', 24000., 1400., 23900.,  0., 19.688415  , 6)
		#  (6, 1, 'RCP', 24000., 1400., 23900.,  0., 18.88773933, 6)]

		# The third header gives information about the pixel positions - these are also read in by srtreduce so ignore them here
		# print(inputfits[3].header)
		# print(inputfits[3].columns.names)
		# print(np.shape(inputfits[3].data))
		# print(inputfits[3].data)
		# ['id', 'xOffset', 'yOffset', 'relativePower']
		# (7,)
		# [(0,  0.        ,  0.        , 1.  ) (1,  0.00033355, -0.00057773, 0.97)
		#  (2, -0.00033355, -0.00057773, 0.99) (3, -0.0006671 ,  0.        , 0.97)
		#  (4, -0.00033355,  0.00057773, 0.95) (5,  0.00033355,  0.00057773, 0.97)
		#  (6,  0.0006671 ,  0.        , 0.97)]

		# The fourth header has the data ... this is what we want here!
		# ['time', 'raj2000', 'decj2000', 'az', 'el', 'par_angle', 'derot_angle', 'flag_cal', 'flag_track', 'weather', 'Ch0']
		# ... except az/el and ra/dec are the wrong way around, so flip them
		# flag_cal and flag_track seem to be empty
		# data is in 'Ch0'
		# print(inputfits[4].header)
		# print(inputfits[4].columns.names)
		# print(np.shape(inputfits[4].data))
		self.rotang = np.asarray(np.degrees(inputfits[4].data['derot_angle']))
		self.parang = np.asarray(np.degrees(inputfits[4].data['par_angle']))
		self.ra = np.asarray(np.degrees(inputfits[4].data.field(3)))
		self.dec = np.asarray(np.degrees(inputfits[4].data.field(4)))
		self.az = np.asarray(np.degrees(inputfits[4].data.field(1)))
		self.el = np.asarray(np.degrees(inputfits[4].data.field(2)))
		self.humidity = inputfits[4].data['weather'][0]
		self.temperature = inputfits[4].data['weather'][1]
		self.pressure = inputfits[4].data['weather'][2]
		self.data = np.asarray(inputfits[4].data.field(-1), dtype=np.float64)

		# Header 5 is 'ANTENNA TEMP TABLE' - but seems to be empty? So ignore.
		# print(inputfits[5].header)
		# print(inputfits[5].columns.names)
		# print(np.shape(inputfits[5].data))
		# print(inputfits[5].data)

		# Header 6 I don't understand?
		# ['time', 'SRP_TX', 'SRP_TY', 'SRP_TZ', 'SRP_RX', 'SRP_RY', 'SRP_RZ', 'GFR_RZ', 'M3R_RZ']
		# print(inputfits[6].header)
		# print(inputfits[6].columns.names)
		# print(np.shape(inputfits[6].data))
		# print(inputfits[6].data)
		# exit()

		# ... and we're done with reading in!
		inputfits.close()
		return

	def write_srt_fits(filename, prefix):
		ra_col = fits.Column(name='ra',format='D',array=self.ra)
		dec_col = fits.Column(name='dec',format='D',array=self.dec)
		col1 = fits.Column(name='ch1',format='D',array=timestream)
		cols = fits.ColDefs([ra_col, dec_col, col1])
		hdu = fits.BinTableHDU.from_columns(cols)

		hdr = fits.Header()
		hdr['INSTRUM'] = 'SRT-PERS'
		hdr['VERSION'] = '0.0'
		hdr['OBSNAME'] = prefix
		primary_hdu = fits.PrimaryHDU(header=hdr)
		hdul = fits.HDUList([primary_hdu, hdu])
		hdul.writeto(filename,overwrite=True)
		return 0
