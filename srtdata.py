#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
from scipy.optimize import curve_fit

class srtdata():
	def __init__(self,srtreduce, filename=''):
		# Directories
		self.filename = filename
		self.ra = []
		self.dec = []
		self.az = []
		self.el = []
		self.mjd = []
		self.pa = []
		self.rotang = []
		self.datacube = []
		self.srtreduce = srtreduce
		self.read_srt_fits_file()

	def read_srt_fits_file(self):
		print(self.filename)
		inputfits = fits.open(self.filename)

		if not self.srtreduce.is_set:
			self.srtreduce.set_params(inputfits)
		# exit()
		print(inputfits.info())
		print(inputfits[0].header)
		# 0  PRIMARY       1 PrimaryHDU      71   ()
		# 1  SECTION TABLE    1 BinTableHDU     29   7R x 7C   [J, 8A, D, J, D, D, D]
		# 2  RF INPUTS     1 BinTableHDU     34   14R x 9C   [J, J, 8A, D, D, D, D, D, J]
		# 3  FEED TABLE    1 BinTableHDU     21   7R x 4C   [J, D, D, D]
		# 4  DATA TABLE    1 BinTableHDU     40   84R x 11C   [1D, 1D, 1D, 1D, 1D, 1D, 1D, 1D, 1J, 3D, 4096J]
		# 5  ANTENNA TEMP TABLE    1 BinTableHDU     32   83R x 7C   [2D, 2D, 2D, 2D, 2D, 2D, 2D]
		# 6  SERVO TABLE    1 BinTableHDU     38   83R x 9C   [D, D, D, D, D, D, D, D, D]

		# print(inputfits[1].header)
		print(inputfits[1].columns.names)
		print(np.shape(inputfits[1].data[0]))
		print(inputfits[1].data)
		# ['id', 'type', 'sampleRate', 'bins', 'frequency', 'bandWidth', 'flux']
		#(7,)
		#[(0, 'stokes', 3000., 1024, 0., 1500., 0.)
	 	# (1, 'stokes', 3000., 1024, 0., 1500., 0.)
	 	# (2, 'stokes', 3000., 1024, 0., 1500., 0.)
	 	# (3, 'stokes', 3000., 1024, 0., 1500., 0.)
	 	# (4, 'stokes', 3000., 1024, 0., 1500., 0.)
	 	# (5, 'stokes', 3000., 1024, 0., 1500., 0.)
	 	# (6, 'stokes', 3000., 1024, 0., 1500., 0.)]


		print(inputfits[2].header)
		print(inputfits[2].columns.names)
		print(np.shape(inputfits[2].data))
		print(inputfits[2].data)

		print(inputfits[3].header)
		print(inputfits[3].columns.names)
		print(np.shape(inputfits[3].data))
		print(inputfits[3].data)

		print(inputfits[4].header)
		print(inputfits[4].columns.names)
		print(np.shape(inputfits[4].data))
		# print(inputfits[4].data)

		print(inputfits[5].header)
		print(inputfits[5].columns.names)
		print(np.shape(inputfits[5].data))
		print(inputfits[5].data)

		print(inputfits[6].header)
		print(inputfits[6].columns.names)
		print(np.shape(inputfits[6].data))
		print(inputfits[6].data)
		# print(inputfits[2].columns.names)
		# print(inputfits[3].columns.names)
		# # print(inputfits[3].data)
		# print(inputfits[4].columns.names)
		# print(inputfits[5].columns.names)

		# print(inputfits[2].data)
		# print(inputfits[3].data)
		# print(inputfits[4].data.field(3)[0])
		self.ra = np.asarray(np.degrees(inputfits[4].data['raj2000']))
		self.dec = np.asarray(np.degrees(inputfits[4].data['decj2000']))
		self.az = np.asarray(np.degrees(inputfits[4].data['az']))
		self.el = np.asarray(np.degrees(inputfits[4].data['el']))
		self.pa = np.asarray(np.degrees(inputfits[4].data['par_angle']))
		print(self.pa)
		exit()
		self.rotang = np.asarray(np.degrees(inputfits[4].data['derot_angle']))
		# self.ra = np.asarray(np.degrees(inputfits[4].data.field(3)))
		# self.dec = np.asarray(np.degrees(inputfits[4].data.field(4)))
		# self.az = np.asarray(np.degrees(inputfits[4].data.field(1)))
		# self.el = np.asarray(np.degrees(inputfits[4].data.field(2)))
		# print(ra)
		# print(dec)
		# print(az)
		# print(el)
		# exit()
		# ra = np.degrees(inputfits[4].data.field(1))
		# dec = np.degrees(inputfits[4].data.field(2))
		# az = np.degrees(inputfits[4].data.field(3))
		# el = np.degrees(inputfits[4].data.field(4))
		self.data = np.asarray(inputfits[4].data.field(-1), dtype=np.float64)
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
