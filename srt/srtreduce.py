#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Galactic, ICRS, FK5
from srtobservation import *
from srtdata import *


class srtreduce():
	def __init__(self,basedir="", outdir='', project='', band='K', num_beams = 0):
		# Directories
		self.basedir = basedir
		self.outdir = outdir
		os.makedirs(self.outdir, exist_ok=True)

		# Basic telescope information
		self.is_set = False
		self.telescope = EarthLocation(lat=0.689283579621821*u.rad, lon=0.161358481873679*u.rad, height=650*u.m)
		self.num_beams = num_beams
		self.pix_xoff = np.zeros(0)
		self.pix_yoff = np.zeros(0)
		self.pix_rel = np.zeros(0)
		self.nbins = 0
		self.bw = 0
		self.lo = 0
		self.rotang = 0.0

		# Analysis options
		self.project = project
		self.band = band
		# Constants
		self.k = 1.380e-23 # Boltzman constant m2 kg s-2 K-1
		self.c = 2.9979e8

		# Stable version control
		self.ver = 0.0

	def set_params(self, inputfits):
		self.telescope = EarthLocation(lat=inputfits[0].header['SiteLatitude']*u.rad, lon=inputfits[0].header['SiteLongitude']*u.rad, height=inputfits[0].header['SiteHeight']*u.m)
		if self.project == '':
			self.project = inputfits[0].header['Project_Name'].strip()
			# print(self.project)
		if self.num_beams == 0:
			self.num_beams = inputfits[0].header['BEAMS']
			# print(self.num_beams)
		for i in range(0,self.num_beams):
			# print(inputfits[3].data[i])
			self.pix_xoff = np.append(self.pix_xoff, inputfits[3].data[i][1])
			self.pix_yoff = np.append(self.pix_yoff, inputfits[3].data[i][2])
			self.pix_rel = np.append(self.pix_rel, inputfits[3].data[i][3])
		# print(inputfits[3].data)
		# print(self.pix_xoff)
		# print(self.pix_yoff)
		# print(self.pix_rel)
		self.nbins = inputfits[1].data[0]['bins']
		self.bw = inputfits[1].data[0]['bandWidth']
		self.lo = inputfits[2].data[0]['localOscillator']
		self.calc_freqs()
		self.rotang = inputfits[3].header['DEWUSER']
		if self.rotang < -360.0:
			self.rotang = 0.0
		return

	def calc_freqs(self):
		self.freqs = self.lo + np.arange(0,self.bw,self.bw/self.nbins)
		# print(len(self.freqs))
		# print(self.freqs)
		return

	def get_pixpos(self, pixnum):
		return self.pix_xoff[pixnum], self.pix_yoff[pixnum], self.pix_rel[pixnum]
