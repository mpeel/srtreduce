#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
import astropy.units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, Galactic, ICRS, FK5
from srt_observation import *
from srt_data import *

class srt_reduce():
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
		self.pix_rel = np.zeros(0) # Relative power
		self.attenuations = np.zeros(0)
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

	# Read in basic parameters from an already open data fits file
	# NOTE: this assumes they are the same for all observations that you'll work on with this instance.
	def set_params(self, inputfits):
		# Overall telescope and receiver properties
		self.telescope = EarthLocation(lat=inputfits[0].header['SiteLatitude']*u.rad, lon=inputfits[0].header['SiteLongitude']*u.rad, height=inputfits[0].header['SiteHeight']*u.m)
		if self.project == '':
			self.project = inputfits[0].header['Project_Name'].strip()
		if self.num_beams == 0:
			self.num_beams = inputfits[0].header['BEAMS']
		# Frequency parameters
		self.nbins = inputfits[1].data[0]['bins']
		# self.bw = inputfits[1].data[0]['bandWidth'] # This reads in 1500?
		self.bw = inputfits[2].data[0]['bandWidth'] # This reads in 1400?
		self.lo = inputfits[2].data[0]['localOscillator']
		# Pixel positions and rotation
		for i in range(0,self.num_beams):
			self.attenuations = np.append(self.attenuations, inputfits[2].data[i*2]['attenuation'])
			self.pix_xoff = np.append(self.pix_xoff, inputfits[3].data[i][1])
			self.pix_yoff = np.append(self.pix_yoff, inputfits[3].data[i][2])
			self.pix_rel = np.append(self.pix_rel, inputfits[3].data[i][3])
		self.rotang = inputfits[3].header['DEWUSER']
		if self.rotang < -360.0:
			self.rotang = 0.0
		return

	def calc_freqs(self):
		# NOTE: Is this right, or need to offset by 0.5*bw/nbins?
		self.freqs = self.lo + np.arange(0,self.bw,self.bw/self.nbins)
		return

	def get_pixpos(self, pixnum):
		return self.pix_xoff[pixnum], self.pix_yoff[pixnum], self.pix_rel[pixnum]
