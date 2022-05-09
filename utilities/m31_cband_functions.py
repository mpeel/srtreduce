#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Functions for analysing the C-Band M31 map from SRT
# Main code: m31_cband.py
# 
# Version history:
#
# 17-Apr-2020  M. Peel       Started

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
from photutils import EllipticalAperture, EllipticalAnnulus, aperture_photometry

def getmap(infile):
	hdul = fits.open(infile)
	header = hdul[0].header
	if np.shape(hdul[0].data)[0] == 1:
		data = np.asarray(hdul[0].data[0][0])
	else:
		data = np.asarray(hdul[0].data)
	return header, data

def plotmap(header,data,outfile,pltrange=(0,0)):
	# wcs = WCS(header)
	# slices = ('x','y',0,0)

	# ax = plt.subplot(projection=wcs,slices=slices)
	# plt.axis('off')
	plt.subplot()
	plt.tight_layout()
	if pltrange != (0,0):
		plt.imshow(data,clim=pltrange)
	else:
		plt.imshow(data)
	# plt.title(galaxy + ' convolved')
	cbar = plt.colorbar()
	# cbar.set_label('mJy')
	plt.savefig(outfile)
	plt.clf()
	return

def smoothmap(data,fwhm,pixelsize):

	gauss = Gaussian2DKernel(np.abs(fwhm/pixelsize))
	result = convolve_fft(data, gauss)
	return result

def aperflux(inmap, freq, res_arcmin, lon, lat, aper_inner_radius, aper_outer_radius1, aper_outer_radius2, units, column=0, dopol=False, nested=False, noise_model=0, abratio=1.0, angle=0.0, silent=False,quickplot=''):
	positions = (lon,lat)
	aperture = EllipticalAperture(positions, a=int(np.abs(aper_inner_radius)),b=int(np.abs(aper_inner_radius*abratio)),theta=angle*np.pi/180.0)
	annulus_aperture = EllipticalAnnulus(positions, a_in=int(np.abs(aper_outer_radius1)), a_out=int(np.abs(aper_outer_radius2)), b_out=int(np.abs(aper_outer_radius2*abratio)),theta=angle*np.pi/180.0)
	aperture_mask = aperture.to_mask(method='center').multiply(inmap)
	annulus_aperture_mask = annulus_aperture.to_mask(method='center')
	annulus_aperture_mask_data = annulus_aperture_mask.multiply(inmap)
	plt.imshow(aperture_mask)
	plt.colorbar()
	plt.savefig(quickplot+'_aperture.png')
	plt.clf()
	plt.imshow(annulus_aperture_mask_data)
	plt.colorbar()
	plt.savefig(quickplot+'_annulus.png')
	plt.clf()
	inner = aperture_photometry(inmap, aperture)[0]['aperture_sum']
	print(np.shape(annulus_aperture_mask))
	outer = np.median(annulus_aperture_mask_data[annulus_aperture_mask.data>0])*aperture.area
	# outer = aperture_photometry(inmap, annulus_aperture)['aperture_sum']
	return inner,outer, inner-outer
