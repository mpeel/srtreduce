#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os
from scipy.optimize import curve_fit
# from srt_functions import *
from photutils import CircularAperture, aperture_photometry, CircularAnnulus

def calc_beam_area(beam1,beam2):
	return (np.pi * (beam1*np.pi/180.0)*(beam2*np.pi/180.0))/(4.0*np.log(2.0))
def do_aperture(filename):
	positions = [(64., 64.)]
	aperture = CircularAperture(positions, r=25.)
	annulus_aperture = CircularAnnulus(positions, r_in=25, r_out=35)

	m51high = fits.open(filename)
	print(m51high.info())
	print(m51high[0].header['CDELT1']*60)
	dataset = m51high[0].data[0]
	print(dataset)
	print(np.shape(dataset))
	dataset[~np.isfinite(dataset)]=0.0
	m51sum = np.nansum(dataset[dataset>-0.0001])
	print(m51sum)
	print(np.nanmax(dataset))
	print(np.nanmin(dataset))

	phot_table = aperture_photometry(dataset, aperture)
	print(phot_table)
	bg_table = aperture_photometry(dataset, annulus_aperture)
	print(bg_table)

	beamsize = 0.8
	beamarea = calc_beam_area(beamsize/60,beamsize/60)
	print(beamarea)
	pixelsize = 0.25
	pixelarea = ((pixelsize/60)*(np.pi/180.0))**2
	print(pixelarea)
	print(pixelarea/beamarea)
	print(m51sum*pixelarea/beamarea)
	value = m51sum*pixelarea/beamarea
	return value

n6946_low = do_aperture('/Users/mpeel/Documents/SRT/nearbygal/N6946_18.6GHZ/SPLAT_N6946_SWT_I_AVERAGE_DENOISE.FITS')

# m51_high = do_aperture('/Users/mpeel/Documents/SRT/nearbygal/SPLAT_M51high_SWT_I_AVERAGE_DENOISE.FITS')
# m51_low = do_aperture('/Users/mpeel/Documents/SRT/nearbygal/SPLAT_M51low_SWT_I_AVERAGE_DENOISE.FITS')
# print(np.log10(m51_high/m51_low)/np.log10(24.6/18.6))
