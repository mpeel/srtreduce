#!/usr/bin/env python
# -*- coding: utf-8  -*-
#
# Analyse the C-Band M31 map from SRT
# 
# Version history:
#
# 17-Apr-2020  M. Peel       Started

import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from m31_cband_functions import *

# Configuration
basedir = "/Users/mpeel/Documents/SRT_M31/"
indir = basedir+'LATEST_CBAND/'
outdir = basedir+'analyse/'

rawmaps = ['SRT_6156.FITS','SRT_6469.FITS','SRT_6781.FITS','SRT_7094.FITS','SPLAT_6000_7250_AVERAGE_DIRECT_STACK_I_SWT_BLANK.FITS']
corrmaps = ['SRT_6156_SUB_BASELEV_CORR.FITS','SRT_6469_SUB_BASELEV_CORR.FITS','SRT_6781_SUB_BASELEV_CORR.FITS','SRT_7094_SUB_BASELEV_CORR.FITS','SRT_6000_7250_SUB_BASELEV_CORR.FITS']

beamsize = 2.7 # arcmin
pixsize = 0.9 # arcmin
beamarea = np.pi*beamsize**2/(4.0*np.log(2.0))
beamtopix = beamarea*(pixsize**2)
# print(beamarea)
# print(beamtopix)
# exit()
output = 33.0 # arcmin
aper_r1 = 60.0
aper_r2 = 60.0
aper_r3 = 80.0
# Options for what to run
plotmaps = True

if plotmaps:
	for i in range(0,len(rawmaps)):
		filename = rawmaps[i]
		header,data=getmap(indir+filename)
		pixelsize = header['CDELT1']*60.0
		rescale = np.abs(beamarea*(pixelsize**2))
		data = data/beamtopix
		plotmap(header,data,outdir+filename.replace('.FITS','.png'),pltrange=(0,0.005))
		# print(pixelsize)
		smooth = smoothmap(data,np.sqrt(output**2-beamsize**2),pixelsize)
		plotmap(header,smooth,outdir+filename.replace('.FITS','_smooth.png'))

		lon = np.shape(data)[0]/2
		lat = np.shape(data)[1]/2
		aper_inner_radius = int(aper_r1 / pixelsize)
		aper_outer_radius1 = int(aper_r2 / pixelsize)
		aper_outer_radius2 = int(aper_r3 / pixelsize)
		# print(aper_inner_radius)
		print(aperflux(data, 5.0, beamsize, lon, lat, aper_inner_radius, aper_outer_radius1, aper_outer_radius2, 'JyPix', abratio=0.7, angle=-45.0,quickplot=outdir+filename.replace('.FITS','aperture')))
		# print(aperflux(smooth, 5.0, beamsize, lon, lat, aper_inner_radius, aper_outer_radius1, aper_outer_radius2, 'JyBeam', abratio=0.7, angle=-45.0,quickplot=outdir+filename.replace('.FITS','smoothaperture')))
		# exit()
		data_raw = data.copy()

		filename = corrmaps[i]
		header,data=getmap(indir+filename)
		pixelsize = header['CDELT1']*60.0
		rescale = np.abs(beamarea/(pixelsize**2))
		data = data/beamtopix
		plotmap(header,data,outdir+filename.replace('.FITS','.png'),pltrange=(0,0.005))
		smooth = smoothmap(data,np.sqrt(output**2-beamsize**2),pixelsize)
		plotmap(header,smooth,outdir+filename.replace('.FITS','_smooth.png'))
		
		lon = np.shape(data)[0]/2
		lat = np.shape(data)[1]/2
		aper_inner_radius = int(aper_r1 / pixelsize)
		aper_outer_radius1 = int(aper_r2 / pixelsize)
		aper_outer_radius2 = int(aper_r3 / pixelsize)
		print(aper_inner_radius)
		print(aperflux(data, 5.0, beamsize, lon, lat, aper_inner_radius, aper_outer_radius1, aper_outer_radius2, 'JyPix', abratio=0.7, angle=-45.0,quickplot=outdir+filename.replace('.FITS','aperture')))
		# print(aperflux(smooth, 5.0, beamsize, lon, lat, aper_inner_radius, aper_outer_radius1, aper_outer_radius2, 'JyBeam', abratio=0.7, angle=-45.0,quickplot=outdir+filename.replace('.FITS','smoothaperture')))

		diff = data_raw-data
		diff[~np.isfinite(diff)] = 0.0
		plotmap(header,diff,outdir+filename.replace('.FITS','_diff.png'),pltrange=(-0.005,0.005))
		print(np.sum(diff[diff > 0.001]))
		# exit()