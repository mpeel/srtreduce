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
	print(navg)
	print(npoints//navg)
	for i in range(0,npoints//navg):
		start = i*navg
		end = ((i+1)*navg)
		print('start:' + str(start) + ', end: ' + str(end) + ', navg: ' + str(navg))
		if option == 0:
			data[start:end] = data[start:end] - np.median(data[start:end])
		else:
			xline = range(0,len(data[start:end]))
			if len(xline) != 0:
				A,B=curve_fit(linfit,xline,data[start:end])[0]
				# print(A,B)
				data[start:end] = data[start:end] - (A*xline+B)
	if navg % npoints:
		if option == 0:
			data[(npoints//navg)*navg-1:] = data[(npoints//navg)*navg-1:] - np.median(data[(npoints//navg)*navg-1:])
		else:
			xline = range(0,len(data[(npoints//navg)*navg-1:]))
			if len(xline) != 0:
				A,B=curve_fit(linfit,xline,data[(npoints//navg)*navg-1:])[0]
				# print(A,B)
				data[(npoints//navg)*navg-1:] = data[(npoints//navg)*navg-1:] - (A*xline+B)

	return data
