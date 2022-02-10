#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import os


ext = 'radec.txt'
numext = 1
inputlist = os.listdir('/Users/mpeel/Documents/git/astrocode/srt')
filelist = [f for f in inputlist if ext in f]
print(filelist)
filelist.sort()
print(filelist)
i = 0
for file in filelist:
	ra, dec = np.loadtxt(file,unpack=True)
	plt.plot(ra, dec,',',label=file.replace('_radec.txt',''))
	fig = plt.gcf()
	fig.set_size_inches(10, 8)
	plt.xlim([50,62])
	plt.ylim([27,36])
	plt.title('Summary RA/Dec coverage')
	plt.xlabel('RA [deg]')
	plt.ylabel('Dec [deg]')
	plt.legend(loc='upper right',fontsize='x-small')
	plt.tight_layout()
	plt.savefig('summary_ra_dec_'+'{0:02d}'.format(i)+'_ind.png')
	i+=1
	plt.clf()
