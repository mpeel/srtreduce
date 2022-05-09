#!/usr/bin/env python
# This class handles reading in an SRT observation and analysing it
# Mike Peel

from srtfunctions import *
from srtdata import *
import os

class srtobservation():
	def __init__(self,srtreduce='',foldername=""):
		# Directories
		self.foldername = foldername
		self.filelist = []

		# Telescope
		self.srtreduce = srtreduce

		# Data
		self.data = []

		# Auto-import the filelist and data
		self.get_filelist()
		self.read_data()

	def get_filelist(self):
		inputlist = os.listdir(self.foldername)
		for filename in inputlist:
			if 'summary' not in filename:
				self.filelist.append(filename)
		self.filelist.sort()
		return

	def read_data(self):
		for filename in self.filelist:
			self.data.append(srtdata(self.srtreduce, self.foldername+filename))
			print('Hi')
			print(self.data[0].filename)
			print(self.data[0].ra)
			exit()
