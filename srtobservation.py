from srtdata import *
import os

class srtobservation():
	def __init__(self,foldername=""):
		# Directories
		self.foldername = foldername
		self.filelist = []
		self.get_filelist()

	def get_filelist(self):
		inputlist = os.listdir(self.foldername)
		for filename in inputlist:
			if 'summary' not in filename:
				self.filelist.append(filename)
		self.filelist.sort()
		return
