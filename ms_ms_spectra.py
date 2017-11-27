# -*- coding: utf-8 -*-
# @Author: Charles Starr
# @Date:   2016-09-09 13:37:14
# @Last Modified by:   Charles Starr
# @Last Modified time: 2017-11-27 15:03:15

# The purpose of this module is to handle MS/MS data from a variety of 
# input formats. It currently supports mgf and mzml filetypes. The 
# MassExperiment object handles the unpacking of the file and creates
# an MSMSSpectrum object for each MS/MS shot in the file. The pymzml
# package is used to handle mzml files.

import pymzml


class MassExperiment(object):
	# This class will be a container for all MS/MS spectra acquired in
	# a single LC injection


	def __init__(self, ms_datafile, filetype):

		self.ms_datafile = ms_datafile
		self.filetype = filetype
		self.ms_ms_spectra = self.parse_datafile()

	def parse_datafile(self):
		# Selects a parsing fuction based on filetype.

		if self.filetype == 'mgf':
			ms_ms_spectra = self.parse_mgf()

		if self.filetype == 'mzml':
			ms_ms_spectra = self.parse_mzml()
			
		return ms_ms_spectra

	def parse_mgf(self):
		# Instructions for unpacking an mgf file.

		ms_ms_list = []

		with open(self.ms_datafile) as ms_ms_data:
			for line in ms_ms_data:
				if line == 'BEGIN IONS\n':

					scan_number = int(
						ms_ms_data.next().strip('\n').split('.')[1]
						)

					retention_time = float(
						ms_ms_data.next().strip('\n').split('=')[1]
						)
					
					mz_intensity = ms_ms_data.next().strip('\n').split('=')[1]

					# Sometimes intensity is not included with mz
					try:
						pep_mz, intensity = map(float, mz_intensity.split())
					except:
						pep_mz = float(mz_intensity)
						intensity = 0.0

					charge = int(ms_ms_data.next().strip('\n')[7])

					peaks = []
					while True:
						peak_line = ms_ms_data.next().strip('\n')
						if peak_line != 'END IONS':
							peaks.append(
								tuple(map(float, peak_line.split(' '))))
						else:
							ms_ms_list.append(MSMSSpectrum(
								scan_number, retention_time, pep_mz, 
								intensity, charge, peaks)
								)
							break
		
		return ms_ms_list

	def parse_mzml(self):
		# Instructions for parsing an mzml file.

		msrun = pymzml.run.Reader(self.ms_datafile, extraAccessions=[('MS:1000042', ['value'])])
		ms_ms_list = []

		for i, spectrum in enumerate(msrun,1):
			if spectrum['ms level'] == 2:
				scan_number = i
				retention_time = spectrum['MS:1000016']
				pep_mz = spectrum['MS:1000744']
				# Sometimes there is no intensity value.
				try:
					intensity = spectrum['MS:1000042']
				except:
					intensity = 0.0

				charge = spectrum['MS:1000041']
				activation_method = spectrum['MS:1000512'].split('@')[1][:3]
				peaks = spectrum.peaks
				ms_ms_list.append(MSMSSpectrum(
					scan_number, retention_time, pep_mz, 
					intensity, charge, peaks, activation_method)
					)

		return ms_ms_list


class MSMSSpectrum(object):
	# Contains useful information about each individual MS/MS shot


	def __init__(
		self, scan, rt, precursor_mz, precursor_intensity, 
		precursor_charge, peaks, activation_method
		):
		
		self.scan_number = scan
		self.rt = rt
		self.precursor_mz = precursor_mz
		self.precursor_intensity = precursor_intensity
		self.precursor_charge = precursor_charge
		self.peaks = peaks
		self.activation_method = activation_method
		self.precursor_exact_mass = self.calc_precursor()

	def calc_precursor(self):
		# Calculate the exact mass of the MS/MS precursor ion

		exact_mass = (self.precursor_mz 
			* self.precursor_charge 
			- self.precursor_charge
			)
		
		return exact_mass