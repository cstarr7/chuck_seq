# -*- coding: utf-8 -*-
# @Author: Charles Starr
# @Date:   2016-09-09 13:37:14
# @Last Modified by:   Charles Starr
# @Last Modified time: 2016-09-14 01:00:34
import pymzml

class Mass_Experiment(object):
	#this class will be a container for all MS/MS spectra acquired in a single LC injection
	def __init__(self, ms_datafile, filetype):
		self.ms_datafile = ms_datafile
		self.filetype = filetype
		self.ms_ms_spectra = self.parse_datafile()

	def parse_datafile(self):
		#selects a parsing fuction based on filetype
		if self.filetype == 'mgf':
			ms_ms_spectra = self.parse_mgf()
			return ms_ms_spectra
		if self.filetype == 'mzml':
			ms_ms_spectra = self.parse_mzml()
			return ms_ms_spectra

	def parse_mgf(self):
		ms_ms_list = []
		with open(self.ms_datafile) as ms_ms_data:
			for line in ms_ms_data:
				if line == 'BEGIN IONS\n':
					scan_number = int(ms_ms_data.next().strip('\n').split('.')[1]) #split string to get scan number
					retention_time = float(ms_ms_data.next().strip('\n').split('=')[1]) #split string to get retention time
					print scan_number
					mz_intensity = ms_ms_data.next().strip('\n').split('=')[1]
					try:
						pep_mz, intensity = map(float, mz_intensity.split()) #split to get m/z and intensity
					except:
						pep_mz = float(mz_intensity)#sometimes no intensity so we have to deal separately
						intensity = None
					charge = int(ms_ms_data.next().strip('\n')[7]) #get charge, questionable approach
					peaks = [] #empty list to hold upcoming peak information
					while True:
						peak_line = ms_ms_data.next().strip('\n')
						if peak_line != 'END IONS':
							peaks.append(tuple(map(float, peak_line.split(' ')))) #add peak m/z and intensity as tuple
						else:
							ms_ms_list.append(MS_MS_Spectrum(scan_number, retention_time, pep_mz, intensity,
															charge, peaks)) #create spectrum object
							break
		return ms_ms_list

	def parse_mzml(self):
		msrun = pymzml.run.Reader(self.ms_datafile)
		ms_ms_list = []
		for i, spectrum in enumerate(msrun,1):
			if spectrum['ms level'] == 2:
				scan_number = i
				retention_time = spectrum['MS:1000016']
				pep_mz = spectrum['MS:1000744']
				try:
					intensity = spectrum['MS:1000042']
				except:
					intensity = 0
				charge = spectrum['MS:1000041']
				#fragment_type = spectrum['activation']
				peaks = spectrum.peaks
				ms_ms_list.append(MS_MS_Spectrum(scan_number, retention_time, pep_mz, intensity,
															charge, peaks))
		return ms_ms_list

class MS_MS_Spectrum(object):
	#this class will hold peak information for an individual MS/MS shot
	def __init__(self, scan, rt, precursor_mz, precursor_intensity, precursor_charge, peaks):
		self.scan_number = scan
		self.retention_time = rt
		self.precursor_mz = precursor_mz
		self.precursor_intensity = precursor_intensity
		self.precursor_charge = precursor_charge
		self.peaks = peaks
		self.precursor_exact_mass = self.calc_precursor()

	def calc_precursor(self):
		exact_mass = self.precursor_mz * self.precursor_charge - self.precursor_charge
		return exact_mass

'''
test_experiment = Mass_Experiment('etd_reanalyzed.mgf', 'mgf')
for spectrum in test_experiment.ms_ms_spectra:
	print spectrum.scan_number
	print spectrum.retention_time
	print spectrum.precursor_mz
	print spectrum.precursor_intensity
	print spectrum.precursor_charge
	print spectrum.peaks
	print spectrum.precursor_exact_mass
'''