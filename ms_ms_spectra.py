# -*- coding: utf-8 -*-
# @Author: Charles Starr
# @Date:   2016-09-09 13:37:14
# @Last Modified by:   Charles Starr
# @Last Modified time: 2016-09-09 15:09:08

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

	def parse_mgf(self):
		ms_ms_list = []
		with open(self.ms_datafile) as ms_ms_data:
			for line in ms_ms_data:
				if line == 'BEGIN IONS\n':
					scan_number = int(ms_ms_data.next().strip('\n').split('.')[1]) #split string to get scan number
					retention_time = float(ms_ms_data.next().strip('\n').split('=')[1]) #split string to get retention time
					pep_mz, intensity = map(float, ms_ms_data.next().strip('\n').split('=')[1].split()) #split to get m/z and intensity
					charge = int(ms_ms_data.next().strip('\n')[7])
					peaks = []
					while True:
						peak_line = ms_ms_data.next().strip('\n')
						if peak_line != 'END IONS':
							peaks.append(tuple(map(float, peak_line.split(' '))))
						else:
							ms_ms_list.append(MS_MS_Spectra(scan_number, retention_time, pep_mz, intensity,
															charge, peaks))
							break
		return ms_ms_list


class MS_MS_Spectra(object):
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