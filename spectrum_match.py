# -*- coding: utf-8 -*-
# @Author: Charles Starr
# @Date:   2016-09-09 15:14:56
# @Last Modified by:   Charles Starr
# @Last Modified time: 2017-11-28 15:46:05

# This module intergrates a Library object and MassExperiment object
# and attempts to align theoretical ions with real spectra

import simple_library_constructor
import ms_ms_spectra
from pandas import DataFrame
from pandas import ExcelWriter
import numpy
from itertools import chain

class MassMatcher(object):
	# Class will contain the library and experimental spectrum data and
	# pick library members to attempt to match theoretical ions with 
	# experimental ones based on precursor exact mass


	def __init__(self, experiment, library, tolerance, outfile, rt_low, rt_high):

		self.mass_experiment = experiment
		self.peptide_library = library
		self.tolerance = tolerance
		self.rt_low = rt_low
		self.rt_high = rt_high
		self.outfile = outfile
		self.match_attempts = self.scan_ms_ms()
		self.tally_match_attempts()
		self.matches_tofile()

	def scan_ms_ms(self):
		# Iterates through list of ms/ms spectra selecting one at a 
		# time for analysis

		potential_matches = []

		for spectrum in self.mass_experiment.ms_ms_spectra:
			if spectrum.rt >= self.rt_low and spectrum.rt <= self.rt_high:
				new_matches = self.search_library(spectrum.precursor_exact_mass)
				potential_matches.extend([MatchAttempt(
					spectrum, peptide, spectrum.activation_method, self.tolerance)
					for peptide in new_matches]
					)
		
		return potential_matches

	def search_library(self, spectrum_exact_mass):
		# Binary search for library peptides that match precursor

		start = 0
		end = len(self.peptide_library.peptide_list) - 1
		matches = [] # Holds peptides that match

		while start <= end:
			index = (start+end)/2
			library_peptide = self.peptide_library.peptide_list[index]
			
			if abs(library_peptide.exact_mass 
				- spectrum_exact_mass) < self.tolerance:
				matches.append(library_peptide)
				neighbor = self.walk_library(index, spectrum_exact_mass)
				matches.extend(neighbor) 
				break
			
			elif library_peptide.exact_mass - spectrum_exact_mass < 0:
				start = index + 1
			
			else:
				end = index - 1

		return matches

	def walk_library(self, index, spectrum_exact_mass):
		# Finds peptides adjacent to the initially matched library member

		walked_matches = []
		up_index = index + 1
		dn_index = index - 1

		try:
			while abs(self.peptide_library.peptide_list[up_index].exact_mass 
				- spectrum_exact_mass) < self.tolerance:
				walked_matches.append(self.peptide_library.peptide_list[up_index])
				up_index += 1
		except:
			pass
		
		try:
			while abs(self.peptide_library.peptide_list[dn_index].exact_mass 
				- spectrum_exact_mass) < self.tolerance:
				walked_matches.append(self.peptide_library.peptide_list[dn_index])
				dn_index -= 1
		except:
			pass
		
		return walked_matches

	def tally_match_attempts(self):
		# Does simple analysis on peptides that were matched

		for attempt in self.match_attempts:
			match_dic = attempt.library_peptide.match_dict
			match_dic[attempt.activation_method]['match attempts'] += 1
			match_dic[attempt.activation_method]['matches'] += attempt.matches
			match_dic[attempt.activation_method]['match score'] += attempt.match_score
			attempt.library_peptide.scan_ids.append((
				attempt.spectrum.scan_number, attempt.matches,
				attempt.activation_method,
				map(lambda x: "%.2f" % x, [ions[0] for ions in attempt.matched_ions])
				)
			)
			attempt.library_peptide.precursor_intensities.append(
				attempt.spectrum.precursor_intensity
				)

		for peptide in self.peptide_library.peptide_list:
			peptide.scan_ids.sort(key=lambda x: x[1], reverse=True)
			peptide.calc_per_attempt()
			peptide.calc_total_score()

		return

	def matches_tofile(self):
		# Uses pandas dataframe to create a table and output an 
		# excel file for data summary.
		
		sorted_outlist = sorted(
			self.peptide_library.peptide_list, 
			key=lambda x: x.total_score, reverse=True
			)
		# Breaking line length timit to format the output for the program
		outframe = DataFrame({'Sequence':[peptide.sequence for peptide in sorted_outlist],
						'Exact mass':[peptide.exact_mass for peptide in sorted_outlist],
						'Precursor median':[numpy.median(peptide.precursor_intensities) if peptide.precursor_intensities else 0 for peptide in sorted_outlist],
						'Precursor mean':[numpy.mean(peptide.precursor_intensities) if peptide.precursor_intensities else 0 for peptide in sorted_outlist],
						'a ions':[map(lambda x: "%.2f" % x, peptide.ion_series['a_ions']) for peptide in sorted_outlist],
						'b ions':[map(lambda x: "%.2f" % x,peptide.ion_series['b_ions']) for peptide in sorted_outlist],
						'c ions':[map(lambda x: "%.2f" % x,peptide.ion_series['c_ions']) for peptide in sorted_outlist],
						'x ions':[map(lambda x: "%.2f" % x,peptide.ion_series['x_ions']) for peptide in sorted_outlist],
						'y ions':[map(lambda x: "%.2f" % x,peptide.ion_series['y_ions']) for peptide in sorted_outlist],
						'z ions':[map(lambda x: "%.2f" % x,peptide.ion_series['z_ions']) for peptide in sorted_outlist],
						'CID Match attempts':[peptide.match_dict['cid']['match attempts'] for peptide in sorted_outlist],
						'CID Matches':[peptide.match_dict['cid']['matches'] for peptide in sorted_outlist],
						'CID Matches/attempt':[peptide.match_dict['cid']['matches/attempt'] for peptide in sorted_outlist],
						'CID Match Score':[peptide.match_dict['cid']['match score'] for peptide in sorted_outlist],
						'HCD Match attempts':[peptide.match_dict['hcd']['match attempts'] for peptide in sorted_outlist],
						'HCD Matches':[peptide.match_dict['hcd']['matches'] for peptide in sorted_outlist],
						'HCD Matches/attempt':[peptide.match_dict['hcd']['matches/attempt'] for peptide in sorted_outlist],
						'HCD Match Score':[peptide.match_dict['hcd']['match score'] for peptide in sorted_outlist],
						'ETD Match attempts':[peptide.match_dict['etd']['match attempts'] for peptide in sorted_outlist],
						'ETD Matches':[peptide.match_dict['etd']['matches'] for peptide in sorted_outlist],
						'ETD Matches/attempt':[peptide.match_dict['etd']['matches/attempt'] for peptide in sorted_outlist],
						'ETD Match Score':[peptide.match_dict['etd']['match score'] for peptide in sorted_outlist],
						'Total Match Score':[peptide.total_score for peptide in sorted_outlist],
						'Scan IDs':[peptide.scan_ids for peptide in sorted_outlist]})
		
		# Reorders columns in dataframe
		outframe = outframe[
			['Sequence', 'Precursor median', 'Precursor mean', 'Exact mass',
			 'CID Match attempts', 'CID Matches', 'CID Matches/attempt', 'CID Match Score',
			 'HCD Match attempts', 'HCD Matches', 'HCD Matches/attempt', 'HCD Match Score',
			 'ETD Match attempts', 'ETD Matches', 'ETD Matches/attempt', 'ETD Match Score',
			 'Total Match Score', 'Scan IDs', 'a ions', 'b ions', 'c ions', 'x ions', 'y ions', 
			 'z ions'
			 ]
		]

		outframe.set_index('Sequence', inplace=True)
		writer = ExcelWriter(self.outfile)
		outframe.to_excel(writer)
		writer.save()

		return


class MatchAttempt(object):
	# Takes a single library peptide and experimental spectrum and
	# identifies ions that are present in both the theoretical and 
	# real spectra


	def __init__(self, spectrum, library_peptide, activation, tolerance):

		self.spectrum = spectrum
		self.library_peptide = library_peptide
		self.activation_method = activation
		self.tolerance = tolerance
		self.matched_ions = []
		self.ion_search()
		self.matches = len(self.matched_ions)
		self.match_score = self.matches * spectrum.precursor_intensity


	def ion_search(self):
		# Directs search to correct approach based on activation method

		if self.activation_method == 'etd':
			return self.etd_ion_search()

		elif self.activation_method == 'hcd':
			return self.cid_ion_search()

		elif self.activation_method == 'cid':
			return self.cid_ion_search()

		else:
			return

	def etd_ion_search(self):
		# Search for c/z ions in spectra created by etd fragmentation
		# Generates fragments if not already present

		if not self.library_peptide.ion_series['c_ions']:
			self.library_peptide.generate_c_ions()
		if not self.library_peptide.ion_series['z_ions']:
			self.library_peptide.generate_z_ions()

		self.ion_binary_search('c_ions')
		self.ion_binary_search('z_ions')
		
		return

	def cid_ion_search(self):
		# Search for b/y ions in spectra created by CID/HCD fragmentation
		# Generates fragments if not already present

		if not self.library_peptide.ion_series['b_ions']:
			self.library_peptide.generate_b_ions()
		if not self.library_peptide.ion_series['y_ions']:
			self.library_peptide.generate_y_ions()

		self.ion_binary_search('b_ions')
		self.ion_binary_search('y_ions')
		
		return

	def ion_binary_search(self, ion_type):
		# Searches an experimental spectrum ions of a specified type

		for ion in self.library_peptide.ion_series[ion_type]:
			start = 0
			end = len(self.spectrum.peaks) - 1
			
			while start <= end:
				index = (start + end)/2
				if abs(ion - self.spectrum.peaks[index][0]) < self.tolerance:
					if self.spectrum.peaks[index] not in self.matched_ions:
						self.matched_ions.append(self.spectrum.peaks[index])
					else:
						self.matched_ions.extend(self.walk_spectrum(ion, index))
					break
				elif ion - self.spectrum.peaks[index][0] < 0:
					end = index - 1
				else:
					start = index + 1
		
		return

	def walk_spectrum(self, ion, index):
		# In the case that a matched ion has already been found in a 
		# different ion series, check to see if there are any other 
		# ions that would match the criterion

		try:
			up_ion = abs(ion - self.spectrum.peaks[index + 1][0])
			if up_ion <= self.tolerance:
				return [self.spectrum.peaks[index + 1]]
		except:
			pass

		try:
			down_ion = abs(ion - self.spectrum.peaks[index - 1][0])
			if down_ion <= self.tolerance:
				return [self.spectrum.peaks[index -1]]
		except:
			pass
		
		return []












