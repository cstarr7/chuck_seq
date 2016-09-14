# -*- coding: utf-8 -*-
# @Author: Charles Starr
# @Date:   2016-09-09 15:14:56
# @Last Modified by:   Charles Starr
# @Last Modified time: 2016-09-14 01:04:00
import simple_library_constructor
import ms_ms_spectra
from pandas import DataFrame
from pandas import ExcelWriter

class Mass_Matcher(object):
	#class will contain the library and experimental spectrum data and pick library members to 
	#attempt to match theoretical ions with experimental ones based on precursor exact mass
	def __init__(self, experiment, library, frag_type, tolerance, outfile):
		self.mass_experiment = experiment
		self.peptide_library = library
		self.frag_type = frag_type
		self.tolerance = tolerance
		self.outfile = outfile
		self.match_attempts = self.scan_ms_ms()
		self.tally_match_attempts()
		self.matches_tofile()

	def scan_ms_ms(self):
		#iterates through list of ms/ms spectra selecting one at a time for analysis
		potential_matches = []
		for spectrum in self.mass_experiment.ms_ms_spectra:
			new_matches = self.search_library(spectrum.precursor_exact_mass)
			potential_matches.append([Match_Attempt(spectrum, peptide, self.frag_type, self.tolerance) for peptide in new_matches])
		return [match for matches in potential_matches for match in matches]

	def search_library(self, spectrum_exact_mass):
		start = 0
		end = len(self.peptide_library.peptide_list) - 1
		matches = []
		while start <= end:
			index = (start+end)/2
			library_peptide = self.peptide_library.peptide_list[index]
			if abs(library_peptide.exact_mass - spectrum_exact_mass) < self.tolerance: #library peptide within tolerance?
				matches.append(library_peptide) #add id'd peptide to matches
				neighbor_matches = self.walk_library(index, spectrum_exact_mass) #walk up and down list until no matches
				matches.extend(neighbor_matches) #add walked matches to matches
				break
			elif library_peptide.exact_mass - spectrum_exact_mass < 0:
				start = index + 1
			else:
				end = index - 1
		return matches

	def walk_library(self, index, spectrum_exact_mass):
		#walk up and down the peptide list from first identified match to find all matches
		walked_matches = []
		up_index = index + 1
		down_index = index - 1
		while abs(self.peptide_library.peptide_list[up_index].exact_mass - spectrum_exact_mass) < self.tolerance:
			walked_matches.append(self.peptide_library.peptide_list[up_index])
			up_index += 1
		while abs(self.peptide_library.peptide_list[down_index].exact_mass - spectrum_exact_mass) < self.tolerance:
			walked_matches.append(self.peptide_library.peptide_list[down_index])
			down_index -= 1
		return walked_matches

	def tally_match_attempts(self):
		#does simple analysis across runs to get some descriptive parameters
		for attempt in self.match_attempts:
			attempt.library_peptide.match_attempts += 1
			attempt.library_peptide.matches += attempt.matches
			attempt.library_peptide.scan_ids.append((attempt.spectrum.scan_number, attempt.matches))
		for peptide in self.peptide_library.peptide_list:
			peptide.scan_ids.sort(key=lambda x: x[1], reverse=True)
			peptide.calc_per_attempt()

	def matches_tofile(self):
		#uses pandas dataframe to create a table and output an excel file for data summary
		sorted_outlist = sorted(self.peptide_library.peptide_list, key=lambda x: x.match_attempts, reverse=True)
		outframe = DataFrame({'Sequence':[peptide.sequence for peptide in sorted_outlist],
						'Exact mass':[peptide.exact_mass for peptide in sorted_outlist],
						'a ions':[peptide.ion_series['a_ions'] for peptide in sorted_outlist],
						'b ions':[peptide.ion_series['b_ions'] for peptide in sorted_outlist],
						'c ions':[peptide.ion_series['c_ions'] for peptide in sorted_outlist],
						'x ions':[peptide.ion_series['x_ions'] for peptide in sorted_outlist],
						'y ions':[peptide.ion_series['y_ions'] for peptide in sorted_outlist],
						'z ions':[peptide.ion_series['z_ions'] for peptide in sorted_outlist],
						'Match attempts':[peptide.match_attempts for peptide in sorted_outlist],
						'Matches':[peptide.matches for peptide in sorted_outlist],
						'Matches/attempt':[peptide.matches_per_attempt for peptide in sorted_outlist],
						'Scan IDs':[peptide.scan_ids for peptide in sorted_outlist]})
		outframe.set_index('Sequence', inplace=True)
		writer = ExcelWriter(self.outfile)
		outframe.to_excel(writer)
		writer.save()

class Match_Attempt(object):
	#takes a single library peptide and experimental spectrum and ID's common ions
	def __init__(self, spectrum, library_peptide, frag_type, tolerance):
		self.spectrum = spectrum
		self.library_peptide = library_peptide
		self.frag_type = frag_type
		self.tolerance = tolerance
		self.matched_ions = []
		self.ion_search()
		self.matches = len(self.matched_ions)

	def ion_search(self):
		#directs search to correct approach based on fragmentation method
		if self.frag_type == 'ETD':
			return self.etd_ion_search()
		elif self.frag_type == 'HCD':
			return self.cid_ion_search()
		elif self.frag_type == 'CID':
			return self.cid_ion_search()
		else:
			return

	def etd_ion_search(self):
		#search for c and z ions in spectra created by etd fragmentation
		#generates fragments if not already present
		if not self.library_peptide.ion_series['c_ions']:
			self.library_peptide.generate_c_ions()
		if not self.library_peptide.ion_series['z_ions']:
			self.library_peptide.generate_z_ions()
		self.ion_binary_search('c_ions')
		self.ion_binary_search('z_ions')
		return

	def cid_ion_search(self):
		#search for b and b ions in spectra created by etd fragmentation
		#generates fragments if not already present
		if not self.library_peptide.ion_series['b_ions']:
			self.library_peptide.generate_c_ions()
		if not self.library_peptide.ion_series['y_ions']:
			self.library_peptide.generate_z_ions()
		self.ion_binary_search('b_ions')
		self.ion_binary_search('y_ions')
		return

	def ion_binary_search(self, ion_type):
		#searches an experimental spectrum ions of a specified type
		for ion in self.library_peptide.ion_series[ion_type]: #match c_ions with binary search
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
		#in the case that a matched ion has already been found in a different search
		#to see if there are any other ions that would match the criterion
		up_ion = abs(ion - self.spectrum.peaks[index + 1][0])
		down_ion = abs(ion - self.spectrum.peaks[index - 1][0])
		if up_ion <= self.tolerance:
			return [self.spectrum.peaks[index + 1]]
		elif down_ion <= self.tolerance:
			return [self.spectrum.peaks[index -1]]
		else:
			return []












