# -*- coding: utf-8 -*-
# @Author: Charles Starr
# @Date:   2016-09-09 17:49:30
# @Last Modified by:   Charles Starr
# @Last Modified time: 2016-09-12 23:08:41
import simple_library_constructor
import spectrum_match
import ms_ms_spectra


def main(template_csv, ms_datafile, filetype, fragment_type, tolerance, outfile):
	mass_experiment = ms_ms_spectra.Mass_Experiment(ms_datafile, filetype)
	library = simple_library_constructor.Library(template_csv)
	mass_matcher = spectrum_match.Mass_Matcher(mass_experiment, library, fragment_type, tolerance, outfile)
main('arva_csv.csv', 'etd_reanalyzed.mgf', 'mgf', 'ETD', 0.1, 'out.xlsx')

