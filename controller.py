# -*- coding: utf-8 -*-
# @Author: Charles Starr
# @Date:   2016-09-09 17:49:30
# @Last Modified by:   Charles Starr
# @Last Modified time: 2017-11-28 15:28:36

# This module integrates all other elements of chuck_seq

import simple_library_constructor
import spectrum_match
import ms_ms_spectra


def main(template_csv, ms_datafile, filetype, tolerance, outfile, rt_low = 0.0, rt_high = 200.0):

	mass_experiment = ms_ms_spectra.MassExperiment(ms_datafile, filetype)
	library = simple_library_constructor.Library(template_csv)
	mass_matcher = spectrum_match.MassMatcher(mass_experiment, library, tolerance, outfile, rt_low, rt_high)

main('phd_csv.csv', '1458_Sample_7.mzML', 'mzml', 0.1, 'PHD_scoretest.xlsx')

