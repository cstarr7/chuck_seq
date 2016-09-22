# -*- coding: utf-8 -*-
# @Author: Charles Starr
# @Date:   2016-09-09 17:49:30
# @Last Modified by:   Charles Starr
# @Last Modified time: 2016-09-21 17:27:18

# This module integrates all other elements of chuck_seq

import simple_library_constructor
import spectrum_match
import ms_ms_spectra


def main(template_csv, ms_datafile, filetype, tolerance, outfile):

	mass_experiment = ms_ms_spectra.MassExperiment(ms_datafile, filetype)
	library = simple_library_constructor.Library(template_csv)
	mass_matcher = spectrum_match.MassMatcher(mass_experiment, library, tolerance, outfile)

main('arva_csv.csv', '1444_1_cid_de3mscan3.mzML', 'mzml', 0.1, '1444_1_cid_de3mscan3.xlsx')

