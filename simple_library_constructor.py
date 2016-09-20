# -*- coding: utf-8 -*-
"""
Created on Mon Oct 05 15:50:24 2015

@author: Charlie
"""

# The purpose of this module is to create a combinatorial peptide 
# library from user defined library template. The template should 
# be a CSV file with positions defined by row and variable residues
# defined in the columns within each row. The program understands the
# 20 naturally occurring amino acids as well as the character X, which
# indicates that no amino acid will be present at the position. The 
# library object instance creates peptide objects based on 
# deconvolution of the template. Several attributes for each peptide
# are instantiated for use when analyzing MS/MS spectra.

import csv
from atomic_mass_data import *
from itertools import chain


class Library(object):
    

    def __init__(self, residue_csv):

        self.parent_sequence = self.parent_constructor(residue_csv) 
        self.peptide_list = self.library_builder(residue_csv)        
        
    def library_builder(self, residue_csv):
    # Residue data will be a tuple of lists. Each list in the tuple will
    # correspond to a site in the peptide and contain potential residues 
    # the Library instance will contain all sequences built from the
    # template by duplicating the library at a given site, adding each
    # potential residue to one of the copies, and rejoining the 
    # libraries after each iteration.

        master_library_string = [''] # hold library sequences as string
        master_library_index = [] # hold library as peptide objects

        with open(residue_csv) as csv_file:
            reader = csv.reader(csv_file, dialect='excel')

            for residues in reader:
                new_residues = [
                    residue for residue 
                    in residues 
                    if residue != '' 
                    ]

                temp_libs = self.master_replicator(
                    new_residues, master_library_string
                    )

                master_library_string = self.residue_appender(
                    new_residues, temp_libs
                    )

        for peptide_string in master_library_string:
            master_library_index.append(Peptide(peptide_string))

        master_library_index.sort(reverse=False)
        return master_library_index
        
    def master_replicator(self, residues, library):
        # This function will replicate the master_library the number of
        # times necessary based on the number of varied residues at a
        # given position

        temp_libs = []

        for i in range(len(residues)):
            temp_libs.append(library[:])

        return temp_libs
    
    def residue_appender(self, residues, temp_libs):
        # Adds residues to growing library for each copy created in 
        # at a given position.

        new_master = []

        for temp_copy, residue in zip(temp_libs, residues):
            for sequence in temp_copy:
                if residue != 'X':
                    sequence += residue
                new_master.append(sequence)

        return new_master # collapse temporary libs
    
    def parent_constructor(self, residue_CSV):
        # Builds an returns the sequence of the parent peptide by 
        # looking at the first value in each row and concatonating.

        parent_sequence = ''

        with open(residue_CSV) as csv_file:
            reader = csv.reader(csv_file, dialect='excel')
            for residues in reader:
                parent_sequence += residues[0]
        
        return parent_sequence
        

class Peptide(object):

    
    def __init__(self, sequence):

        self.sequence = sequence
        self.exact_mass = self.exact_mass_calculator()
        
        self.ion_series = {
            'a_ions':[],
            'b_ions':[],
            'c_ions':[],
            'x_ions':[],
            'y_ions':[],
            'z_ions':[]
            }

        self.match_dict = {} # Hold match info for each activation type
        for activation_method in ['etd', 'hcd', 'cid']:
            self.match_dict[activation_method] = {
                'match attempts':0, 
                'matches':0, 
                'matches/attempt':0
                }

        self.scan_ids = [] # MS/MS scans matches to peptide
    
    def __lt__(self, other):

        return self.exact_mass < other.exact_mass

    def residue_mass(self, sequence):
        # Handles consistent masses for each sequence. First calculates
        # backbone mass based on length of sequence, then appends each
        # residue mass.

        exact_mass = (len(sequence) - 1) * termini_map['bond']
        
        for amino_acid in sequence:
            exact_mass += amino_acid_map[amino_acid]
        
        return exact_mass

    def exact_mass_calculator(self):
        # Calculates the mw weight for a peptide represented as a 
        # string of amino acids

        exact_mass = self.residue_mass(self.sequence)
        exact_mass += termini_map['normal_n'] # add N-terminus mass
        exact_mass += termini_map['c_amide'] # add amidated C-terminus
        
        return exact_mass

    def generate_a_ions(self):
        # Calculates the mass of each a_ion and adds it to the list

        frags = [self.sequence[:a] for a in range(1,len(self.sequence)+1)]

        for subseq in frags:
            exact_mass = self.residue_mass(subseq)
            exact_mass += termini_map['normal_n']
            exact_mass += termini_map['a_terminus']
            self.ion_series['a_ions'].append(exact_mass)
        
        return

    def generate_b_ions(self):
        # Calculates mass of each b_ion

        frags = [self.sequence[:b] for b in range(1,len(self.sequence)+1)]

        for subseq in frags:
            exact_mass = self.residue_mass(subseq)
            exact_mass += termini_map['normal_n']
            exact_mass += termini_map['b_terminus']
            self.ion_series['b_ions'].append(exact_mass)
        
        return

    def generate_c_ions(self):
        # Calculates mass of each c_ion

        frags = [self.sequence[:c] for c in range(1,len(self.sequence)+1)]

        for subseq in frags:
            exact_mass = self.residue_mass(subseq)
            exact_mass += termini_map['normal_n']
            exact_mass += termini_map['c_terminus']
            self.ion_series['c_ions'].append(exact_mass)
        self.ion_series['c_ions'].pop(-1) # remove nonexistant c-ion
        
        return

    def generate_x_ions(self):
        # Calculates the mass of each x ion

        frags = [self.sequence[x:] for x in range(0,len(self.sequence))]
        
        for subseq in frags:
            exact_mass = self.residue_mass(subseq)
            exact_mass += termini_map['c_amide']
            exact_mass += termini_map['x_terminus']
            self.ion_series['x_ions'].append(exact_mass)
        self.ion_series['x_ions'].pop(0) # remove nonexistant x-ion
        
        return

    def generate_y_ions(self):
        # Calculates the mass of each y ion

        frags = [self.sequence[y:] for y in range(0,len(self.sequence))]
        
        for subseq in frags:
            exact_mass = self.residue_mass(subseq)
            exact_mass += termini_map['c_amide']
            exact_mass += termini_map['y_terminus']
            self.ion_series['y_ions'].append(exact_mass)
        
        return

    def generate_z_ions(self):
        # Calculates the mass of each z ion

        frags = [self.sequence[z:] for z in range(0,len(self.sequence))]
        for subseq in frags:
            exact_mass = self.residue_mass(subseq)
            exact_mass += termini_map['c_amide']
            exact_mass += termini_map['z_terminus']
            self.ion_series['z_ions'].append(exact_mass)
        
        return

    def calc_per_attempt(self):
        # Calculate the number of spectrum matches per match attempt

        for key in self.match_dict.iterkeys():
            if self.match_dict[key]['match attempts']:
                self.match_dict[key]['matches/attempt'] = (float(
                    self.match_dict[key]['matches'])
                    /self.match_dict[key]['match attempts']
                    )

        return


    