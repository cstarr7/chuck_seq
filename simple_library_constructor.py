# -*- coding: utf-8 -*-
"""
Created on Mon Oct 05 15:50:24 2015

@author: Charlie
"""

#the purpose of this program is to create a combinatorial peptide library using
#a list of potentially varied residues

import csv
from atomic_mass_data import *
#global dictionaries that contain static mass information for sequence construction

class Library(object):
    
    def __init__(self,residue_CSV):
        self.parent_sequence = self.parentConstructor(residue_CSV) 
        self.peptide_list = self.libraryBuilder(residue_CSV)        
        
    def libraryBuilder(self,residue_CSV):
    #residue data will be a tuple of lists. each list in the tuple will contain
    #correspond to a site in the peptide and contain potential residues 
    #the masterLibrary will contain all sequences built from the library    
        master_library_string = ['']
        master_library_index = []
        with open(residue_CSV) as csv_file:
            reader = csv.reader(csv_file, dialect='excel')
            for residues in reader:
                new_residues = [residue for residue in residues if residue != '' ]
                #getting a copy number of the masterLibrary equal to the number of residues
                #at the current site        
                tempLibs = self.masterReplicator(new_residues, master_library_string)
                master_library_string = self.residueAppender(new_residues, tempLibs)
                #masterLibrary = itertools.chain(*updatedLibs)
        for peptide_string in master_library_string:
            master_library_index.append(Peptide(peptide_string))
        master_library_index.sort(reverse=False)
        return master_library_index
        
    def masterReplicator(self,residues,library):
        #this function will replicate the masterLibrary the number of times necessary
        #based on the number of varied residues at a given position
        listOfLists = []
        for i in range(len(residues)):
            listOfLists.append(library[:])
        return listOfLists
    
    def residueAppender(self,residues,tempLib):
        #adds residues to growing library    
        copyTemp = tempLib[:]
        newMaster = []
        for tempLists, residue in zip(copyTemp, residues):
            for sequence in tempLists:
                if residue != 'X':
                    sequence = sequence + residue
                newMaster.append(sequence)
        return newMaster
    
    def parentConstructor(self, residue_CSV):
        #builds an returns the sequence of the parent peptide
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
        self.ion_series ={'a_ions':[],
                        'b_ions':[],
                        'c_ions':[],
                        'x_ions':[],
                        'y_ions':[],
                        'z_ions':[]}
        self.match_attempts = 0
        self.matches = 0
        self.matches_per_attempt = 0
        self.scan_ids = []
    
    def __lt__(self, other):
        return self.exact_mass < other.exact_mass

    def residue_mass(self, subseq):
        #handles consistent masses for each sequence
        exact_mass = (len(subseq) - 1) * termini_map['bond'] #backbone weight (w/o termini)
        for amino_acid in subseq:
            exact_mass += amino_acid_map[amino_acid]
        return exact_mass

    def exact_mass_calculator(self):
        #calculates the mw weight for a peptide represented as a string of amino acids
        exact_mass = self.residue_mass(self.sequence)
        exact_mass += termini_map['normal_n'] #add N-terminus weight
        exact_mass += termini_map['c_amide'] #add amidated C-terminus weight
        return exact_mass

    def generate_a_ions(self):
        #calculates the mass of each a_ion and adds it to the list
        for subseq in [self.sequence[:a] for a in range(1,len(self.sequence)+1)]:
            exact_mass = self.residue_mass(subseq)
            exact_mass += termini_map['normal_n']
            exact_mass += termini_map['a_terminus']
            self.ion_series['a_ions'].append(exact_mass)
        return

    def generate_b_ions(self):
        #calculates mass of each b_ion
        for subseq in [self.sequence[:b] for b in range(1,len(self.sequence)+1)]:
            exact_mass = self.residue_mass(subseq)
            exact_mass += termini_map['normal_n']
            exact_mass += termini_map['b_terminus']
            self.ion_series['b_ions'].append(exact_mass)
        return

    def generate_c_ions(self):
        #calculates mass of each c_ion
        for subseq in [self.sequence[:c] for c in range(1,len(self.sequence)+1)]:
            exact_mass = self.residue_mass(subseq)
            exact_mass += termini_map['normal_n']
            exact_mass += termini_map['c_terminus']
            self.ion_series['c_ions'].append(exact_mass)
        self.ion_series['c_ions'].pop(-1)
        return

    def generate_x_ions(self):
        #calculates the mass of each x ion
        for subseq in [self.sequence[x:] for x in range(0,len(self.sequence))]:
            exact_mass = self.residue_mass(subseq)
            exact_mass += termini_map['c_amide']
            exact_mass += termini_map['x_terminus']
            self.ion_series['x_ions'].append(exact_mass)
        self.ion_series['x_ions'].pop(0)
        return

    def generate_y_ions(self):
        #calculates the mass of each y ion
        for subseq in [self.sequence[y:] for y in range(0,len(self.sequence))]:
            exact_mass = self.residue_mass(subseq)
            exact_mass += termini_map['c_amide']
            exact_mass += termini_map['y_terminus']
            self.ion_series['y_ions'].append(exact_mass)
        return

    def generate_z_ions(self):
        #calculates the mass of each z ion
        for subseq in [self.sequence[z:] for z in range(0,len(self.sequence))]:
            exact_mass = self.residue_mass(subseq)
            exact_mass += termini_map['c_amide']
            exact_mass += termini_map['z_terminus']
            self.ion_series['z_ions'].append(exact_mass)
        return

    def calc_per_attempt(self):
        if self.match_attempts:
            self.matches_per_attempt = float(self.matches)/self.match_attempts

'''
def gen_library(infile):
    library = Library(infile)
    return library

library = gen_library("arva_csv.csv")
for peptide in library.peptide_list:
    peptide.generate_a_ions()
    peptide.generate_b_ions()
    peptide.generate_c_ions()
    peptide.generate_x_ions()
    peptide.generate_y_ions()
    peptide.generate_z_ions()
    print peptide.sequence
    print peptide.a_series
    print peptide.b_series
    print peptide.c_series
    print peptide.x_series
    print peptide.y_series
    print peptide.z_series
'''
    