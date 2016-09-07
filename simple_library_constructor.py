# -*- coding: utf-8 -*-
"""
Created on Mon Oct 05 15:50:24 2015

@author: Charlie
"""

#the purpose of this program is to create a combinatorial peptide library using
#a list of potentially varied residues

import csv

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
        master_library_index.sort(reverse=True)
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
    
    def __lt__(self, other):
        return self.exact_mass < other.exact_mass
        

    def exact_mass_calculator(self):
        #calculates the mw weight for a peptide represented as a string of amino acids
        peptide_exact_mass = 0
        amino_acid_index = {'R': 100.087472,
                            'Q': 72.044939,
                            'F': 91.054775,
                            'Y': 107.04969,
                            'W': 130.065674,
                            'K': 72.081324,
                            'G': 1.007825,
                            'A': 15.023475,
                            'H': 81.045273,
                            'S': 31.01839,
                            'P': 41.039125,
                            'E': 73.028955,
                            'D': 59.013305,
                            'T': 45.03404,
                            'C': 46.994185,
                            'M': 75.025485,
                            'L': 57.070425,
                            'N': 58.029289,
                            'I': 57.070425,
                            'V': 43.054775,
                            'Z': 98.081324}
        peptide_exact_mass += (len(self.sequence) - 1) * 56.013639 #backbone weight (w/o termini)
        peptide_exact_mass += 16.018724 #add N-terminus weight
        peptide_exact_mass += 57.021464 #add amidated C-terminus weight
        for amino_acid in self.sequence:
            peptide_exact_mass += amino_acid_index[amino_acid]
        return peptide_exact_mass

def gen_library(infile):
    library = Library(infile)
    return library

library = gen_library("arva_csv.csv")
for peptide in library.peptide_list:
    print peptide.sequence
    print peptide.exact_mass
    