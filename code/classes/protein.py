import random
from .amino import Amino

class Protein():

    def __init__(self, **kwargs):
        self.source_file = kwargs.get('file', None)
        self.source_string = kwargs.get('string', None)
        self.source_protein = kwargs.get('protein', None)
        self.source_new_amino = kwargs.get('new_amino', None)
        self.aminos = []
        self.load_aminos()
        self.set_first_amino()
        self.score = 0
 

    def load_aminos(self):

        # Loading amino acids from existing protein
        if self.source_protein:

            # No new fold, just copying another protein
            if self.source_new_amino is None:
                self.aminos = self.source_protein.get_aminos()
                return self.aminos

            # Copying a protein + adding an amino acid
            else:
                self.aminos = self.source_protein.get_aminos()
                self.add_amino(self.source_new_amino)
                return self.aminos

        # Making the protein from a string
        if self.source_file:

            # open source file and read first line
            with open(self.source_file, 'r') as infile:
                protein_string = infile.readline()

        elif self.source_string:
            protein_string = self.source_string
    
        assert protein_string is not None

        # transform string into list of chars
        amino_list = [element for element in protein_string]

        # iterate over every amino
        for amino_type in amino_list:
            self.add_amino(amino_type)
                
        return self.aminos


    def set_first_amino(self):
        
        # set coordinate and previous amino of first amino to zero
        if self.aminos[0] is not None:
            self.aminos[0].set_coordinate(0, 0)
            self.aminos[0].set_previous_amino(0)

        return self.aminos


    def get_score(self):
        return self.score


    def set_score(self, score):
        self.score = score
    

    def get_aminos(self):
        return self.aminos


    def add_amino(self, amino_type):
        # initialize amino
        amino = Amino(len(self.aminos), amino_type)

        # add amino to protein-structure
        self.aminos.append(amino)


    def get_amino(self, coordinate):
        
        # returns amino that is present on a certain coordinate
        for amino in self.aminos:
            if amino.coordinate == coordinate:
                return amino
                
        return None


    def to_string(self):
        repr = ""
        for amino in self.aminos:
            repr += f" [{amino.fold}] {amino.type}"

        return repr