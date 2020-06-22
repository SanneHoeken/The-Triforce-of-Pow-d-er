import random
from .amino import Amino

class Protein():

    def __init__(self, dimensionality=2, **kwargs):
        self.source_file = kwargs.get('file', None)
        self.source_string = kwargs.get('string', None)
        self.source_protein = kwargs.get('protein', None)
        self.source_new_amino = kwargs.get('new_amino', None)
        self.D = dimensionality
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
                self.source_string = protein_string

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

            if self.D == 2:
                self.aminos[0].set_coordinate((0, 0))
            elif self.D == 3:
                self.aminos[0].set_coordinate((0, 0, 0))
                
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


    def to_string_with_previous(self):
        repr = ""
        for amino in self.aminos:
            repr += f" [{amino.fold} | prev: {amino.previous_amino}] {amino.type}"

        return repr


    def to_string_with_coord(self):
        repr = ""
        for amino in self.aminos:
            repr += f"{amino.type} ({amino.coordinate}) | "

        return repr


    def calculate_coordinate(self, fold, coordinate):
        """
        Method that takes the fold and coordinate of an amino
        Returns the coordinate of the next coordinate following the fold
        """
        co = list(coordinate)
        
        # adjust x-value following the fold
        if fold == 1 or fold == -1:
            co[0] += fold

        # adjust y-value following the fold    
        elif fold == 2:
            co[1] += 1

        elif fold == -2:
            co[1] -= 1

        # adjust z-value following the fold
        elif fold == 3:
            co[2] += 1
        
        elif fold == -3:
            co[2] -= 1
                
        return tuple(co)

    
    def calculate_score(self):
        """
        Returns a protein's stability score
        """ 
        score = 0

        # checks for every H- or C- amino in protein if not-connected neighbor is H-amino
        for amino in self.aminos:
            
            if amino.coordinate is not None:

                if amino.type == 'H' or amino.type == 'C':

                    # initializes not-connected directions if protein is 2D
                    if len(amino.coordinate) == 2:
                        free_folds = list({-2, 2, -1, 1} - {amino.previous_amino, amino.fold})
                    
                    # initializes not-connected directions if protein is 3D
                    elif len(amino.coordinate) == 3:
                        free_folds = list({-3, 3, -2, 2, -1, 1} - {amino.previous_amino, amino.fold})
                        
                    # checks for every not-connected direction if H- or C- amino is present
                    for fold in free_folds:
                        
                        # calculates coordinate of not-connected direction
                        next_coordinate = self.calculate_coordinate(fold, amino.coordinate)

                        # checks for hydrophobe/cysteine neighbor and changes score
                        neighbor = self.get_amino(next_coordinate)
                        if neighbor is not None:
                            if neighbor.type == 'P':
                                score -= 0
                            elif amino.type == 'H' and neighbor.type == 'H':
                                score -= 1
                            elif amino.type == 'H' and neighbor.type == 'C':
                                score -= 1
                            elif amino.type == 'C' and neighbor.type == 'H':
                                score -= 1
                            elif amino.type == 'C' and neighbor.type == 'C':
                                score -= 5            
                    
        return int(0.5 * score)