import random
from .amino import Amino

class Protein():

    # Note: Op dit moment moeten de aminos meteen gedurende het aanmaak van de eiwit worden geladen. Dat betekent dat we geen mogelijkheid
    # hebben om een nieuwe eiwit aan te maken en pas later aminos te laden. Ik weet niet of het uitmaakt, misschien niet, maar
    # ik wilde het wel opmerken.
    def __init__(self, **kwargs):
        self.source_file = kwargs.get('file', None)
        self.source_string = kwargs.get('string', None)
        self.source_protein = kwargs.get('protein', None)
        self.source_new_amino = kwargs.get('new_amino', None)
        self.aminos = []
        self.load_aminos()
        self.score = 0
 
 
    def get_score(self):
        return self.score
    
    def get_aminos(self):
        return self.aminos

    def load_aminos(self):

        # Loading amino acids from existing protein
        if self.source_protein:

            # No new fold, just copying another protein
            if self.source_new_amino is None:
                self.aminos = self.source_protein.get_aminos
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



    def add_amino(self, amino_type):
        # initialize amino
        amino = Amino(len(self.aminos), amino_type)

        # add amino to protein-structure
        self.aminos.append(amino)



    def get_fold(self, x, y):
        
        # get possible folds
        possible_folds = self.get_possible_folds(x, y)

        # choose random fold
        # later version: optimal choice
        if len(possible_folds) == 0:
            print("Protein folding resulted in dead end") 
            return 0
                
        fold = random.choice(possible_folds)
        
        return fold



    def get_possible_folds(self, x, y):
        
        possible_folds = []
        
        # dict with all neighbouring coordinates 
        coordinates = {1: (x+1, y), -1: (x-1, y), 2:(x, y+1), -2: (x, y-1)} 

        # add possible folds to list
        for key, value in coordinates.items(): 
            if self.is_not_occupied(value):
                possible_folds.append(key)

        return possible_folds


    def is_not_occupied(self, coordinate):
    
        # return True if coordinate is not occupied, else False
        return all([amino_object.coordinate != coordinate for amino_object in self.aminos])


    def calculate_coordinate(self, fold, x, y):
     
        # compute x-value following the fold
        if fold == 1 or fold == -1:
            x_tmp = x + fold

        # keep current x-value if no move along x-axis 
        else:
            x_tmp = x

        # compute y-value following the fold    
        if fold == 2 or fold == -2:
            y_tmp = y + 0.5 * fold
        
        # keep current y-value if no move along y-axis
        else:
            y_tmp = y
        
        return int(x_tmp), int(y_tmp)


    def calculate_score(self):
        
        score = 0

        # checks for every H-amino in protein if not-connected neighbor is H-amino
        for amino in self.aminos:
            
            if amino.type == 'H':
                x, y = amino.coordinate

                # initializes not-connected directions
                free_folds = list({-2, 2, -1, 1} - {amino.previous_amino, amino.fold})
                
                # checks for every not-connected direction if H-amino is present
                for fold in free_folds:
                    next_coordinate = self.calculate_coordinate(fold, x, y)

                    # checks for hydrophobe neighbor
                    neighbor = self.get_aminotype(next_coordinate)
                    if neighbor == 'H':
                        score -= 1

        return 0.5 * score


    def get_aminotype(self, coordinate):
        
        # returns amino that is present on a certain coordinate
        for amino in self.aminos:
            if amino.coordinate == coordinate:
                return amino.type
                
        return None