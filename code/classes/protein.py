import random
from .amino import Amino

class Protein():

    def __init__(self, source_file):
        self.aminos = self.load_protein(source_file)
        self.fold_protein()
        self.score = self.calculate_score()


    def load_protein(self, source_file):

        protein = []

        # open source file and read first line
        with open(source_file, 'r') as infile:
            protein_string = infile.readline()

        # transform string into list of chars
        amino_list = [element for element in protein_string]

        # iterate over every amino
        for i, amino_type in enumerate(amino_list):

            # initialize amino
            amino = Amino(i, amino_type)

            # add amino to protein-structure
            protein.append(amino)

        return protein


    def fold_protein(self):
        
        # set first coordinate to (0,0) and occupied fold to 0
        x = 0
        y = 0
        previous_amino = 0
        
        # iterate over every amino
        for amino in self.aminos:

            # set amino's coordinates
            amino.set_coordinate(x, y)

            # set amino's occupied fold
            amino.set_previous_amino(previous_amino)
            
            # generate and set fold
            fold = self.get_fold(x, y)
            amino.set_fold(fold)
                
            # compute next coordinate following the fold
            new_x, new_y = self.calculate_coordinate(fold, x, y)

            # set next coordinate values
            x = new_x
            y = new_y

            # set previous amino to inverse fold
            previous_amino = -fold


    def get_fold(self, x, y):
        
        # get possible folds
        possible_folds = self.get_possible_folds(x, y)

        # choose random fold
        # later version: optimal choice
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