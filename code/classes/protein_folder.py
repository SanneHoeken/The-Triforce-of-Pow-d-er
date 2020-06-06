import random
from .amino import Amino
from .protein import Protein

# The goal of this class will be to take a protein and try different folding possibilities until it is satisfied with the folding
class ProteinFolder(): # Note: Is it necessary to make this a class?

    def __init__(self, protein = None):
        self.origin_protein = protein
        self.proteins_tree = None # Todo: look for good existing pythons packages that implement trees in a simple / fast way
        self.finished_folded_protein = None

    
    def set_origin_protein(self, protein):
        self.origin_protein = protein

    # Here we can implement our optimization algorithm
    # For example, instead of just randomly folding one protein, we could make a loop (or whatever other way to create a tree) and
    # create a new protein every time, copy of the previous one but with one different folding, until we find a satisfying protein.
    # (That's just an example, we still need to decide what type of algorithm we use, how to prune etc. :))
    def fold(self, protein_to_fold = None, fold_position = 0):

        protein = protein_to_fold

        # If no argument is passed: uses origin protein
        if protein == None:
            assert self.origin_protein
            protein = self.origin_protein

        # set first coordinate to (0,0) and occupied fold to 0
        x = 0
        y = 0
        previous_amino = 0

        # repeat folding until all aminos are folded
        while fold_position < (len(protein.get_aminos()) - 2):    

            # iterate over every amino
            for amino in protein.get_aminos()[fold_position:]:

                # set amino's coordinates
                amino.set_coordinate(x, y)

                # set amino's occupied fold
                amino.set_previous_amino(previous_amino)
                
                # generate fold
                fold = self.get_fold(protein, x, y)
                
                # break if protein ran into dead end
                if fold == 0:
                    
                    # reset all amino values
                    for amino in protein.get_aminos()[:fold_position+1]:
                        amino.reset_amino()
                    
                    # reset fold position
                    fold_position = 0
                    break  
                
                # set fold and move fold position 
                amino.set_fold(fold)
                fold_position += 1
                    
                # compute next coordinate following the fold
                new_x, new_y = self.calculate_coordinate(fold, x, y)

                # set next coordinate values
                x = new_x
                y = new_y

                # set previous amino to inverse fold
                previous_amino = -fold
            
        self.finished_folded_protein = protein


    def get_fold(self, protein, x, y):
        
        # get possible folds
        possible_folds = self.get_possible_folds(protein, x, y)

        # choose random fold
        # later version: optimal choice
        if len(possible_folds) == 0:
            print("Protein folding resulted in dead end") 
            return 0
                
        fold = random.choice(possible_folds)
        
        return fold


    def get_possible_folds(self, protein, x, y):    
        possible_folds = []
        
        # dict with all neighbouring coordinates 
        coordinates = {1: (x+1, y), -1: (x-1, y), 2:(x, y+1), -2: (x, y-1)} 

        # add possible folds to list
        for key, value in coordinates.items(): 
            if self.is_free_space(protein, value):
                possible_folds.append(key)

        return possible_folds



    def is_free_space(self, protein, coordinate):
    
        # return True if coordinate is not occupied, else False
        return all([amino_object.coordinate != coordinate for amino_object in protein.get_aminos()])



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


    def calculate_score(self, protein = None):
        
        score = 0

        if protein is None:
            protein = self.finished_folded_protein

        assert protein is not None

        # checks for every H-amino in protein if not-connected neighbor is H-amino
        for amino in protein.get_aminos():
            
            if amino.type == 'H':
                x, y = amino.coordinate

                # initializes not-connected directions
                free_folds = list({-2, 2, -1, 1} - {amino.previous_amino, amino.fold})
                
                # checks for every not-connected direction if H-amino is present
                for fold in free_folds:
                    next_coordinate = self.calculate_coordinate(fold, x, y)

                    # checks for hydrophobe neighbor
                    neighbor = protein.get_aminotype(next_coordinate)
                    if neighbor == 'H':
                        score -= 1

        protein.set_score(0.5 * score)