import random
from .amino import Amino
from .protein import Protein

# The goal of this class will be to take a protein and try different folding possibilities until it is satisfied witht he folding
class ProteinFolder():

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
        
        # iterate over every amino
        for amino in protein.get_aminos():

            # set amino's coordinates
            amino.set_coordinate(x, y)

            # set amino's occupied fold
            amino.set_previous_amino(previous_amino)
            
            # generate and set fold
            fold = protein.get_fold(x, y)
            if fold == 0:
                # iets waardoor onze amino niet verdergaat met deze state
                pass
             
            amino.set_fold(fold)
                
            # compute next coordinate following the fold
            new_x, new_y = protein.calculate_coordinate(fold, x, y)

            # set next coordinate values
            x = new_x
            y = new_y

            # set previous amino to inverse fold
            previous_amino = -fold

        



