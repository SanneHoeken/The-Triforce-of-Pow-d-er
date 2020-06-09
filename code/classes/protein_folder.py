import random
from .amino import Amino
from .protein import Protein
from .protein_tree import ProteinTree

# The goal of this class will be to take a protein and try different random folding possibilities until it is satisfied with the folding
class ProteinFolder():

    def __init__(self, protein = None):
        self.origin_protein = protein
        self.finished_folded_protein = None
        self.proteins_tree = None

    
    def set_origin_protein(self, protein):
        self.origin_protein = protein

    
    def fold(self, protein_to_fold = None, fold_position = 0):
        # This method does nothing, please override it in children classes!
        pass

    def calculate_score(self, protein = None):
        
        score = 0

        if protein is None:
            protein = self.finished_folded_protein

        assert protein is not None

        # checks for every H- or C- amino in protein if not-connected neighbor is H-amino
        for amino in protein.get_aminos():
            
            if amino.type == 'H' or amino.type == 'C':
                x, y = amino.coordinate

                # initializes not-connected directions
                free_folds = list({-2, 2, -1, 1} - {amino.previous_amino, amino.fold})
                
                # checks for every not-connected direction if H-amino is present
                for fold in free_folds:
                    next_coordinate = self.calculate_coordinate(fold, x, y)

                    # checks for hydrophobe/cysteine neighbor and changes score
                    neighbor = protein.get_aminotype(next_coordinate)
                    if amino.type == 'H' and neighbor == 'H':
                        score -= 1
                    if amino.type == 'H' and neighbor == 'C':
                        score -= 1
                    if amino.type == 'C' and neighbor == 'H':
                        score -= 1
                    if amino.type == 'C' and neighbor == 'C':
                        score -= 5            
                    

        protein.set_score(int(0.5 * score))
    
    
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