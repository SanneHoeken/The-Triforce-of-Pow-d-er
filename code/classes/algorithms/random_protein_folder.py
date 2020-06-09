import random
from code.classes.amino import Amino
from code.classes.protein import Protein
from code.classes.protein_folder import ProteinFolder

# The goal of this class will be to take a protein and try different random folding possibilities until it is satisfied with the folding
class RandomProteinFolder(ProteinFolder):

    def try_random(self, protein_to_fold=None, iterations=1):
        
        """
        TO CONTINUE
        """
        
        best_score = 1
        protein = protein_to_fold

        # If no argument is passed: uses origin protein
        if protein == None:
            assert self.origin_protein
            protein = self.origin_protein

        # fold protein iterations times 
        for i in range(iterations):
            self.fold(protein)
            self.calculate_score()
            score = protein.get_score()

            # keep protein with lowest score
            if score < best_score:
                best_score = score
                best_protein = Protein(protein=self.finished_folded_protein)
                iteration = i
        
        print(iteration)
        self.finished_folded_protein = best_protein


    def fold(self, protein_to_fold=None, fold_position = 0):

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
        if len(possible_folds) == 0: 
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


