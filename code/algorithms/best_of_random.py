import random
import copy
from code import Amino, Protein, calculate_score

class BestOfRandom():

    def __init__(self, protein, iterations=1):
        self.protein = protein
        self.iterations = iterations
        self.best_score = 1
        self.best_protein = None

    def run(self):

        # fold protein iterations times 
        for i in range(self.iterations):
        
            # make copy of protein, fold the copy, calculate score
            protein_copy = copy.deepcopy(self.protein)
            random_folder = RandomFolder(protein_copy)
            random_folder.fold_protein()
            score = random_folder.protein.get_score()

            # keep protein with lowest score
            if score < self.best_score:
                self.best_score = score
                self.best_protein = random_folder.protein
        
        # set input protein to best configuration
        self.protein = self.best_protein
        self.protein.set_score(self.best_score)


class RandomFolder():

    def __init__(self, protein):
        self.protein = protein
        self.folding = 0


    def fold_protein(self):

        length = len(self.protein.get_aminos())

        # keep folding until protein is completely folded
        while self.folding < length - 1:
            
            # fold amino and start over if amino folding failed
            if self.fold_amino() == 1:
                
                for amino in self.protein.get_aminos()[1:]:
                    
                    # reset all values of aminos
                    amino.reset_amino()
                
                # set fold position back to zero
                self.folding = 0

        score = calculate_score(self.protein)
        self.protein.set_score(score)


    def fold_amino(self):

        # generate fold-values for current amino
        current_amino = self.protein.get_aminos()[self.folding]
        new_values = self.get_values(current_amino)
        
        # start over protein folding if no values are returned
        if new_values is None:
            return 1

        # else fold amino in protein according to the generated values
        amino_folder = FoldAmino(self.protein, self.folding, new_values)
        amino_folder.fold_amino_in_protein()
        self.protein = copy.deepcopy(amino_folder.protein)

        # update fold position
        self.folding += 1


    def get_values(self, amino):
        
        possible_values = self.get_possible_values(amino)
        
        # return None if no possible values
        if len(possible_values) == 0:
            return None
        
        # choose random values
        fold, coordinate = random.choice(possible_values)
        
        return fold, coordinate

    
    def get_possible_values(self, amino):    
        
        possible_values = []

        if amino.coordinate:
            x, y = amino.coordinate
        
            # dict with all possible fold-coordinate combinations 
            values = {1: (x+1, y), -1: (x-1, y), 2:(x, y+1), -2: (x, y-1)} 

            # add possible folds to list
            for fold, coordinate in values.items(): 
                if self.is_free_space(self.protein, coordinate):
                    possible_values.append((fold, coordinate))

        return possible_values


    def is_free_space(self, protein, coordinate):

        # return True if coordinate is not occupied, else False
        return all([amino_object.coordinate != coordinate for amino_object in protein.get_aminos()])


class FoldAmino():

    def __init__(self, protein, position, values):
        self.protein = copy.deepcopy(protein)
        self.fold, self.coordinate = values
        self.position = position
        self.score = 1
    
    def fold_amino_in_protein(self):

        # get amino to fold and its next neighbor
        current_amino = self.protein.aminos[self.position]
        next_amino = self.protein.aminos[self.position + 1]

        # set fold of current amino
        current_amino.set_fold(self.fold)

        # set coordinate and previous amino of next amino
        next_amino.set_coordinate(self.coordinate[0], self.coordinate[1])
        next_amino.set_previous_amino(-self.fold)

        # set score of new protein configuration
        self.score = calculate_score(self.protein)