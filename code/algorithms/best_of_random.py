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
            folder = RandomFolder(protein_copy)
            folder.fold_protein()
            score = protein_copy.get_score()

            # keep protein with lowest score
            if score < self.best_score:
                self.best_score = score
                self.best_protein = protein_copy
        
        # HIER GAAT NOG IETS MIS
        self.protein = self.best_protein
        self.protein.set_score(self.best_score)


class RandomFolder():

    def __init__(self, protein):
        self.protein = protein
        self.folding = 0


    def fold_protein(self):

        aminos = self.protein.get_aminos()

        while self.folding < len(aminos) - 1:
            
            current_amino = aminos[self.folding]
            next_amino = aminos[self.folding + 1]
            self.fold_amino(current_amino, next_amino)

        score = calculate_score(self.protein)
        self.protein.set_score(score)


    def fold_amino(self, current_amino, next_amino):

        new_values = self.get_values(current_amino)
        
        # stop folding if no values are returned
        if new_values is None:
            for amino in self.protein.get_aminos():
                amino.reset_amino()
            self.folding = 0

        fold, coordinate = new_values

        # set fold of current amino
        current_amino.set_fold(fold)

        # set coordinate and previous amino of next amino
        next_amino.set_coordinate(coordinate[0], coordinate[1])
        next_amino.set_previous_amino(-fold)

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