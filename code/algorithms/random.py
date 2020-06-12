import random
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
            random_folder = RandomFolder(self.protein)
            random_folder.fold_protein()
            score = random_folder.protein.get_score()

            # store protein values if score is lower
            if score < self.best_score:
                self.best_score = score
                self.best_protein = [(amino.fold, amino.coordinate) for amino in random_folder.protein.get_aminos()]

            # reset protein
            for amino in self.protein.get_aminos()[1:]:
                amino.reset_amino()

        # set values of input protein to best values
        for i, values in enumerate(self.best_protein):
            self.protein.get_aminos()[i].set_fold(values[0])
            self.protein.get_aminos()[i].set_coordinate(values[1][0], values[1][1])
            self.protein.get_aminos()[i].set_previous_amino(-values[0] if values[0] is not None else None)
        
        self.protein.set_score(self.best_score)


class RandomFolder():

    def __init__(self, protein):
        self.protein = protein
        self.coordinates = set()
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
                
                # set fold position back to zero and clear coordinates
                self.folding = 0
                self.coordinates = set()

        score = calculate_score(self.protein)
        self.protein.set_score(score)


    def fold_amino(self):

        # store current amino and its coordinate
        current_amino = self.protein.get_aminos()[self.folding]
        self.coordinates.add(current_amino.coordinate)

        # generate fold-values for current amino
        new_values = self.get_values(current_amino)
        
        # start over protein folding if no values are returned
        if new_values is None:
            return 1

        # else fold amino in protein according to the generated values
        self.set_values(new_values, self.folding)

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
                if self.is_free_space(coordinate):
                    possible_values.append((fold, coordinate))

        return possible_values


    def is_free_space(self, coordinate):

        # return True if coordinate is not occupied, else False
        return coordinate not in self.coordinates


    def set_values(self, values, id):

        # set fold of current amino
        self.protein.aminos[id].set_fold(values[0])

        # set coordinate en previous amino of next amino
        self.protein.aminos[id + 1].set_coordinate(values[1][0], values[1][1])
        self.protein.aminos[id + 1].set_previous_amino(-values[0])
        
        # return score of new configuration
        return calculate_score(self.protein)