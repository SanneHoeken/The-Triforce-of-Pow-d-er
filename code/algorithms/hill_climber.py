import random
from code import Amino, Protein, calculate_score, calculate_coordinate

class HillClimber():

    def __init__(self, protein, iterations=1):
        self.protein = protein
        self.protein_length = len(self.protein.get_aminos())
        self.iterations = iterations
        self.tmp_archive = []
        self.best_score = 1
    

    def run(self):
        
        for i in range(self.iterations):

            # store values of current protein configuration
            self.tmp_archive = [(amino.fold, amino.coordinate) for amino in self.protein.get_aminos()]
            
            # mutate protein
            self.mutate()

            # undo mutation if new configuration is not valid
            if not self.is_valid_protein():
                self.undo_mutation()

            # update best score if new configuration is better, else undo mutation
            score = calculate_score(self.protein)
            if score < self.best_score:
                self.best_score = score
            else:
                self.undo_mutation()


    def mutate(self):

        # get mutation values: position of amino and new fold
        pos, fold = self.generate_mutation()
        
        x, y = (None, None)

        for amino in self.protein.get_aminos()[pos: self.protein_length]:
            
             # change fold for selected amino
            if amino.id == pos:
                amino.set_fold(fold)

                # update coordinate
                x, y = amino.coordinate

             # change coordinates and previous amino-fold for amino after selected amino
            elif amino.id == pos + 1:
                new_x, new_y = calculate_coordinate(fold, x, y)
                amino.set_previous_amino(-fold)
                amino.set_coordinate(new_x, new_y)

                # update fold and coordinate
                fold = amino.fold
                x, y = new_x, new_y

            # change coordinates for every amino after selected amino
            else:
                new_x, new_y = calculate_coordinate(fold, x, y)
                amino.set_coordinate(new_x, new_y)

                # update fold and coordinate
                fold = amino.fold
                x, y = new_x, new_y
    

    def generate_mutation(self):

        # choose random position
        position = random.randint(0, self.protein_length - 1)

        # choose random fold
        new_fold = random.choice([-2, 2, -1, 1])
        return position, new_fold


    def is_valid_protein(self):

        # returns Ture if no double coordinates in protein
        coordinates = [amino.coordinate for amino in self.protein.get_aminos()]
        
        return len(set(coordinates)) == len(coordinates)


    def undo_mutation(self):
        
        # changes all values of protein to values in temporary archive
        for i, values in enumerate(self.tmp_archive):
            self.protein.get_aminos()[i].set_fold(values[0])
            self.protein.get_aminos()[i].set_coordinate(values[1][0], values[1][1])
            self.protein.get_aminos()[i].set_previous_amino(-values[0] if values[0] is not None else None)
