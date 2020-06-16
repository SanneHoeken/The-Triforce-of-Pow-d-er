import random
from code import calculate_score, calculate_coordinate, Amino, Protein

class HillClimber():

    def __init__(self, protein, iterations=1, mutations_per_iteration=1):
        self.protein = protein
        self.protein_length = len(self.protein.get_aminos())
        self.iterations = iterations
        self.mutations_per_iteration = mutations_per_iteration
        self.archive = []
        self.tmp_archive = []
        self.best_score = self.protein.get_score()
    

    def run(self):
        
        for i in range(self.iterations):

            # store values of current protein configuration
            self.archive = [(amino.fold, amino.coordinate) for amino in self.protein.get_aminos()]
            
            for j in range(self.mutations_per_iteration):

                # store values of current protein configuration
                self.tmp_archive = [(amino.fold, amino.coordinate) for amino in self.protein.get_aminos()]
                
                # mutate protein
                self.mutate()

                # undo mutation if new configuration is not valid
                if not self.is_valid_protein():
                    self.undo_mutation(self.tmp_archive)

            # update best score if new configuration is better, else undo mutation series
            score = calculate_score(self.protein)
            if score < self.best_score:
                self.best_score = score
            else:
                self.undo_mutation(self.archive)
            
        # set score of protein
        self.protein.set_score(self.best_score)


    def mutate(self):

        # get mutation values: position of amino and new fold
        pos, fold = self.generate_mutation()
        
        co = (None, None)

        for amino in self.protein.get_aminos()[pos: self.protein_length]:
            
            # change fold for selected amino
            if amino.id == pos:
                amino.set_fold(fold)

                # update coordinate
                co = amino.coordinate

            # change coordinates and previous amino-fold for amino after selected amino
            elif amino.id == pos + 1:
                new_co = calculate_coordinate(fold, co)
                amino.set_previous_amino(-fold)
                amino.set_coordinate(new_co)

                # update fold and coordinate
                fold = amino.fold
                co = new_co

            # change coordinates for every amino after selected amino
            else:
                new_co = calculate_coordinate(fold, co)
                amino.set_coordinate(new_co)

                # update fold and coordinate
                fold = amino.fold
                co = new_co
    

    def generate_mutation(self):

        # choose random position
        position = random.randint(0, self.protein_length - 1)

        # choose random fold
        if len(self.protein.get_aminos()[0].coordinate) == 2: 
            new_fold = random.choice([-2, 2, -1, 1])

        if len(self.protein.get_aminos()[0].coordinate) == 3: 
            new_fold = random.choice([-3, 3, -2, 2, -1, 1])
        
        return position, new_fold


    def is_valid_protein(self):

        # returns True if no double coordinates in protein
        coordinates = [amino.coordinate for amino in self.protein.get_aminos()]
        
        return len(set(coordinates)) == len(coordinates)


    def undo_mutation(self, previous_configuration):
        
        # changes all values of protein to values of previous configuration
        for i, values in enumerate(previous_configuration):
            self.protein.get_aminos()[i].set_fold(values[0])
            self.protein.get_aminos()[i].set_coordinate(values[1])
            self.protein.get_aminos()[i].set_previous_amino(-values[0] if values[0] is not None else None)
