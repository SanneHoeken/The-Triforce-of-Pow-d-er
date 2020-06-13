import random, math
from code import Amino, Protein, calculate_score, calculate_coordinate
from code.algorithms.hill_climber import HillClimber

class SimulatedAnnealing(HillClimber):

    def __init__(self, protein, iterations=1, mutations_per_iteration=1, temperature=1):

        super().__init__(protein, iterations, mutations_per_iteration)
        self.T0 = temperature
        self.T = temperature


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

            
            score = calculate_score(self.protein)
            probability = math.exp((self.best_score - score) / self.T)

            # update best score if acceptance probability is higher than random probabiity
            # else undo mutation series
            if random.random() < probability:
                self.best_score = score
            else:
                self.undo_mutation(self.archive)

            self.update_temperature()
            
        # set score of protein
        self.protein.set_score(self.best_score)


    def update_temperature(self):
        
        # linear cooling
        self.T = self.T - (self.T0 / self.iterations)
