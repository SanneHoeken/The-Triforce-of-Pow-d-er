import random, math
from code import Amino, Protein, HillClimber

class SimulatedAnnealing(HillClimber):
    """
    The SimulatedAnnealing class that mutates a protein i times by changing n random amino folds.
    Each improvement is kept for the next iteration if the mutated protein is valid.
    Also sometimes accepts mutated configurations that are worse, depending on the current temperature.
    """
    def __init__(self, protein, iterations=1, mutations_per_iteration=1, temperature=1):
        super().__init__(protein, iterations, mutations_per_iteration)
        self.T0 = temperature
        self.T = temperature


    def run(self):
        """
        Runs the Simulated Annealing algorithm for a specific amount of iterations with
        a specified amount of mutations per iteration, and a specified starting temperature.
        """
        for i in range(self.iterations):

            # store values of current protein configuration
            self.archive = [(amino.fold, amino.coordinate, amino.previous_amino) for amino in self.protein.get_aminos()]
            
            # consecultively mutates protein specified times
            mutation_count = 0
            
            while mutation_count < self.mutations_per_iteration:
                
                # store values of current protein configuration
                self.tmp_archive = [(amino.fold, amino.coordinate, amino.previous_amino) for amino in self.protein.get_aminos()]
                
                # mutate protein
                self.mutate()

                # undo mutation if new configuration is not valid
                if not self.is_valid_protein():
                    self.undo_mutation(self.tmp_archive)
                else:
                    mutation_count += 1

            # calculate score of mutated configuration
            score = self.protein.calculate_score()

            # calculate probability of accepting the mutated configuration
            probability = math.exp((self.best_score - score) / self.T)
            
            # update best score if acceptance probability is higher than random probabiity
            # else undo mutation series
            if random.random() < probability:
                self.best_score = score
            else:
                self.undo_mutation(self.archive)

            # update temperature
            self.update_temperature()

        # set score of protein
        self.protein.set_score(self.best_score)


    def update_temperature(self):
        """
        This function implements a linear cooling scheme. Temperature will become zero 
        after all iterations passed to the run() method have passed.
        """
        self.T = self.T - (self.T0 / self.iterations)
