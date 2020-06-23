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

            # Store values of current protein configuration
            self.archive = [(amino.get_fold(), amino.get_coordinate(), amino.get_previous_amino()) for amino in self.protein.get_aminos()]
            
            # Consecultively mutates protein specified times
            mutation_count = 0
            
            while mutation_count < self.mutations_per_iteration:
                
                # Store values of current protein configuration
                self.tmp_archive = [(amino.get_fold(), amino.get_coordinate(), amino.get_previous_amino()) for amino in self.protein.get_aminos()]
                
                # Mutate protein
                self.mutate()

                # Undo mutation if new configuration is not valid
                if not self.is_valid_protein():
                    self.undo_mutation(self.tmp_archive)
                else:
                    mutation_count += 1

            # Calculate score of mutated configuration
            score = self.protein.calculate_score()

            # Calculate probability of accepting the mutated configuration
            ### COMMENT TOEVOEGEN: even uitleggen dat dit een heuristiek is en wat de invloed hiervan op de probability is?
            probability = math.exp((self.best_score - score) / self.T)
            
            # Update best score if acceptance probability is higher than random probabiity
            # Else undo mutation series
            if random.random() < probability:
                self.best_score = score
            else:
                self.undo_mutation(self.archive)

            # Update temperature
            self.update_temperature()

        # Set score of protein
        self.protein.set_score(self.best_score)


    def update_temperature(self):
        """
        This function implements a linear cooling scheme. Temperature will become zero 
        after all iterations passed to the run() method have passed.
        """
        self.T = self.T - (self.T0 / self.iterations)
