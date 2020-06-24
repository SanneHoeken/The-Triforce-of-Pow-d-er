import random
from code import Amino, Protein
from code.algorithms.greedy import Greedy

"""
This algorithm is based on the work of:
Chen, M., & Huang, W. Q. (2005). A branch and bound algorithm for the protein folding problem in the HP lattice model. Genomics, proteomics & bioinformatics, 3(4), 225-230.
"""

class BranchBound(Greedy):
    """
    The BranchBound class that recursively folds a protein by searching for the best possibilities. 
    In the search method, the potency of the partial configuration is assessed for each H- and C-amino. 
    A negative evaluation leads to pruning of that configuration. No pruning is done at P-aminos. 
    The potential of a configuration is determined by two variables: the average score of a configuration of 
    a certain length and the best score for a configuration of a certain length. The score of a configuration 
    is compared with these two variables. If the score is better than the best score, the partial configuration 
    is not pruned. A score worse than the average score is pruned with a given probability p1. A score better 
    than the average value but worse than the best score is pruned with a probability p2.
    """
    def __init__(self, protein, p1=0.8, p2=0.5):
        self.protein = protein
        self.p1 = p1
        self.p2 = p2
        self.all_scores = {}
        self.best_scores = {}
        self.best_protein = None
        self.best_score = 1
        self.finished_folded_protein = self.protein


    def run(self):
        """
        Runs the search method. After searching, the best values and score are set in the protein.
        """
        self.search(0)
        
        # Set best protein values to protein
        for i, values in enumerate(self.best_protein):
            self.protein.get_aminos()[i].set_fold(values[0])
            self.protein.get_aminos()[i].set_coordinate(values[1])
            self.protein.get_aminos()[i].set_previous_amino(values[2])
        
        # Set best score to protein
        self.protein.set_score(self.best_score)


    def search(self, id):
        """
        Searching depht first by recursively visiting every configuration (except for the pruned configurations)
        """
        
        # Get possible values for current amino (specified with given id)
        possible_values = self.get_possible_values(self.protein.get_aminos()[id])
        
        # Iterate over all possible values 
        for values in possible_values:
            
            # Set values to protein
            self.set_values(values, id)

            # Calculate score of new protein configuration
            score = self.protein.calculate_score()

            # If last amino is reached, stores best protein, continues loop
            if id == len(self.protein.get_aminos()) - 2:
                if score < self.best_score:
                    self.best_score = score
                    self.best_protein = [(amino.get_fold(), amino.get_coordinate(), amino.get_previous_amino()) for amino in self.protein.get_aminos()]
                continue
            
            # Add score to all scores in current depth
            if id not in self.all_scores:
                self.all_scores[id] = []
            
            self.all_scores[id].append(score)
            
            # If current depth not in best_score dict, add it
            if id not in self.best_scores:
                self.best_scores[id] = 1
            
            # Store current amino
            amino = self.protein.get_aminos()[id]
            
            # Check for pruning if amino is H or C
            if amino.type == 'H' or amino.type == 'C':
                
                bestscore = self.best_scores[id]
                meanscore = sum(self.all_scores[id]) / len(self.all_scores[id])

                # Continue searching if score is as good as/better than best score
                if score <= bestscore:
                    self.best_scores[id] = score
                    self.search(id + 1)
                
                # Continue with probility 1 - p1 if score is higher than meanscore
                elif score > meanscore: 
                    if random.random() > self.p1: 
                        self.search(id + 1)

                # Continue with probability 1 - p2 if score is higher than best score but lower than meanscore
                elif score > bestscore and score <= meanscore: 
                    if random.random() > self.p2:
                        self.search(id + 1)

                else:
                    # Reset values of last amino folding
                    self.reset_values(id)
            
            # Continue searching if amino is P
            else:
                self.search(id + 1)


    def is_free_space(self, coordinate):
        """
        Checks whether coordinate is already in proteins occupied coordinates
        """
        # Return True if coordinate is not occupied, else False
        return coordinate not in [amino.get_coordinate() for amino in self.protein.get_aminos()]
