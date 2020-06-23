import random
from code import Amino, Protein
from code.algorithms.greedy import Greedy

class BranchBound(Greedy):

    def __init__(self, protein, p1=0.8, p2=0.5):
        self.protein = protein
        self.p1 = p1
        self.p2 = p2
        self.all_scores = {}
        self.best_scores = {}
        self.best_protein = None
        self.best_score = 1


    def run(self):
        self.search(0)
        
        # Set best protein values to protein
        for i, values in enumerate(self.best_protein):
            self.protein.get_aminos()[i].set_fold(values[0])
            self.protein.get_aminos()[i].set_coordinate(values[1])
            self.protein.get_aminos()[i].set_previous_amino(values[2])
        
        # Set best score to protein
        self.protein.set_score(self.best_score)


    def search(self, id):
        
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
            
            # Reset values of last amino folding
            self.reset_values(id)
            
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
                
                # Continue searching if score is as good as/better than best score
                bestscore = self.best_scores[id]
                if score <= bestscore:
                    self.best_scores[id] = score
                    self.set_values(values, id)
                    self.search(id + 1)
                
                # Continue with probility p1 if score is higher than meanscore
                meanscore = sum(self.all_scores[id]) / len(self.all_scores[id])
                if score > meanscore: 
                    if random.random() > self.p1: 
                        self.set_values(values, id)
                        self.search(id + 1)

                # Continue with probability p2 if score is higher than best score but lower than meanscore
                if score > bestscore and score <= meanscore: 
                    if random.random() > self.p2:
                        self.set_values(values, id)
                        self.search(id + 1)
            
            # Continue searching if amino is P
            else:
                self.set_values(values, id)
                self.search(id + 1)


    def is_free_space(self, coordinate):
        """
        Checks whether coordinate is already in proteins occupied coordinates
        """
        # Return True if coordinate is not occupied, else False
        return coordinate not in [amino.get_coordinate() for amino in self.protein.get_aminos()]
