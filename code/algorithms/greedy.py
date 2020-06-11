import random
import copy
from code import Amino, Protein, RandomFolder, FoldAmino, calculate_score

class BestGreedy():

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
            greedy_folder = Greedy(protein_copy)
            greedy_folder.fold_protein()
            score = greedy_folder.protein.get_score()

            # keep protein with lowest score
            if score < self.best_score:
                self.best_score = score
                self.best_protein = greedy_folder.protein
        
        # set input protein to best configuration
        self.protein = self.best_protein
        self.protein.set_score(self.best_score)    


class Greedy(RandomFolder):

    def get_values(self, amino):
        
        possible_values = self.get_possible_values(amino)  
        
        # return None if no possible values
        if len(possible_values) == 0:
            return None

        best_values = []
        best_score = 1

        # get best of possible values
        for values in possible_values:

            # fold amino according to fold values
            amino_folder = FoldAmino(self.protein, self.folding, values)
            amino_folder.fold_amino_in_protein()
            
            # store values if score is better/best
            score = amino_folder.score
            if score < best_score:
                best_score = score
                best_values = [values]
            elif score == best_score:
                best_values.append(values)
        
        # return one of best values
        return random.choice(best_values)
