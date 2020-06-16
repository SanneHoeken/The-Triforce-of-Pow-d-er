import random
from code import calculate_score, Amino, Protein
from code.algorithms.random import RandomFolder

class BestGreedy():

    def __init__(self, protein, iterations=1):
        self.protein = protein
        self.iterations = iterations
        self.best_score = 1
        self.best_protein = None

    def run(self):

        # fold protein iterations times 
        for i in range(self.iterations):
            greedy_folder = Greedy(self.protein)
            greedy_folder.fold_protein()
            score = greedy_folder.protein.get_score()

            # store protein values if score is lower
            if score < self.best_score:
                self.best_score = score
                self.best_protein = [(amino.fold, amino.coordinate) for amino in greedy_folder.protein.get_aminos()]

            # reset protein
            for amino in self.protein.get_aminos()[1:]:
                amino.reset_amino()

        # set values of input protein to best values
        for i, values in enumerate(self.best_protein):
            self.protein.get_aminos()[i].set_fold(values[0])
            self.protein.get_aminos()[i].set_coordinate(values[1][0], values[1][1])
            self.protein.get_aminos()[i].set_previous_amino(-values[0] if values[0] is not None else None)
        
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
            score = self.set_values(values, self.folding)

            # store values if score is better/best
            if score < best_score:
                best_score = score
                best_values = [values]
            elif score == best_score:
                best_values.append(values)

            # undo amino folding
            self.reset_values(self.folding)       
        
        # return one of best values
        return random.choice(best_values)


    def reset_values(self, id):

        # reset fold of current amino
        self.protein.aminos[id].set_fold(None)

        # reset coordinate and previous amino of next amino
        self.protein.aminos[id + 1].reset_amino()
