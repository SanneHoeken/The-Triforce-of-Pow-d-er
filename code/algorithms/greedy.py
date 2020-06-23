import random
from code import Amino, Protein
from code.algorithms.random import RandomFolder

class BestGreedy():
    """
    The BestGreedy Class that greedily folds a protein i times and saves the best configuration
    """
    def __init__(self, protein, iterations=1):
        self.protein = protein
        self.iterations = iterations
        self.best_score = 1
        self.best_protein = None
        self.finished_folded_protein = self.protein

    def run(self):
        """
        Folds protein greedily for every iteration, stores its values if configuration has a lower score,
        and then resets protein. After the iterations, the best values are set in the protein.
        """
        # fold protein iterations times 
        for i in range(self.iterations):
            greedy_folder = Greedy(self.protein)
            greedy_folder.fold_protein()

            # calculate score
            score = greedy_folder.protein.get_score()

            # update protein values en score if score is lower
            if score < self.best_score:
                self.best_score = score
                self.best_protein = [(amino.fold, amino.coordinate, amino.previous_amino) for amino in greedy_folder.protein.get_aminos()]

            # reset protein
            for amino in self.protein.get_aminos()[1:]:
                amino.reset_amino()

        # set values of protein to best values
        for i, values in enumerate(self.best_protein):
            self.protein.get_aminos()[i].set_fold(values[0])
            self.protein.get_aminos()[i].set_coordinate(values[1])
            self.protein.get_aminos()[i].set_previous_amino(values[2])
        
        self.protein.set_score(self.best_score)    


class Greedy(RandomFolder):
    """
    The Greedy Class that greedily assigns possible folding values to every amino one by one.
    Starts over if folding results in dead end (= no more available values for an amino to fold).
    """
    def get_values(self, amino):
        """
        Retrieves values for amino by choosing the best values out of possible values
        """
        # retrieve all possible values for amino
        possible_values = self.get_possible_values(amino)  
        
        # return None if no possible values
        if len(possible_values) == 0:
            return None

        best_values = []
        best_score = 1

        # get best of possible values
        for values in possible_values:
    
            # fold amino according to fold values
            # store score of new folded configuration
            score = self.set_values(values, self.folding)

            # replace best values if score is lower
            if score < best_score:
                best_score = score
                best_values = [values]
            
            # add best values if score is as low as lowest
            elif score == best_score:
                best_values.append(values)

            # undo amino folding
            self.reset_values(self.folding)       
        
        # return one of best values (randomly chosen)
        return random.choice(best_values)


    def reset_values(self, id):
        """
        Resets fold of amino with specified id, resets coordinate and 
        previous fold of that amino's next neighbor 
        """
        # reset fold of current amino
        self.protein.aminos[id].set_fold(None)

        # reset coordinate and previous amino of next amino
        self.protein.aminos[id + 1].reset_amino()
