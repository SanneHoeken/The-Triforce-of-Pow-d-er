import copy
from code import Protein, Amino

class BruteForce():
    """
    The BruteForce Class that searches depht first for the best protein configuration
    """
    def __init__(self, protein):
        self.protein = protein
        self.best_score = 1
        self.best_state = None
        self.symmetry = set()
    
    
    def run(self):
        """
        Runs the depth first algorithm until all possible states are visited.
        Folds protein according to best state and sets protein's score.
        WARNING: runtime of long proteins (>14 aminos) is really big!
        """
        self.depth_first()
        self.protein_score(self.best_state)
        self.protein.set_score(self.best_score)


    def depth_first(self):
        """
        A Depth First algorithm that builds a stack of lists of folds with a unique assignment of folds for each instance.
        Checks for every state (list of folds) in the maximum depth if the protein configuration according to that state
        is better.  
        """
        # initialize depth and stack
        depth = len(self.protein.get_aminos()) - 1
        stack = [[]]
        
        while len(stack) > 0:

            # get next state from the list of states
            state = stack.pop()

            # check state if state is in maximum depth
            if len(state) == depth:
                
                # fold protein according to list of folds
                # store score of new folded protein
                score = self.protein_score(state)

                # update best state and best score of score is lower
                if score < self.best_score:
                    self.best_score = copy.deepcopy(score)
                    self.best_state = copy.deepcopy(state)
                
                # reset protein
                for amino in self.protein.get_aminos()[1:]:
                    amino.reset_amino()

            # create all possible child-states and adds them to list of states 
            if len(state) < depth:
                for i in [-1, 1, -2, 2]:
                    child = copy.deepcopy(state)
                    child.append(i)
                    stack.append(child)


    def protein_score(self, folds):
        """
        Method that takes protein and list of folds and returns score of folded protein
        """
        # set initial coordinate and previous amino
        coordinate = self.protein.get_aminos()[0].get_coordinate()
        previous_amino = self.protein.get_aminos()[0].get_previous_amino()

        # iterate over every amino in protein
        for amino in self.protein.get_aminos():
            
            # set amino's coordinates
            amino.set_coordinate(coordinate)

            # return 1 if new amino position is not valid
            if not self.is_valid_protein():
                return 1

            # set amino's occupied fold
            amino.set_previous_amino(previous_amino)
            
            if amino.id < len(folds):
                # set fold and move fold position 
                fold = folds[amino.id]
                amino.set_fold(fold)
                    
                # compute next coordinate following the fold
                new_coordinate = self.protein.calculate_coordinate(fold, coordinate)

                # set next coordinate values
                coordinate = new_coordinate

                # set previous amino to inverse fold
                previous_amino = -fold

        # calculate score of folded protein
        score = self.protein.calculate_score()
        
        return score


    def is_valid_protein(self):
        """
        Checks whether protein configuration is valid
        """
        # returns True if no double coordinates in protein
        coordinates = [amino.get_coordinate() for amino in self.protein.get_aminos() if amino.get_coordinate() is not None]
        
        return len(set(coordinates)) == len(coordinates)