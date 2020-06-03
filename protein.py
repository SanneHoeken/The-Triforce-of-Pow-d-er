import random
from amino import Amino

class Protein():

    def __init__(self, amino_list):
        self.aminos = amino_list
        self.folds = []
        self.coordinates = []
        self.score = 0
        self.load_structure(amino_list)

    def load_structure(self, amino_list):
        
        x = 0
        y = 0
        
        for i, amino in enumerate(amino_list):
            
            amino = Amino(i, amino)
            
            self.coordinates.append((x, y))
            amino.set_coordinate(x, y)
            
            fold = random.randint(-2, 3)
            self.folds.append(fold)
            amino.set_fold = fold

            if fold == 1 or fold == -1:
                x += fold
            
            if fold == 2 or fold == -2:
                y += 0.5 * fold

    def set_score(self):
        pass