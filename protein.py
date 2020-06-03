import random
from amino import Amino

class Protein():

    def __init__(self, amino_list):
        self.aminos = amino_list
        self.length = len(amino_list)
        self.folds = []
        self.coordinates = []
        self.load_structure(amino_list)
        self.score = self.calculate_score()

    def load_structure(self, amino_list):
        
        x = 0
        y = 0
        
        for i, amino in enumerate(amino_list):
            
            amino = Amino(i, amino)
            
            self.coordinates.append((x, y))
            amino.set_coordinate(x, y)
            
            while True:

                fold = self.get_fold()
                self.folds.append(fold)
                amino.set_fold = fold
                
                x_tmp, y_tmp = self.calculate_coordinate(fold, x, y)

                if (x_tmp, y_tmp) not in self.coordinates:
                    x = x_tmp
                    y = y_tmp
                    
                    break

    
    def get_fold(self):
        
        return random.choice([-2, -1, 1, 2])

    
    def calculate_coordinate(self, fold, x, y):
        
        if fold == 1 or fold == -1:
            x_tmp = x + fold
        else:
            x_tmp = x
            
        if fold == 2 or fold == -2:
            y_tmp = y + 0.5 * fold
        else:
            y_tmp = y
        
        return x_tmp, y_tmp

    def calculate_score(self):
        
        return 0