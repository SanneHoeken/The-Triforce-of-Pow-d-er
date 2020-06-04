import random
from amino import Amino

class Protein():

    def __init__(self, amino_list):
        self.aminos = []
        self.length = len(amino_list)
        self.load_structure(amino_list)
        self.score = self.calculate_score()

    def load_structure(self, amino_list):
        
        # set first coordinate to (0,0) and occupied fold to 0
        x = 0
        y = 0
        occupied_fold = 0
        
        # iterate over every amino
        for i, amino_type in enumerate(amino_list):

            # initialize amino
            amino = Amino(i, amino_type)
            
            # set amino's coordinates
            amino.set_coordinate(x, y)

            # set amino's occupied fold
            amino.set_occupied_fold(occupied_fold)
            
            # fold amino towards a free coordinate
            while True:
                
                # generate and set fold
                fold = self.get_fold()
                amino.set_fold(fold)
                
                # compute next coordinate following the fold
                x_tmp, y_tmp = self.calculate_coordinate(fold, x, y)

                # keep coordinate and fold if it is not occupied in protein
                if all([amino_object.coordinate != (x_tmp, y_tmp) for amino_object in self.aminos]):
                    break

            # add amino to protein-structure
            self.aminos.append(amino)
            
            # set next coordinate values
            x = x_tmp
            y = y_tmp

            # set next occupied fold to inverse fold
            occupied_fold = -fold
    

    def get_fold(self):
        
        # generate random fold
        return random.choice([-2, -1, 1, 2])


    def calculate_coordinate(self, fold, x, y):
     
        # compute x-value following the fold
        if fold == 1 or fold == -1:
            x_tmp = x + fold

        # keep current x-value if no move along x-axis 
        else:
            x_tmp = x

        # compute y-value following the fold    
        if fold == 2 or fold == -2:
            y_tmp = y + 0.5 * fold
        
        # keep current y-value if no move along y-axis
        else:
            y_tmp = y
        
        return int(x_tmp), int(y_tmp)


    def calculate_score(self):
        
        score = 0

        # checks for every H-amino in protein if not-connected neighbor is H-amino
        for amino in self.aminos:
            
            if amino.type == 'H':
                x, y = amino.coordinate

                # initializes not-connected directions
                free_folds = list({-2, 2, -1, 1} - {amino.occupied_fold, amino.fold})
                
                # checks for every not-connected direction if H-amino is present
                for fold in free_folds:
                    next_coordinate = self.calculate_coordinate(fold, x, y)

                    #TO CONTINUE
                    neighbor = self.get_amino(next_coordinate)
                    if neighbor == 'H':
                        score -= 1

        return score

    def get_amino(self, coordinate):
        
        # returns amino that is present on a certain coordinate

        #TO DO
        return 'H'