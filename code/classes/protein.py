import random
from .amino import Amino

class Protein():

    def __init__(self, source_file):
        self.aminos = self.load_protein(source_file)
        self.fold_protein()
        self.score = self.calculate_score()

    def load_protein(self, source_file):

        protein = []

        # open source file and read first line
        with open(source_file, 'r') as infile:
            protein_string = infile.readline()

        # transform string into list of chars
        amino_list = [element for element in protein_string]

        # iterate over every amino
        for i, amino_type in enumerate(amino_list):

            # initialize amino
            amino = Amino(i, amino_type)

            # add amino to protein-structure
            protein.append(amino)

        return protein

    def fold_protein(self):
        
        # set first coordinate to (0,0) and occupied fold to 0
        x = 0
        y = 0
        occupied_fold = 0
        
        # iterate over every amino
        for amino in self.aminos:

            # set amino's coordinates
            amino.set_coordinate(x, y)

            # set amino's occupied fold
            amino.set_occupied_fold(occupied_fold)
            
            # fold amino towards a free coordinate
            # Charlotte note: during prog1 I had been told what while True wasn't really good practice.
            # -> Why not "while not all" + the conditionline 62? (Or better, a do{} while but not sure that exists in Python)
            while True:
                
                # generate and set fold
                fold = self.get_fold()
                amino.set_fold(fold)
                
                # compute next coordinate following the fold
                x_tmp, y_tmp = self.calculate_coordinate(fold, x, y)

                # keep coordinate and fold if it is not occupied in protein
                if all([amino_object.coordinate != (x_tmp, y_tmp) for amino_object in self.aminos]):
                    break
            
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
        # Note Charlotte: mag ik het proberen te implementeren? Ik begrijp je code helemaal maar om het helemaal vast in mijn hoofd
        # te krijgen zou ik het fijn vinden om het even proberen te modificeren en gebruiken!
        return 'H'