import random
from code.classes.amino import Amino
from code.classes.protein import Protein
from code.classes.protein_folder import ProteinFolder

# The goal of this class will be to take a protein and try different random folding possibilities until it is satisfied with the folding
class RandomProteinFolder(ProteinFolder):

    def fold(self, protein_to_fold = None, fold_position = 0):

        protein = protein_to_fold

        # If no argument is passed: uses origin protein
        if protein == None:
            assert self.origin_protein
            protein = self.origin_protein

        # set first coordinate to (0,0) and occupied fold to 0
        x = 0
        y = 0
        previous_amino = 0

        # repeat folding until all aminos are folded
        while fold_position < (len(protein.get_aminos()) - 2):    

            # iterate over every amino
            for amino in protein.get_aminos()[fold_position:]:

                # set amino's coordinates
                amino.set_coordinate(x, y)

                # set amino's occupied fold
                amino.set_previous_amino(previous_amino)
                
                # generate fold
                fold = self.get_fold(protein, x, y)
                
                # break if protein ran into dead end
                if fold == 0:
                    
                    # reset all amino values
                    for amino in protein.get_aminos()[:fold_position+1]:
                        amino.reset_amino()
                    
                    # reset fold position
                    fold_position = 0
                    break  
                
                # set fold and move fold position 
                amino.set_fold(fold)
                fold_position += 1
                    
                # compute next coordinate following the fold
                new_x, new_y = self.calculate_coordinate(fold, x, y)

                # set next coordinate values
                x = new_x
                y = new_y

                # set previous amino to inverse fold
                previous_amino = -fold
            
        self.finished_folded_protein = protein