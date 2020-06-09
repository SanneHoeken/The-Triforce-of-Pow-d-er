import queue
import copy
from code.classes.amino import Amino
from code.classes.protein import Protein
from code.classes.algorithms.random_protein_folder import RandomProteinFolder
from code.algorithms.calculate_coordinate import calculate_coordinate
from code.visualisations import visualize

class Experiment():
    """
    Experimenting with symmetry
    """
    
    def get_symmetry(self, folds):
        rotate_90 = {-2: 1, 2: -1, -1: -2, 1: 2}
        y_mirror = {-2: 2, 2: -2, -1: -1, 1: 1}
        x_mirror = {-2: -2, 2: 2, -1: 1, 1: -1}
        
        s1 = folds[:]
        s2 = [rotate_90[fold] for fold in s1]
        s3 = [rotate_90[fold] for fold in s2]
        s4 = [rotate_90[fold] for fold in s3]
        s1_x = [x_mirror[fold] for fold in s1]
        s1_y = [y_mirror[fold] for fold in s1]
        s2_x = [x_mirror[fold] for fold in s2]
        s2_y = [y_mirror[fold] for fold in s2]

        return s1, s2, s3, s4, s1_x, s1_y, s2_x, s2_y


    def fold_protein(self, protein, folds):
        """Takes protein and list of folds
        returns folded part of the protein"""

        folded_protein = Protein(protein=protein)
        folded_protein.aminos = protein.get_aminos()[:len(folds) + 1]

        x = 0
        y = 0
        previous_amino = 0

        for amino in folded_protein.get_aminos():
            
            # set amino's coordinates
            amino.set_coordinate(x, y)

            # set amino's occupied fold
            amino.set_previous_amino(previous_amino)
            
            if amino.id < len(folds):
                # set fold and move fold position 
                fold = folds[amino.id]
                amino.set_fold(fold)
                    
                # compute next coordinate following the fold
                new_x, new_y = calculate_coordinate(fold, x, y)

                # set next coordinate values
                x = new_x
                y = new_y

                # set previous amino to inverse fold
                previous_amino = -fold
        
        return folded_protein

    
    def breadth_first(self, protein):
        # implementation of breadth first
        depth = len(protein.get_aminos()) - 1
        q = queue.Queue()
        q.put([])
        
        while not q.empty():
            state = q.get()
            if len(state) < depth:
                for i in [-1, 1, -2, 2]:
                    child = copy.deepcopy(state)
                    child.append(i)
                    q.put(child)


    def depth_first(self, protein):
        # implementation of depth first
        depth = len(protein.get_aminos()) - 1
        stack = [[]]
        
        while len(stack) > 0:
            state = stack.pop()
            if len(state) < depth:
                for i in [-1, 1, -2, 2]:
                    child = copy.deepcopy(state)
                    child.append(i)
                    stack.append(child)


# TEST CODE
from code.visualisations import visualize

if __name__ == "__main__":
    string = "HHPHPH"

    protein = Protein(string=string)
    folds = [2, 1, 2, 1, -2]
    
    # get symmetric configurations
    symmetric_configs = Experiment().get_symmetry(folds)

    # plot all symmetric configurations
    for config in symmetric_configs:
        folded = Experiment().fold_protein(protein, config)
        visualize.visualize_protein(folded)