import random
import copy
from code import Amino, Protein, calculate_coordinate, calculate_score


def fold(protein, iterations=1):

    best_score = 1

    # fold protein iterations times 
    for i in range(iterations):
        
        # make copy of protein, fold the copy, calculate score
        protein_copy = copy.deepcopy(protein)
        fold_candidate(protein_copy)
        score = calculate_score(protein_copy)

        # keep protein with lowest score
        if score < best_score:
            best_score = score
            best_protein = protein_copy
            iteration = i
    
    print(f'Best result at iteration: {iteration} (out of {iterations})')
    protein = best_protein
    protein.set_score(best_score)
    
    return protein
    

def fold_candidate(protein_to_fold, fold_position = 0):

    protein_candidate = protein_to_fold

    # set first coordinate to (0,0) and occupied fold to 0
    x = 0
    y = 0
    previous_amino = 0

    # repeat folding until all aminos are folded
    while fold_position < (len(protein_candidate.get_aminos()) - 2):
        
        # iterate over every amino
        for amino in protein_candidate.get_aminos()[fold_position:]:

            # set amino's coordinates
            amino.set_coordinate(x, y)

            # set amino's occupied fold
            amino.set_previous_amino(previous_amino)
            
            # generate fold
            fold = get_fold(protein_candidate, x, y)
            
            # break if protein ran into dead end
            if fold == 0:
                
                # reset all amino values
                for amino in protein_candidate.get_aminos()[:fold_position+1]:
                    amino.reset_amino()
                
                # reset fold position
                fold_position = 0
                break  
            
            # set fold and move fold position 
            amino.set_fold(fold)
            fold_position += 1
                
            # compute next coordinate following the fold
            new_x, new_y = calculate_coordinate(fold, x, y)

            # set next coordinate values
            x = new_x
            y = new_y

            # set previous amino to inverse fold
            previous_amino = -fold


def get_fold(protein, x, y):
    # get possible folds
    possible_folds = get_possible_folds(protein, x, y)

    # choose random fold
    if len(possible_folds) == 0: 
        return 0
            
    fold = random.choice(possible_folds)
    
    return fold


def get_possible_folds(protein, x, y):    
    possible_folds = []
    
    # dict with all neighbouring coordinates 
    coordinates = {1: (x+1, y), -1: (x-1, y), 2:(x, y+1), -2: (x, y-1)} 

    # add possible folds to list
    for key, value in coordinates.items(): 
        if is_free_space(protein, value):
            possible_folds.append(key)

    return possible_folds


def is_free_space(protein, coordinate):

    # return True if coordinate is not occupied, else False
    return all([amino_object.coordinate != coordinate for amino_object in protein.get_aminos()])