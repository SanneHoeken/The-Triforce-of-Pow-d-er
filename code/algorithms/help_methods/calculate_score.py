from code.algorithms.help_methods.calculate_coordinate import calculate_coordinate

def calculate_score(protein):
    """
    Method that takes a protein and returns its stability score
    """ 
    score = 0

    # checks for every H- or C- amino in protein if not-connected neighbor is H-amino
    for amino in protein.get_aminos():
        
        if amino.coordinate is not None:

            if amino.type == 'H' or amino.type == 'C':

                # initializes not-connected directions if protein is 2D
                if len(amino.coordinate) == 2:
                    free_folds = list({-2, 2, -1, 1} - {amino.previous_amino, amino.fold})
                
                # initializes not-connected directions if protein is 3D
                elif len(amino.coordinate) == 3:
                    free_folds = list({-3, 3, -2, 2, -1, 1} - {amino.previous_amino, amino.fold})
                    
                # checks for every not-connected direction if H- or C- amino is present
                for fold in free_folds:
                    
                    # calculates coordinate of not-connected direction
                    next_coordinate = calculate_coordinate(fold, amino.coordinate)

                    # checks for hydrophobe/cysteine neighbor and changes score
                    neighbor = protein.get_amino(next_coordinate)
                    if neighbor is not None:
                        if neighbor.type == 'P':
                            score -= 0
                        elif amino.type == 'H' and neighbor.type == 'H':
                            score -= 1
                        elif amino.type == 'H' and neighbor.type == 'C':
                            score -= 1
                        elif amino.type == 'C' and neighbor.type == 'H':
                            score -= 1
                        elif amino.type == 'C' and neighbor.type == 'C':
                            score -= 5            
                
    return int(0.5 * score)