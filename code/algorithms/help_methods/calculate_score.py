from code.algorithms.help_methods.calculate_coordinate import calculate_coordinate

def calculate_score(protein):
        
    score = 0

    # checks for every H- or C- amino in protein if not-connected neighbor is H-amino
    for amino in protein.get_aminos():
        
        if amino.coordinate is not None:

            if amino.type == 'H' or amino.type == 'C':
                x, y = amino.coordinate

                # initializes not-connected directions
                free_folds = list({-2, 2, -1, 1} - {amino.previous_amino, amino.fold})
                
                # checks for every not-connected direction if H-amino is present
                for fold in free_folds:
                    next_coordinate = calculate_coordinate(fold, x, y)

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