def get_symmetry(self, folds):
    """
    Takes a list of folds and returns all lists of folds that are symmteric to input list
    """
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

    return s2, s3, s4, s1_x, s1_y, s2_x, s2_y