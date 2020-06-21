
class Amino():

    def __init__(self, position, amino_type):
        self.id = position
        self.type = amino_type
        self.fold = None
        self.coordinate = None
        self.previous_amino = None

    def set_fold(self, fold):
        """
        Meaning of the fold:
        1   = positive step in the first dimension (X-axis direction)
        -1  = negative step in the first dimension (X-axis direction)
        2   = positive step in the second dimension (Y-axis direction)
        -2  = negative step in the second dimension (Y-axis direction)
        3   = positive step in the third dimension (Z-axis direction)
        -3  = negative step in the third dimension (Z-axis direction)
        """
        self.fold = fold

    def set_previous_amino(self, previous_amino):
        """
        = negative of fold of previous amino
        = fold direction which is occupied by the previous amino
        """
        self.previous_amino = previous_amino

    def set_coordinate(self, coordinate):
        self.coordinate = coordinate
    
    def reset_amino(self):
        """
        Sets the fold, coordinate and previous amino values to None
        """
        self.fold = None
        self.coordinate = None
        self.previous_amino = None
