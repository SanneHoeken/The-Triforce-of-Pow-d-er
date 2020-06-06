
class Amino():

    def __init__(self, position, amino_type):
        self.id = position
        self.type = amino_type
        self.fold = None  # in 2D-grid -2, 2, -1 or 1
        self.coordinate = None
        self.previous_amino = None

    def set_previous_amino(self, previous_amino):
        self.previous_amino = previous_amino

    def set_fold(self, fold):
        self.fold = fold

    def set_coordinate(self, x, y):
        self.coordinate = (x, y)
    
    def reset_amino(self):
        self.fold = None
        self.coordinate = None
        self.previous_amino = None
