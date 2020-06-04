
class Amino():

    def __init__(self, position, amino_type):
        self.position = position
        self.type = amino_type
        self.fold = 0
        self.coordinate = (0, 0)
        self.occupied_fold = 0

    def set_occupied_fold(self, occupied_fold):
        self.occupied_fold = occupied_fold

    def set_fold(self, fold):
        self.fold = fold

    def set_coordinate(self, x, y):
        self.coordinate = (x, y)