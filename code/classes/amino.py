
class Amino():

    def __init__(self, position, amino_type):
        self.id = position
        self.type = amino_type
        self.fold = None
        self.coordinate = None
        self.occupied_fold = None

    def set_occupied_fold(self, occupied_fold):
        self.occupied_fold = occupied_fold

    def set_fold(self, fold):
        self.fold = fold

    def set_coordinate(self, x, y):
        self.coordinate = (x, y)

    def __repr__(self):
        return f"Amino({self.id})"