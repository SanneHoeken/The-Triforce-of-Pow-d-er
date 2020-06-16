
class Amino():

    def __init__(self, position, amino_type):
        self.id = position
        self.type = amino_type
        self.fold = None
        self.coordinate = None
        self.previous_amino = None

    def set_fold(self, fold):
        """
        1   = positieve stap in de eerste dimensie (X-as richting)
        -1  = negatieve stap in de eerste dimensie (X-as richting)
        2   = positieve stap in de tweede dimensie (Y-as richting)
        -2  = negatieve stap in de tweede dimensie (Y-as richting)
        3   = positieve stap in de derde dimensie (Z-as richting)
        -3  = negatieve stap in de derde dimensie (Z-as richting)
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
        self.fold = None
        self.coordinate = None
        self.previous_amino = None
