from .protein import Protein

class ProteinTree():

    def __init__(self, origin):
        self.protein = origin
        self.next_possibilities = []

    def add_possibility(self, protein_tree):
        # Checks that the argument is of the correct class, so as not to add a wrong object to the list of possibilities
        assert isinstance(protein_tree, ProteinTree)

        self.next_possibilities.append(protein_tree)

    def is_leaf(self, protein_tree):

        # Checks that the argument is of the correct class
        assert isinstance(protein_tree, ProteinTree)

    

        return (len(protein_tree.next_possibilities) == 0)