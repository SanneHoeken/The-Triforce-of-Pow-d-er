from .protein import Protein
from .amino import Amino

# Note: I haven't tried using this class yet!
class ProteinTree():

    def __init__(self, origin, depth = 0):
        self.amino = origin
        self.next_amino = []
        self.depth = depth

    def add_amino(self, amino):
        # Checks that the argument is of the correct class, so as not to add a wrong object to the list of possibilities
        if isinstance(amino, ProteinTree):
            self.next_amino.append(amino)

        if isinstance(amino, Amino):
            amino_tree = ProteinTree(amino, self.depth + 1)
            self.next_amino.append(amino_tree)


    def is_leaf(self, protein_tree):

        # Checks that the argument is of the correct class
        assert isinstance(protein_tree, ProteinTree)

        return (len(self.next_amino) == 0)


    def to_string(self):
        tree_repr = self.amino
        
        for amino in self.next_amino:
            tree_repr = tree_repr + "\n" + (self.depth * "\t") + "|_" + amino.to_string()

        return tree_repr