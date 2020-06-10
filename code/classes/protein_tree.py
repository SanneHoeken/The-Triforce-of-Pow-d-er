from .protein import Protein
from .amino import Amino

# Note: I haven't tried using this class yet!
class ProteinTree():

    def __init__(self, origin, parent = None, depth = 0):
        """
        ProteinTree class
        Init takes origin (Protein object) as argument.
        Optionnal arguments: parent (ProteinTree object, default None), depth (integer, default 0)
        Other attributes: 
            - next_amino: list of ProteinTree objects, default is []
            - score: default is 0
        """
        self.current_protein = origin
        self.parent = None
        self.next_amino = []
        self.depth = depth
        self.score = 0

    def add_amino(self, amino):
        """
        Adds amino (ProteinTree object) as node to the list attribute next_amino.
        """
        # Checks that the argument is of the correct class, so as not to add a wrong object to the list of possibilities
        if isinstance(amino, ProteinTree):
            self.next_amino.append(amino)

        if isinstance(amino, Amino):
            amino_tree = ProteinTree(amino, self.depth + 1)
            self.next_amino.append(amino_tree)


    def is_leaf(self, protein_tree):
        """
        Returns True if the argument is a lead, False otherwise.
        Argument: ProteinTree object.
        """
        # Checks that the argument is of the correct class
        assert isinstance(protein_tree, ProteinTree)

        return (len(self.next_amino) == 0)


    def to_string(self):
        """
        Returns a (to be improved) string representation of the tree.
        Arguments: none.
        """
        current_amino = self.current_protein.get_aminos().pop()
        fold = current_amino.fold if current_amino.fold is not None else 0
        tree_repr = current_amino.type + " (" + str(fold) + "), score = " + str(self.score)
        
        for amino in self.next_amino:
            tree_repr = tree_repr + "\n" + (self.depth * " ") + "|_" + amino.to_string()

        return tree_repr