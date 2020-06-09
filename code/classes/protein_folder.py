import random
from .amino import Amino
from .protein import Protein
from .protein_tree import ProteinTree

# The goal of this class will be to take a protein and try different random folding possibilities until it is satisfied with the folding
class ProteinFolder():

    def __init__(self, protein = None):
        self.origin_protein = protein
        self.finished_folded_protein = None
        self.proteins_tree = None

    
    def set_origin_protein(self, protein):
        self.origin_protein = protein

    
    def fold(self, protein_to_fold = None, fold_position = 0):
        # This method does nothing, please override it in children classes!
        pass