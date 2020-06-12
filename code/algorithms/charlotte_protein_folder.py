import copy
import logging

from code.classes.amino import Amino
from code.classes.protein import Protein
from code.algorithms.help_methods.calculate_score import calculate_score
from code.algorithms.help_methods.calculate_coordinate import calculate_coordinate
from code.classes.protein_tree import ProteinTree

class CharlotteProteinFolder():

    def __init__(self, protein):
        """
        Takes a Protein object as argument.
        """
        self.source_protein = protein
        self.finished_folded_protein = None
        first_protein = Protein(string=protein.get_aminos()[0].type)
        first_protein.aminos[0].set_coordinate(0, 0)
        self.first_node = ProteinTree(first_protein)
        self.pruning_depth = round(len(protein.get_aminos()))
    

    def fold(self, fold_position = 0):
        """
        Stores folded protein in attribute finished_folded_protein.
        Folds protein using breadth first search:
            1. lists possibilities, create node for each (with parent and depth)
            2. add nodes to list of nodes
            3. visit non visited nodes

        No pruning yet.
        Returns nothing.
        """
        # logging.basicConfig(filename='charlotte.log',level=logging.DEBUG)
        
        # set first coordinate to (0,0) and occupied fold to 0
        non_visited_nodes = []
        non_visited_nodes.append(self.first_node)
        visited_nodes = []
        best_node = self.first_node
        
        # Goes through the queue of non_visited_nodes
        while len(non_visited_nodes) > 0:
            
            node = non_visited_nodes.pop()
            
            if node.depth >= len(self.source_protein.get_aminos()) - 1:
                continue

            assert isinstance(node, ProteinTree)
        
            # Retrieves the current protein and its last added amino
            protein = node.current_protein
            current_amino = protein.get_aminos()[node.depth]
            # logging.debug(f'Entering new non visited node with depth {node.depth}, current protein: {protein.to_string()}, current score: {node.score}.')
        
            if current_amino is not None:
                if node.depth == 0:
                    current_amino.set_coordinate(0, 0)
        
                x, y = current_amino.coordinate
    
                folds = self.get_possible_folds(protein, x, y) if node.depth > 0 else [1]
                # logging.debug(f'Possible folds: {folds}')
                
                # Goes through all possible folds from current protein
                for fold in folds:
                    new_amino = Amino(node.depth + 1, self.source_protein.get_aminos()[node.depth + 1].type)
                    new_amino.set_fold(fold)
                    new_amino.previous_amino = -fold
                                    
                    # Computes new coordinate for the newly created amino after fold
                    new_x, new_y = calculate_coordinate(new_amino.fold, x, y)
                    new_amino.set_coordinate(new_x, new_y)
                    # logging.debug(f'\t Trying fold {fold} with new amino {new_amino.type}. New coordinates: {new_x, new_y}')
                    
                    # Copies current protein and add new amino at the end
                    new_protein = copy.deepcopy(protein)
                    new_protein.aminos.append(new_amino)
                    
                    # Creates new node for the current protein
                    new_node = ProteinTree(new_protein, node, node.depth + 1)
                    new_node.score = calculate_score(new_protein)
                    node.next_amino.append(new_node)
                    # logging.debug(f'\t Score: {new_node.score}. Adding to non_visited_nodes, which now contains {len(non_visited_nodes) + 1} elements.')
                     
                    # Adds node to queue
                    non_visited_nodes.append(new_node)

                    # If score has improved, update best node
                    # !! To do: what if after pruning, the best node ends up leading to an impossible protein?
                    # -> Keep track of "best nodes"? Go back in tree?
                    if new_node.score <= best_node.score and new_node.depth >= best_node.depth:
                        # logging.debug(f"Changing best node, new best node is at depth {best_node.depth}, new best score is {best_node.score}")
                        best_node = new_node

                # Updates queue and archive
                visited_nodes.append(node)
        
        # Picks best node
        protein = best_node.current_protein
                     
        self.finished_folded_protein = protein


    def get_random_fold(self, protein, x, y):
        """
        Returns random possible fold.
        Arguments: protein, x, y
        """
        # get possible folds
        possible_folds = self.get_possible_folds(protein, x, y)

        # choose random fold
        if len(possible_folds) == 0: 
            return 0
                
        fold = random.choice(possible_folds)
        
        return fold


    def get_possible_folds(self, protein, x, y):  
        """
        Returns list of possible folds.
        Arguments: protein, x, y.
        """  
        possible_folds = []
        
        # dict with all neighbouring coordinates 
        coordinates = {1: (x+1, y), -1: (x-1, y), 2:(x, y+1), -2: (x, y-1)} 

        # add possible folds to list
        for key, value in coordinates.items(): 
            if self.is_free_space(protein, value):
                possible_folds.append(key)

        return possible_folds


    def is_free_space(self, protein, coordinate):
        """
        Returns True if coordinate is not occupied, else False.
        """
        return all([amino_object.coordinate != coordinate for amino_object in protein.get_aminos()])


    def set_score(self, protein=None):
        """
        Calculates protein score and sets score attribute.
        Argument: protein. (Default: attribute finished_folded_protein)
        """
        
        if protein == None:
            protein = self.finished_folded_protein

        score = calculate_score(protein)
        protein.set_score(score)