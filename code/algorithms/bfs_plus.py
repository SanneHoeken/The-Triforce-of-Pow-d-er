import copy
import logging
from datetime import datetime

from code.classes.amino import Amino
from code.classes.protein import Protein
from code.algorithms.help_methods.calculate_score import calculate_score
from code.algorithms.help_methods.calculate_coordinate import calculate_coordinate
from code.classes.protein_tree import ProteinTree

class BFSPlus():

    def __init__(self, protein, pruning_depth = 8, pruning_distance_factor = 4, max_queue_size = 2000):
        """
        Takes a Protein object as argument.
        """
        self.source_protein = protein
        self.finished_folded_protein = None
        first_protein = Protein(string=protein.get_aminos()[0].type)
        first_protein.aminos[0].set_coordinate((0, 0)) 
        self.first_node = ProteinTree(first_protein)
        self.pruning_depth = pruning_depth # Advice: 8
        self.relevance_score = 0
        self.node_count = 0
        self.pruning_distance_factor = pruning_distance_factor # Advice: 4
        self.max_queue_size = max_queue_size # Advice: 2000



    def fold(self, fold_position = 0):
        """
        Stores folded protein in attribute finished_folded_protein.
        Folds protein using breadth first search:
            1. lists possibilities, create node for each (with parent and depth)
            2. add nodes to list of nodes
            3. visit non visited nodes

        Pruning by relevance starts from the beginning (relevance is calculated internally based on best node and depth)
        Other prunings start at pruning depth.
            - Pruning distance factor is used to prune nodes at every depth of the tree: only nodes every depth * pruning_distance_factor are kept.

        Max queue size allows only a certain a mount of nodes to be stored (and processed).
            - If max queue size is too small: nodes on the right of the tree will never be used (can't be stored in the queue).
            - If maz queue size is too large: more nodes will have to be processed and total processing time will be increased.
        Returns nothing.
        """
        protein_size = len(self.source_protein.get_aminos())
        non_visited_nodes = []
        non_visited_nodes.append(self.first_node)
        visited_nodes = []
        good_but_pruned = []
        best_node = self.first_node
        
        # Goes through the queue of non_visited_nodes
        while len(non_visited_nodes) > 0:
            
            node = non_visited_nodes.pop(0)
            
            # Depth is too high, we shouldn't go further, skip
            if node.depth >= protein_size - 1:
                continue

            # Score is too high and irrelevant, skip
            if node.score > self.relevance_score and node.depth > self.pruning_depth:
                continue

            # Retrieves the current protein and its last added amino
            protein = node.current_protein
            current_amino = protein.get_aminos()[node.depth]
            
            if current_amino is not None:
                if node.depth == 0:
                    current_amino.set_coordinate((0, 0))
        
                x, y = current_amino.coordinate
    
                folds = self.get_possible_folds(protein, x, y) if node.depth > 0 else [1]
                next_amino = self.source_protein.aminos[node.depth + 1]
                
                # Goes through all possible folds from current protein
                for fold in folds:
                    new_amino = Amino(node.depth + 1, next_amino.type)
                    current_amino.set_fold(fold)
                    new_amino.previous_amino = 0 - fold
                                    
                    # Computes new coordinate for the newly created amino after fold
                    new_x, new_y = calculate_coordinate(current_amino.fold, (x, y))
                    new_amino.set_coordinate((new_x, new_y))
                    
                    # Copies current protein and add new amino at the end
                    new_protein = copy.deepcopy(protein)
                    new_protein.aminos.append(new_amino)
                    
                    curr_score = calculate_score(new_protein)
                    
                    # Pruning: only adds 1 node per pruning_distance
                    if curr_score <= self.relevance_score:
                        
                        # Creates new node for the current protein
                        new_node = ProteinTree(new_protein, node, node.depth + 1, self.node_count + 1)
                        new_node.score = curr_score
                        self.node_count = self.node_count + 1
                        node.next_amino.append(new_node)
                            
                        if (len(non_visited_nodes) < self.max_queue_size) and (len(non_visited_nodes) == 0 or self.node_count >= (non_visited_nodes[-1].id + (node.depth * self.pruning_distance_factor))):
                            
                            # Adds node to queue
                            non_visited_nodes.append(new_node)

                            # If score has improved, update best node
                            if new_node.score < best_node.score or (new_node.score == best_node.score and new_node.depth == protein_size):
                                best_node = new_node

                                # Calculates new relevance score (this heuristic doesn't come from anywhere else and was invented by the student who wrote this code - it can probably be improved)
                                if best_node.depth > 0 and node.depth >= self.pruning_depth:
                                    self.relevance_score = best_node.score + 1 + (self.source_protein.source_string.count('P') * 4/protein_size)
                                    
                        else:

                            # Add node to the queue of other good nodes that have been pruned.
                            if len(good_but_pruned) < self.max_queue_size:
                                good_but_pruned.append(new_node)
                            
                            # Remove oldest good but pruned node to add a new one.
                            elif len(good_but_pruned) == self.max_queue_size:
                                good_but_pruned.pop(0)
                                good_but_pruned.append(new_node)
                            
            # Updates queue and archive
            visited_nodes.append(node)

            if len(non_visited_nodes) == 0 and best_node.depth <= protein_size:
                non_visited_nodes.append(good_but_pruned.pop())
                        
        # Picks best node
        protein = best_node.current_protein
                     
        self.finished_folded_protein = protein
        


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