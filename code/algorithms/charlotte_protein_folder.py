import copy
import logging
from datetime import datetime

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
        first_protein.aminos[0].set_coordinate((0, 0))
        self.first_node = ProteinTree(first_protein)
        self.pruning_depth = 8 #round(len(protein.get_aminos()) / 2)
        self.relevance_score = 0
        self.relevance_score_heur = "best_node.score + (self.source_protein.source_string.count('P') * 4/len(self.source_protein.source_string))"
        self.pruning_distance = 0
        self.pruning_distance_heur = "node.depth * 4"
        self.node_count = 0
        self.max_queue_size =2000

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
        now = datetime.now()
        current_time = now.strftime("%d%m_%H%M%S")
        #logging.basicConfig(filename="charlotte_" + current_time + ".log",level=#logging.deBUG)
        
        # set first coordinate to (0,0) and occupied fold to 0
        non_visited_nodes = []
        non_visited_nodes.append(self.first_node)
        visited_nodes = []
        good_but_pruned = []
        best_node = self.first_node
        
        # Goes through the queue of non_visited_nodes
        while len(non_visited_nodes) > 0:
            
            node = non_visited_nodes.pop(0)
            
            if node.depth >= len(self.source_protein.get_aminos()) - 1:
                #logging.debug(f'Node depth is { node.depth }, this is too high: skip.')
                continue

            if node.score > self.relevance_score and node.depth > self.pruning_depth:
                #logging.debug(f'Score is not interesting anymore (relevance score = { self.relevance_score }), skipping to next node.')
                continue

            assert isinstance(node, ProteinTree)
        
            # Retrieves the current protein and its last added amino
            protein = node.current_protein
            current_amino = protein.get_aminos()[node.depth]
            #logging.debug(f'Entering new non visited node with depth {node.depth}, current protein: {protein.to_string()}, current score: {node.score}.')
        
            if current_amino is not None:
                if node.depth == 0:
                    current_amino.set_coordinate((0, 0))
        
                x, y = current_amino.coordinate
    
                folds = self.get_possible_folds(protein, x, y) if node.depth > 0 else [1]
                next_amino = self.source_protein.get_aminos()[node.depth + 1]
                #logging.debug(f'Possible folds: {folds}')
                
                # Goes through all possible folds from current protein
                for fold in folds:
                    new_amino = Amino(node.depth + 1, next_amino.type)
                    current_amino.set_fold(fold)
                    new_amino.previous_amino = 0 - fold
                                    
                    # Computes new coordinate for the newly created amino after fold
                    new_x, new_y = calculate_coordinate(current_amino.fold, (x, y))
                    new_amino.set_coordinate((new_x, new_y))
                    #logging.debug(f'\t Trying fold {fold} with new amino {new_amino.type}. New coordinates: {new_x, new_y}')
                    
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
                            
                        if (len(non_visited_nodes) < self.max_queue_size) and (len(non_visited_nodes) == 0 or self.node_count >= (non_visited_nodes[-1].id + (node.depth * 4))):
                            #logging.debug(f'\t Score: {new_node.score}. Adding to non_visited_nodes, which now contains {len(non_visited_nodes) + 1} elements.')
                            
                            # Adds node to queue
                            non_visited_nodes.append(new_node)

                            # If score has improved, update best node
                            # !! To do: what if after pruning, the best node ends up leading to an impossible protein?
                            # -> Keep track of "best nodes"? Go back in tree?
                            if new_node.score < best_node.score or (new_node.score == best_node.score and new_node.depth == len(self.source_protein.aminos)):
                                best_node = new_node

                                if best_node.depth > 0 and node.depth >= self.pruning_depth:
                                    self.relevance_score = best_node.score + 1 + (self.source_protein.source_string.count('P') * 4/len(self.source_protein.source_string))
                                    #good_but_pruned.clear()
                                #logging.debug(f"Changing best node, new best node is at depth {best_node.depth}, new best score is {best_node.score}, new relevance score is {self.relevance_score}")
                                #logging.debug(f"Current protein = {new_protein.to_string_with_coord()}")
                        else:
                            if len(good_but_pruned) < self.max_queue_size:
                                good_but_pruned.append(new_node)
                            elif len(good_but_pruned) == self.max_queue_size:
                                good_but_pruned.pop(0)
                                good_but_pruned.append(new_node)
                            #logging.debug(f"Added node to good_but_pruned, current size = { len(good_but_pruned) }")



            # Updates queue and archive
            visited_nodes.append(node)

            if len(non_visited_nodes) == 0 and best_node.depth <= len(self.source_protein.aminos):
                non_visited_nodes.append(good_but_pruned.pop())
                #logging.debug("Taking node from good_but_pruned")
            #else: 
                #logging.debug(f"Non visited nodes size = { len(non_visited_nodes) }, best_node depth = { best_node.depth }, current depth = { node.depth }, origin size = { len(self.source_protein.aminos) }")
        
        # Picks best node
        protein = best_node.current_protein
                     
        self.finished_folded_protein = protein
        print(f"Depth = {best_node.depth}, original length = { len(self.source_protein.aminos) }, end length = { len(protein.aminos) }")


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