import copy
import logging
import random
from datetime import datetime
import cProfile, pstats, io
import profile

from code.classes.amino import Amino
from code.classes.protein import Protein
from code.classes.protein_tree import ProteinTree

class BBProteinFolder():

    def __init__(self, protein):
        """
        Takes a Protein object as argument.
        """
        self.source_protein = protein
        self.finished_folded_protein = None
        first_protein = Protein(string=protein.get_aminos()[0].type)
        first_protein.aminos[0].set_coordinate((0, 0))
        self.first_node = ProteinTree(first_protein)
        self.total_depth = round(len(protein.get_aminos()))
        self.best_score = 0
        self.dict_best = dict()
        self.dict_count = dict()  
        self.dict_avg = dict()

    def fold(self, fold_position = 0):
        """
        Stores folded protein in attribute finished_folded_protein.
        Folds protein using depth first search:
            1. lists possibilities, create node for each (with parent and depth)
            2. add nodes to list of nodes
            3. visit non visited nodes
        """
        #now = datetime.now()
        #current_time = now.strftime("%d%m_%H%M%S")
        #logging.basicConfig(filename="thomas_" + current_time + ".log",level=logging.DEBUG)
        
        # set first coordinate to (0,0) and occupied fold to 0
        non_visited_nodes = []
        non_visited_nodes.append(self.first_node)
        best_node = self.first_node
        
        # initialise dictionaries
        for amino_object in self.source_protein.get_aminos():
            self.dict_avg[amino_object.id] = 0
            self.dict_best[amino_object.id] = 0
            self.dict_count[amino_object.id] = 0

        # Goes through the queue of non_visited_nodes TODO
        while len(non_visited_nodes) > 0:

            node = non_visited_nodes.pop()

            assert isinstance(node, ProteinTree)
        
            #logging.debug(f'best_node score: {best_node.score}, depth: {node.depth}.')

            # Retrieves the current protein and its last added amino
            protein = node.current_protein
            current_amino = protein.get_aminos()[node.depth]
            
            #logging.debug(f'Entering new non visited node with depth {node.depth}, current protein: {protein.to_string()}, current score: {node.score}.')

            if current_amino is not None: 
                
                if node.depth == 0:
                    current_amino.set_coordinate((0, 0))
                
                x, y = current_amino.coordinate
    
                folds = self.get_possible_folds(protein, x, y) if node.depth > 0 else [1]
                #logging.debug(f'Possible folds: {folds}')

             
                if node.depth == len(self.source_protein.get_aminos()) - 1:
                    continue

                # Goes through all possible folds from current protein
                for fold in folds:
                    
                    new_amino = Amino(node.depth + 1, self.source_protein.get_aminos()[node.depth + 1].type)
                    current_amino.set_fold(fold)
                    new_amino.previous_amino = 0 - fold

                    # Computes new coordinate for the newly created amino after fold
                    new_x, new_y = protein.calculate_coordinate(current_amino.fold, (x, y))
                    new_amino.set_coordinate((new_x, new_y))
                    #logging.debug(f'\t Trying fold {fold} with new amino {new_amino.type}. New coordinates: {new_x, new_y}')
                    
                    # Copies current protein and add new amino at the end
                    new_protein = copy.deepcopy(protein)
                    new_protein.aminos.append(new_amino)
                    
                    curr_score = new_protein.calculate_score()
                    
                    # data for candidate node
                    depth_best_score = self.dict_best[node.depth +1]
                    depth_avg_score = self.dict_avg[node.depth + 1]
                    depth_count = self.dict_count[node.depth + 1]

                    if self.is_viable(new_amino.type, depth_best_score, depth_avg_score, depth_count, curr_score) == 0:

                        # Creates new node for the current protein
                        new_node = ProteinTree(new_protein, node, node.depth + 1)
                        new_node.score = curr_score
                        current_protein = new_node.current_protein

                        node.next_amino.append(new_node)
                        #logging.debug(f'\t Score: {new_node.score}. Adding to non_visited_nodes, which now contains {len(non_visited_nodes) + 1} elements.')
                                
                        new_count = self.dict_count[new_node.depth] + 1
                        new_avg = round((self.dict_avg[new_node.depth] * self.dict_count[new_node.depth]) + current_protein.calculate_score()) / new_count
                        self.dict_avg[new_node.depth] = new_avg
                        self.dict_count[new_node.depth] = new_count

                        if current_protein.calculate_score() <= self.dict_best[new_node.depth]:
                            self.dict_best[new_node.depth] = current_protein.calculate_score()
                            #print(f'depth: {new_node.depth}, score: {self.total_depth}')
                        if new_node.depth == self.total_depth - 1:
                            best_node = new_node
                            best_node.score = new_node.score
                            #print('best_node updated')

                        # Adds node to queue
                        non_visited_nodes.append(new_node)

        # Picks best node
        #print(best_node.score)
        protein = best_node.current_protein
                     
        self.finished_folded_protein = protein
        print(f"Protein score: {protein.calculate_score()}")
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

        score = protein.calculate_score()
        protein.set_score(score)
    
    def is_viable(self, aminotype, depth_best_score, depth_avg_score, depth_count, curr_score):
        if aminotype == 'P':
            return 0 
        if curr_score <= depth_best_score: 
            return 0
        if curr_score <= depth_avg_score and self.cointoss(2):
            return 0
        if curr_score > depth_avg_score and self.cointoss(1):
            return 0
        return 1

        # returns T/F depending on given threshold and chance
    def cointoss(self, threshold):
        random_value = random.randint(0,101)
        if random_value < threshold:
            return 1
        return 0
