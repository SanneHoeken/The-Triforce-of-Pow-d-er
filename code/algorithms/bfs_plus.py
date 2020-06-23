import copy
from code import Amino, Protein, ProteinTree

class BFSplus():

    def __init__(self, protein, dimension = 2, pruning_depth = None, pruning_distance_factor = None, max_queue_size = None, max_savestack_size = None):
        """
        Takes a Protein object as argument.
        """
        self.source_protein = protein
        self.finished_folded_protein = None

        self.dimension = dimension
        
        # Sets default parameters, which give decent results in practice
        default_param = {}
        default_param['pruning_depth_2'] = 8
        default_param['pruning_depth_3'] = 4
        default_param['pruning_distance_2'] = 4
        default_param['pruning_distance_3'] = 6
        default_param['max_queue_size_2'] = 2000
        default_param['max_queue_size_3'] = 100000
        default_param['max_savestack_size_2'] = 2000
        default_param['max_savestack_size_3'] = 4000

        # Initializes attributes to passed arguments, or if None, to default parameters
        self.pruning_depth = pruning_depth if pruning_depth is not None else default_param[f'pruning_depth_{dimension}']
        self.pruning_distance_factor = pruning_distance_factor if pruning_distance_factor is not None else default_param[f'pruning_distance_{dimension}']
        self.max_queue_size = max_queue_size if max_queue_size is not None else default_param[f'max_queue_size_{dimension}']
        self.max_savestack_size = max_savestack_size if max_savestack_size is not None else default_param[f'max_savestack_size_{dimension}']
        self.relevance_score = 0
        self.node_count = 0

        # Initializes first node of the tree
        self.first_node = self.create_first_node()


    def fold(self):
        """
        Folds protein using breadth first search + 4 different types of prunings.
        Stores folded protein in attribute finished_folded_protein.
        """
        protein_total_size = len(self.source_protein.get_aminos())
        
        non_visited_nodes = []
        non_visited_nodes.append(self.first_node)
        visited_nodes = []
        good_but_pruned = []
        best_node = self.first_node
        
        # Goes through the queue of non_visited_nodes
        while len(non_visited_nodes) > 0:
            
            node = non_visited_nodes.pop(0)
            current_depth = node.get_depth()
            current_score = node.get_score()
            
            # Only continue if the current node isn't a leaf
            if current_depth < protein_total_size - 1:

                # If the score is not interesting and this node is past pruning depth: prunes
                if current_score > self.relevance_score and current_depth > self.pruning_depth:
                    continue

                # Retrieves the current protein and its last added amino
                protein = node.current_protein
                current_amino = protein.get_aminos()[current_depth]
                
                if current_amino is not None:
                    if current_depth == 0:
                        self.init_amino_coordinates(current_amino)

                    # If current node is first nodes, forces folding in one direction (to avoid the creation of symmetrical proteins)
                    # Otherwise, gets a list of possible folds from current protein
                    folds = self.get_possible_folds(protein, current_amino.get_coordinate()) if current_depth > 0 else [1]
                    next_amino = self.source_protein.get_aminos()[current_depth + 1]
                    
                    # Goes through all possible folds from current protein
                    for fold in folds:
                        new_amino = Amino(current_depth + 1, next_amino.get_type())
                        current_amino.set_fold(fold)
                        new_amino.set_previous_amino(0 - fold)
                                        
                        # Computes new coordinate for the newly created amino after fold
                        if self.dimension == 3:
                            new_x, new_y, new_z = protein.calculate_coordinate(fold, current_amino.get_coordinate())
                            new_amino.set_coordinate((new_x, new_y, new_z))
                        else:
                            new_x, new_y = protein.calculate_coordinate(fold, current_amino.get_coordinate())
                            new_amino.set_coordinate((new_x, new_y))
                        
                        # Copies current protein and add new amino at the end
                        new_protein = copy.deepcopy(protein)
                        new_protein.append_amino(new_amino)
                        
                        new_score = new_protein.calculate_score()
                        
                        # Relevance pruning
                        if new_score <= self.relevance_score:
                
                            # Creates new node for the new protein
                            new_node = ProteinTree(new_protein, node, current_depth + 1, self.node_count + 1)
                            new_node.set_score(new_score)
                            self.node_count = self.node_count + 1
                            node.add_amino(new_node)

                            # Check for depth and distance pruning: if no pruning, add node to the queue and check score.
                            if (len(non_visited_nodes) < self.max_queue_size) and (len(non_visited_nodes) == 0 or self.node_count >= (non_visited_nodes[-1].get_id() + (current_depth * self.pruning_distance_factor))):
                                
                                # Adds node to queue
                                non_visited_nodes.append(new_node)
                                best_score = best_node.get_score()

                                # If score has improved, update best node
                                if (
                                    new_score < best_score or
                                    (new_score <= best_score and new_node.get_depth() == protein_total_size)
                                    ):
                                    best_node = new_node

                                    # If pruning depth is past: relevance will be used to prune and thus needs to be updated.
                                    if best_node.get_depth() > 0 and current_depth >= self.pruning_depth:
                                        self.update_relevance(new_score, current_depth)

                            # The node has been pruned, but did have a good score: we store it somewhere safe.        
                            else:
                                # Add node to the queue of other good nodes that have been pruned.
                                if len(good_but_pruned) < self.max_savestack_size:
                                    good_but_pruned.append(new_node)

            # Maximal depth has been reached: settles on best node and protein, stops search.
            elif current_depth == (protein_total_size - 1) and current_score == best_node.get_score():
                best_node = node
                protein = node.current_protein
                    
                self.finished_folded_protein = protein
                break

            # Updates queue and archive
            visited_nodes.append(node)

            # The queue of selected nodes is empty, but the end of the protein couldn't be reached
            # Takes a node from the stack of pruned nodes that had a decent score and start again from there
            if len(non_visited_nodes) == 0 and best_node.get_depth() <= protein_total_size:
                if len(good_but_pruned) > 0:
                    retrieved_node = good_but_pruned.pop()
                    non_visited_nodes.append(retrieved_node)
                    best_node = retrieved_node
            
            # Reaches the maximal depth with best node: stops here.
            if (best_node.get_depth() == protein_total_size - 1 or
                node.get_depth() == protein_total_size - 1 and node.get_parent() == best_node):
                break 
        
        # Picks best node
        protein = best_node.get_current_protein()
                     
        self.finished_folded_protein = protein



    def get_possible_folds(self, protein, current_coordinates):  
        """
        Returns list of possible folds.
        """  
        possible_folds = []
        
        if self.dimension == 3:
            x, y, z = current_coordinates

            # Dict with all neighbouring coordinates 
            coordinates = {1: (x+1, y, z), -1: (x-1, y, z), 2:(x, y+1, z), -2: (x, y-1, z), 3: (x, y, z+1), -3: (x, y, z-1)} 

        else:
            x, y = current_coordinates

            # Dict with all neighbouring coordinates 
            coordinates = {1: (x+1, y), -1: (x-1, y), 2:(x, y+1), -2: (x, y-1)} 

        # Add possible folds to list
        for key, value in coordinates.items(): 
            if self.is_free_space(protein, value):
                possible_folds.append(key)

        return possible_folds



    def create_first_node(self):
        """
        Creates the origin ProteinTree node for BFSplus.
        Returns ProteinTree node.
        """
        first_protein = Protein(string=self.source_protein.get_aminos()[0].get_type(), dimensionality=self.dimension)

        if (self.dimension == 3):
            first_protein.get_aminos()[0].set_coordinate((0, 0, 0))
        else:
            first_protein.get_aminos()[0].set_coordinate((0, 0))    
        
        return ProteinTree(first_protein)

    

    def init_amino_coordinates(self, amino):
        """
        Initializes an amino to the 0 position, depending on the dimension in which the algorithm is started.
        """
        if (self.dimension == 3):
            amino.set_coordinate((0, 0, 0))
        else:
            amino.set_coordinate((0, 0))



    def is_free_space(self, protein, coordinate):
        """
        Returns True if coordinate is not occupied, else False.
        """
        return all([amino_object.get_coordinate() != coordinate for amino_object in protein.get_aminos()])



    def set_score(self, protein=None):
        """
        Calculates protein score and sets score attribute.
        Uses attribute finished_folded_protein as default protein.
        """
        
        if protein == None:
            protein = self.finished_folded_protein

        score = protein.calculate_score()
        protein.set_score(score)


    
    def update_relevance(self, score, depth):
        """
        Calculates a relevance score and stores it in self.relevance_score.
        This heuristic doesn't come from anywhere else and was invented by the student who wrote this code: it can probably be improved)
        Arguments: new score, current depth.
        Returns: nothing.
        """
        self.relevance_score = score + 3 - self.dimension + (self.source_protein.source_string[0:depth].count('P') * (self.dimension * 2) / len(self.source_protein.get_aminos()))