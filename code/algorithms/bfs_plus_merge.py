import copy
from code import Amino, Protein, ProteinTree

class BFSPlusMerge():

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
        default_param['max_queue_size_3'] = 10000
        default_param['max_savestack_size_2'] = 2000
        default_param['max_savestack_size_3'] = 4000

        ### HIER NOG EEN COMMENT
        self.pruning_depth = pruning_depth if pruning_depth is not None else default_param[f'pruning_depth_{dimension}']
        self.pruning_distance_factor = pruning_distance_factor if pruning_distance_factor is not None else default_param[f'pruning_distance_{dimension}']
        self.max_queue_size = max_queue_size if max_queue_size is not None else default_param[f'max_queue_size_{dimension}']
        self.max_savestack_size = max_savestack_size if max_savestack_size is not None else default_param[f'max_savestack_size_{dimension}']
        self.relevance_score = 0
        self.node_count = 0

        ### HIER NOG EEN COMMENT
        self.first_node = self.create_first_node()


    ### DEZE FOLD FUNCTIE MODULAIR MAKEN; MEERDERE KLEINERE FUNCTIE
    def fold(self, fold_position = 0):
        """
        ### DIT KORTER MAKEN, ALLEEN FUNCTIONELE DINGEN OVER DE CODE BEHANDELEN, GEEN THEORETISCHE

        Folds protein using breadth first search:
            1. lists possibilities, create node for each (with parent and depth)
            2. add nodes to list of nodes
            3. visit non visited nodes

        Pruning by relevance starts from the beginning (relevance is calculated internally based on best node and depth)
        Other prunings start at pruning depth.
            - Pruning distance factor is used to prune nodes at every depth of the tree: only nodes every depth * pruning_distance_factor are kept.

        Max queue size allows only a certain amount of nodes to be stored (and processed).
            - If max queue size is too small: nodes on the right of the tree will never be used (can't be stored in the queue).
            - If maz queue size is too large: more nodes will have to be processed and total processing time will be increased.
        
        Stores folded protein in attribute finished_folded_protein.
        """
        protein_size = len(self.source_protein.aminos)
        
        non_visited_nodes = []
        non_visited_nodes.append(self.first_node)
        visited_nodes = []
        good_but_pruned = []
        best_node = self.first_node
        
        # Goes through the queue of non_visited_nodes
        while len(non_visited_nodes) > 0:
            
            node = non_visited_nodes.pop(0)
            
            ### HIER NOG EEN COMMENT
            if node.depth < protein_size - 1:

                ### HIER NOG EEN COMMENT
                if node.score > self.relevance_score and node.depth > self.pruning_depth:
                    continue

                # Retrieves the current protein and its last added amino
                protein = node.current_protein
                current_amino = protein.aminos[node.depth]
                
                ### HIER NOG EEN COMMENT
                if current_amino is not None:
                    if node.depth == 0:
                        self.init_amino_coordinates(current_amino)

                    ### HIER NOG EEN COMMENT
                    folds = self.get_possible_folds(protein, current_amino.coordinate) if node.depth > 0 else [1]
                    next_amino = self.source_protein.aminos[node.depth + 1]
                    
                    # Goes through all possible folds from current protein
                    for fold in folds:
                        new_amino = Amino(node.depth + 1, next_amino.type)
                        current_amino.set_fold(fold)
                        new_amino.previous_amino = 0 - fold
                                        
                        # Computes new coordinate for the newly created amino after fold
                        if self.dimension == 3:
                            new_x, new_y, new_z = protein.calculate_coordinate(fold, current_amino.coordinate)
                            new_amino.set_coordinate((new_x, new_y, new_z))
                        else:
                            new_x, new_y = protein.calculate_coordinate(current_amino.fold, current_amino.coordinate)
                            new_amino.set_coordinate((new_x, new_y))
                        
                        # Copies current protein and add new amino at the end
                        new_protein = copy.deepcopy(protein)
                        new_protein.aminos.append(new_amino)
                        
                        curr_score = new_protein.calculate_score()
                        
                        # Relevance pruning
                        if curr_score <= self.relevance_score:
                
                            # Creates new node for the new protein
                            new_node = ProteinTree(new_protein, node, node.depth + 1, self.node_count + 1)
                            new_node.score = curr_score
                            self.node_count = self.node_count + 1
                            node.next_amino.append(new_node)

                            ### HIER EEN LOSSE FUNCTIE VAN MAKEN    
                            # Depth and distance pruning
                            if (len(non_visited_nodes) < self.max_queue_size) and (len(non_visited_nodes) == 0 or self.node_count >= (non_visited_nodes[-1].id + (node.depth * self.pruning_distance_factor))):
                                
                                # Adds node to queue
                                non_visited_nodes.append(new_node)

                                # If score has improved, update best node
                                if (
                                    new_node.score < best_node.score or
                                    (new_node.score <= best_node.score and new_node.depth == protein_size)
                                    ):
                                    best_node = new_node

                                    # Calculates new relevance score (this heuristic doesn't come from anywhere else and was invented by the student who wrote this code - it can probably be improved)
                                    if best_node.depth > 0 and node.depth >= self.pruning_depth:
                                        ### HIER EEN LOSSE FUNCTIE VAN MAKEN
                                        self.relevance_score = best_node.score + 3 - self.dimension + self.source_protein.source_string[0:node.depth].count('C') + (self.source_protein.source_string[0:node.depth].count('P') * (self.dimension * 2) / protein_size)

                            # The node has been pruned, but did have a good score: we store it somewhere safe.        
                            else:
                                # Add node to the queue of other good nodes that have been pruned.
                                if len(good_but_pruned) < self.max_savestack_size:
                                    good_but_pruned.append(new_node)

                                ### REMOVE THIS 
                                # Remove oldest good but pruned node to add a new one.
                                # elif len(good_but_pruned) == self.max_savestack_size:
                                #     good_but_pruned.pop(0)
                                #     good_but_pruned.append(new_node)

            # Maximal depth has been reached, let's stop here. 
            ### EXPLAIN WAT "STOP HERE" CONTAINS; e.g. 'update best score and protein'  
            elif node.depth == (protein_size - 1) and node.score == best_node.score:
                best_node = node
                protein = node.current_protein
                    
                self.finished_folded_protein = protein
                break

            # Updates queue and archive
            visited_nodes.append(node)

            # The queue of selected nodes is empty, but the end of the protein couldn't be reached
            # Takes a node from the stack of pruned nodes that had a decent score and start again from there
            if len(non_visited_nodes) == 0 and best_node.depth <= protein_size:
                if len(good_but_pruned) > 0:
                    retrieved_node = good_but_pruned.pop()
                    non_visited_nodes.append(retrieved_node)
                    best_node = retrieved_node
            
            # Reaches the maximal depth with best node: stops here.
            if (best_node.depth == protein_size - 1 or
                node.depth == protein_size - 1 and node.parent == best_node):
                break 
        
        # Picks best node
        protein = best_node.current_protein
                     
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
        first_protein = Protein(string=self.source_protein.aminos[0].type, dimensionality=self.dimension)

        if (self.dimension == 3):
            first_protein.aminos[0].set_coordinate((0, 0, 0))
        else:
            first_protein.aminos[0].set_coordinate((0, 0))    
        
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
        return all([amino_object.coordinate != coordinate for amino_object in protein.aminos])



    def set_score(self, protein=None):
        """
        Calculates protein score and sets score attribute.
        Takes attribute finished_folded_protein as default protein
        """
        
        if protein == None:
            protein = self.finished_folded_protein

        score = protein.calculate_score()
        protein.set_score(score)