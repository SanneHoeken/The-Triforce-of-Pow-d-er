import random
from code import Amino, Protein

class BestOfRandom():
    """
    The BestOfRandom Class that randomly folds a protein i times and saves the best configuration
    """
    def __init__(self, protein, iterations=1):
        self.protein = protein
        self.iterations = iterations
        self.best_score = 1
        self.best_protein = None
        self.finished_folded_protein = self.protein    

    def run(self):
        """
        Folds protein randomly for every iteration, stores its values if configuration has a lower score,
        and then resets protein. After the iterations, the best values are set in the protein.
        """
        # Fold protein iterations times
        for i in range(self.iterations):
            random_folder = RandomFolder(self.protein)
            random_folder.fold_protein()

            # Calculate score
            score = random_folder.protein.get_score()

            # Update protein values en score if score is lower
            if score < self.best_score:
                self.best_score = score
                self.best_protein = [(amino.get_fold(), amino.get_coordinate(), amino.get_previous_amino()) for amino in random_folder.protein.get_aminos()]

            # Reset protein
            for amino in self.protein.get_aminos()[1:]:
                amino.reset_amino()

        # Set values of protein to best values
        for i, values in enumerate(self.best_protein):
            self.protein.get_aminos()[i].set_fold(values[0])
            self.protein.get_aminos()[i].set_coordinate(values[1])
            self.protein.get_aminos()[i].set_previous_amino(values[2])
        
        self.protein.set_score(self.best_score)


class RandomFolder():
    """
    The RandomFolder Class that randomly assigns possible folding values to every amino one by one.
    Starts over if folding results in dead end (= no more available values for an amino to fold).
    """
    def __init__(self, protein):
        self.protein = protein
        self.coordinates = set()
        self.folding = 0


    def fold_protein(self):
        """
        Folds protein by folding every amino one by one
        Starts over if amino folding failed due to no available values
        """
        length = len(self.protein.get_aminos())

        # Keep folding until protein is completely folded
        while self.folding < length - 1:
            
            # Fold amino and start over if amino folding failed
            if self.fold_amino() == 1:
                
                for amino in self.protein.get_aminos()[1:]:
                    
                    # Reset all values of aminos
                    amino.reset_amino()
                
                # Set fold position back to zero and clear coordinates
                self.folding = 0
                self.coordinates = set()

        # Calculate score and set in protein
        score = self.protein.calculate_score()
        self.protein.set_score(score)


    def fold_amino(self):
        """
        Fold amino by retrieving and setting fold  of current amino,
        and coordinate and previous amino of that amino's next neighbor.
        """
        # Store current amino and its coordinate
        current_amino = self.protein.get_aminos()[self.folding]
        self.coordinates.add(current_amino.get_coordinate())

        # Retrieve values for current amino and its neighbor
        new_values = self.get_values(current_amino)
        
        # return 1 if no values are returned
        if new_values is None:
            return 1

        # set retrieved values in current amino and its neighbor in protein
        self.set_values(new_values, self.folding)

        # update fold position
        self.folding += 1


    def get_values(self, amino):
        """
        Retrieves values for amino by randomly choosing values out of possible values
        """
        # Retrieve all possible values for amino
        possible_values = self.get_possible_values(amino)
        
        # Return None if no possible values
        if len(possible_values) == 0:
            return None
        
        # Choose random values
        fold, coordinate = random.choice(possible_values)
        
        return fold, coordinate

    
    def get_possible_values(self, amino):
        """
        Retrieves all possible values for an amino
        """    
        possible_values = []

        # Retrieve values for 2D-amino
        if len(amino.get_coordinate()) == 2:
            x, y = amino.get_coordinate()
        
            # Dict with all possible fold - next-coordinate combinations 
            values = {1: (x+1, y), -1: (x-1, y), 2:(x, y+1), -2: (x, y-1)} 

        # Retrieve values for 3D amino
        if len(amino.get_coordinate()) == 3:
            x, y, z = amino.get_coordinate()

            # Dict with all possible fold - next-coordinate combinations 
            values = {1: (x+1, y, z), -1: (x-1, y, z), 2:(x, y+1, z), -2: (x, y-1, z), 3: (x, y, z+1), -3: (x, y, z-1)} 

        # Add values to list if coordinate is not already occupied
        for fold, coordinate in values.items(): 
            if self.is_free_space(coordinate):
                possible_values.append((fold, coordinate))

        return possible_values


    def is_free_space(self, coordinate):
        """
        Checks whether coordinate is already in proteins occupied coordinates
        """
        # Return True if coordinate is not occupied, else False
        return coordinate not in self.coordinates


    def set_values(self, values, id):
        """
        Sets fold of amino with specified id, sets coordinate and 
        previous fold of that amino's next neighbor 
        """
        # Set fold of current amino
        self.protein.get_aminos()[id].set_fold(values[0])

        # Set coordinate and previous amino of next amino
        self.protein.get_aminos()[id + 1].set_coordinate(values[1])
        self.protein.get_aminos()[id + 1].set_previous_amino(-values[0])
        
        # Return score of new configuration
        return self.protein.calculate_score()