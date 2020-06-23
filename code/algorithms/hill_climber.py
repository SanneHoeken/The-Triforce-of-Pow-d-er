import random
from code import Amino, Protein

class HillClimber():
    """
    The HillClimber class that mutates a protein i times by changing n random amino folds. 
    Each improvement is kept for the next iteration if the mutated protein is valid.
    """
    def __init__(self, protein, iterations=1, mutations_per_iteration=1):
        self.protein = protein
        self.protein_length = len(self.protein.get_aminos())
        self.iterations = iterations
        self.mutations_per_iteration = mutations_per_iteration
        self.archive = []
        self.tmp_archive = []
        self.best_score = self.protein.get_score()
        self.finished_folded_protein = self.protein
    

    def run(self):
        """
        Runs the Hillclimber algorithm for a specific amount of iterations with
        a specified amount of mutations per iteration.
        """
        
        for i in range(self.iterations):

            # Store values of current protein configuration
            self.archive = [(amino.fold, amino.coordinate, amino.previous_amino) for amino in self.protein.get_aminos()]
            
            # Consecultively mutates protein specified times
            mutation_count = 0
            
            while mutation_count < self.mutations_per_iteration:

                # Store values of current protein configuration
                self.tmp_archive = [(amino.fold, amino.coordinate, amino.previous_amino) for amino in self.protein.get_aminos()]
                
                # Mutate protein
                self.mutate()

                # Undo mutation if new configuration is not valid
                if not self.is_valid_protein():
                    self.undo_mutation(self.tmp_archive)
                else:
                    mutation_count += 1

            # Update best score if new configuration is better, else undo mutation series
            score = self.protein.calculate_score()
            if score < self.best_score:
                self.best_score = score
            else:
                self.undo_mutation(self.archive)
            
        # Set score of protein
        self.protein.set_score(self.best_score)


    def mutate(self):
        """
        Mutates a protein by changing a random amino fold and adjusting all
        amino values in the protein to that change.
        """
        # Get mutation values: position of amino and new fold
        pos, fold = self.generate_mutation()
        
        co = (None, None)

        # Iterates over every amino starting at the mutation's position
        for amino in self.protein.get_aminos()[pos: self.protein_length]:
            
            # Change fold for selected amino
            if amino.id == pos:
                amino.set_fold(fold)

                # Update coordinate
                co = amino.coordinate

            # Change coordinates and previous amino for amino after selected amino
            elif amino.id == pos + 1:
                new_co = self.protein.calculate_coordinate(fold, co)
                amino.set_previous_amino(-fold)
                amino.set_coordinate(new_co)

                # Update fold and coordinate
                fold = amino.fold
                co = new_co

            # Change coordinates for every amino after selected amino
            elif amino.id > pos + 1:
                new_co = self.protein.calculate_coordinate(fold, co)
                amino.set_coordinate(new_co)

                # Update fold and coordinate
                fold = amino.fold
                co = new_co
    

    def generate_mutation(self):
        """
        Retrieves randomly an amino position and a new fold
        """
        # Choose random position of amino
        position = random.randint(0, self.protein_length - 2)

        # Set new fold to amino's current fold
        new_fold = self.protein.get_aminos()[position].fold

        # Retrieve new fold if protein is 2D
        if len(self.protein.get_aminos()[0].coordinate) == 2: 
            
            # Choose new fold randomly until amino's fold is changed
            while self.protein.get_aminos()[position].fold == new_fold:
                new_fold = random.choice([-2, 2, -1, 1])

        # Retrieve new fold if protein is 3D
        elif len(self.protein.get_aminos()[0].coordinate) == 3: 
            
            # Choose new fold randomly until amino's fold is changed
            while self.protein.get_aminos()[position].fold == new_fold:
                new_fold = random.choice([-3, 3, -2, 2, -1, 1])

        return position, new_fold


    def is_valid_protein(self):
        """
        Checks whether protein configuration is valid
        """
        # Returns True if no double coordinates in protein
        coordinates = [amino.coordinate for amino in self.protein.get_aminos()]
        
        return len(set(coordinates)) == len(coordinates)


    def undo_mutation(self, previous_configuration):
        """
        Sets all values of protein to values of specified previous configuration
        """
        for i, values in enumerate(previous_configuration):
            self.protein.get_aminos()[i].set_fold(values[0])
            self.protein.get_aminos()[i].set_coordinate(values[1])
            self.protein.get_aminos()[i].set_previous_amino(values[2])
