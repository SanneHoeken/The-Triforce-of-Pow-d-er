from code import Protein, Amino, Timer
from code import BestGreedy, BestOfRandom, HillClimber, SimulatedAnnealing
import matplotlib.pyplot as plt
import csv

if __name__ == "__main__":
    
    time = Timer()

    algorithms = []
    scores = []
    length = []
    cysteine = []
    times = []
    dimensionalities = []

    # open source file and read first line
    with open('data/all_proteins.txt', 'r') as infile:
        proteins = infile.readlines()

    for i, line in enumerate(proteins):
        
        line = line.replace('\n', '')

        """
        RANDOM
        """
        algorithm = "Random"
        
        # intialize protein
        protein = Protein(string=line)

        # fold protein and set score with algorithm
        folder = BestOfRandom(protein, iterations=1000)

        time.start() 

        folder.run()

        time.stop() 

        # get score
        score = folder.protein.get_score()

        algorithms.append(algorithm)
        length.append(len(line))
        cysteine.append('with C' if 'C' in line else 'without C')
        dimensionalities.append(2)
        scores.append(score)
        times.append(time.get_time())

        """
        GREEDY
        """
        algorithm = "Greedy"
        
        # intialize protein
        protein = Protein(string=line)

        # fold protein and set score with algorithm
        folder = BestGreedy(protein, iterations=1000)

        time.start() 

        folder.run()

        time.stop() 

        # get score
        score = folder.protein.get_score()

        algorithms.append(algorithm)
        length.append(len(line))
        cysteine.append('with C' if 'C' in line else 'without C')
        dimensionalities.append(2)
        scores.append(score)
        times.append(time.get_time())

        """
        HILLCLIMBER
        """
        algorithm = "Hillclimber"
        
        # intialize protein
        protein = Protein(string=line)

        # fold protein and set score with algorithm
        greedy = BestGreedy(protein, iterations=1000)
        greedy.run()
        folder = HillClimber(greedy.protein, iterations=10000, mutations_per_iteration=3)

        time.start() 

        folder.run()

        time.stop() 

        # get score
        score = folder.protein.get_score()

        algorithms.append(algorithm)
        length.append(len(line))
        cysteine.append('with C' if 'C' in line else 'without C')
        dimensionalities.append(2)
        scores.append(score)
        times.append(time.get_time())

        """
        SIMULATING ANNEALING
        """
        algorithm = "Simulating Annealing"
        
        # intialize protein
        protein = Protein(string=line)

        # fold protein and set score with algorithm
        greedy = BestGreedy(protein, iterations=1000)
        greedy.run()
        folder = SimulatedAnnealing(greedy.protein, iterations=10000, mutations_per_iteration=3, temperature=1)

        time.start() 

        folder.run()

        time.stop() 

        # get score
        score = folder.protein.get_score()

        algorithms.append(algorithm)
        length.append(len(line))
        cysteine.append('with C' if 'C' in line else 'without C')
        dimensionalities.append(2)
        scores.append(score)
        times.append(time.get_time())

        """
        RANDOM 3D
        """
        algorithm = "Random"
        
        # intialize protein
        protein = Protein(string=line, dimensionality=3)

        # fold protein and set score with algorithm
        folder = BestOfRandom(protein, iterations=1000)

        time.start() 

        folder.run()

        time.stop() 

        # get score
        score = folder.protein.get_score()

        algorithms.append(algorithm)
        length.append(len(line))
        cysteine.append('with C' if 'C' in line else 'without C')
        dimensionalities.append(3)
        scores.append(score)
        times.append(time.get_time())

        """
        GREEDY 3D
        """
        algorithm = "Greedy"
        
        # intialize protein
        protein = Protein(string=line, dimensionality=3)

        # fold protein and set score with algorithm
        folder = BestGreedy(protein, iterations=1000)

        time.start() 

        folder.run()

        time.stop() 

        # get score
        score = folder.protein.get_score()

        algorithms.append(algorithm)
        length.append(len(line))
        cysteine.append('with C' if 'C' in line else 'without C')
        dimensionalities.append(3)
        scores.append(score)
        times.append(time.get_time())

        """
        HILLCLIMBER 3D
        """
        algorithm = "Hillclimber"
        
        # intialize protein
        protein = Protein(string=line, dimensionality=3)

        # fold protein and set score with algorithm
        greedy = BestGreedy(protein, iterations=1000)
        greedy.run()
        folder = HillClimber(greedy.protein, iterations=10000, mutations_per_iteration=3)

        time.start() 

        folder.run()

        time.stop() 

        # get score
        score = folder.protein.get_score()

        algorithms.append(algorithm)
        length.append(len(line))
        cysteine.append('with C' if 'C' in line else 'without C')
        dimensionalities.append(3)
        scores.append(score)
        times.append(time.get_time())

        """
        SIMULATING ANNEALING 3D
        """
        algorithm = "Simulating Annealing"
        
        # intialize protein
        protein = Protein(string=line, dimensionality=3)

        # fold protein and set score with algorithm
        greedy = BestGreedy(protein, iterations=1000)
        greedy.run()
        folder = SimulatedAnnealing(greedy.protein, iterations=10000, mutations_per_iteration=3, temperature=1)

        time.start() 

        folder.run()

        time.stop() 

        # get score
        score = folder.protein.get_score()

        algorithms.append(algorithm)
        length.append(len(line))
        cysteine.append('with C' if 'C' in line else 'without C')
        dimensionalities.append(3)
        scores.append(score)
        times.append(time.get_time())


    with open('data/all_stats.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['algorithms', 'length', 'cysteine', 'dimensionality', 'scores', 'runtime'])
        writer.writerows(zip(algorithms, length, cysteine, dimensionalities, scores, times))