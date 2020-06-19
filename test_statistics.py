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
    iter_list = []

    # open source file and read first line
    with open('data/all_proteins.txt', 'r') as infile:
        proteins = infile.readlines()

    for i, line in enumerate(proteins):
        print(i)
        line = line.replace('\n', '')
        
        """
        RANDOM
        """
        random2d_dic = {0:400000, 1:200000, 2:50000, 3:20000, 4:60000, 5:50000, 6:20000, 7:20000}
        algorithm = "Random"
        
        # intialize protein
        protein = Protein(string=line)

        # fold protein and set score with algorithm
        folder = BestOfRandom(protein, iterations=random2d_dic[i])

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
        iter_list.append(random2d_dic[i])

        """
        GREEDY
        """
        greedy2d_dic = {0:100000, 1:60000, 2:20000, 3:6000, 4:20000, 5:20000, 6:8000, 7:8000}
        algorithm = "Greedy"
        
        # intialize protein
        protein = Protein(string=line)

        # fold protein and set score with algorithm
        folder = BestGreedy(protein, iterations=greedy2d_dic[i])

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
        iter_list.append(greedy2d_dic[i])

        """
        HILLCLIMBER
        """
        hillclimber2d_dic = {0:600000, 1:350000, 2:160000, 3:100000, 4:260000, 5:120000, 6:70000, 7:80000}
        algorithm = "Hillclimber"
        
        # intialize protein
        protein = Protein(string=line)

        # fold protein and set score with algorithm
        greedy = BestGreedy(protein, iterations=1000)
        greedy.run()
        folder = HillClimber(greedy.protein, iterations=hillclimber2d_dic[i], mutations_per_iteration=3)

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
        iter_list.append(hillclimber2d_dic[i])

        """
        SIMULATING ANNEALING
        """
        annealing2d_dic = {0:550000, 1:410000, 2:210000, 3:140000, 4:240000, 5:180000, 6:150000, 7:170000}
        algorithm = "Simulating Annealing"
        
        # intialize protein
        protein = Protein(string=line)

        # fold protein and set score with algorithm
        greedy = BestGreedy(protein, iterations=1000)
        greedy.run()
        folder = SimulatedAnnealing(greedy.protein, iterations=annealing2d_dic[i], mutations_per_iteration=3, temperature=2)

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
        iter_list.append(annealing2d_dic[i])    

        """
        RANDOM 3D
        """
        random3d_dic = {0:210000, 1:130000, 2:30000, 3:10000, 4:40000, 5:30000, 6:10000, 7:10000}
        algorithm = "Random"
        
        # intialize protein
        protein = Protein(string=line, dimensionality=3)

        # fold protein and set score with algorithm
        folder = BestOfRandom(protein, iterations=random3d_dic[i])

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
        iter_list.append(random3d_dic[i])

        """
        GREEDY 3D
        """
        greedy3d_dic = {0:40000, 1:20000, 2:8000, 3:3000, 4:8000, 5:7000, 6:3000, 7:3000}
        algorithm = "Greedy"
        
        # intialize protein
        protein = Protein(string=line, dimensionality=3)

        # fold protein and set score with algorithm
        folder = BestGreedy(protein, iterations=greedy3d_dic[i])

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
        iter_list.append(greedy3d_dic[i])

        """
        HILLCLIMBER 3D
        """
        hillclimber3d_dic = {0:660000, 1:420000, 2:230000, 3:140000, 4:230000, 5:210000, 6:110000, 7:120000}
        algorithm = "Hillclimber"
        
        # intialize protein
        protein = Protein(string=line, dimensionality=3)

        # fold protein and set score with algorithm
        greedy = BestGreedy(protein, iterations=1000)
        greedy.run()
        folder = HillClimber(greedy.protein, iterations=hillclimber3d_dic[i], mutations_per_iteration=3)

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
        iter_list.append(hillclimber3d_dic[i])

        """
        SIMULATING ANNEALING 3D
        """
        annealing3d_dic = {0:700000, 1:430000, 2:240000, 3:140000, 4:260000, 5:210000, 6:110000, 7:110000}
        algorithm = "Simulating Annealing"
        
        # intialize protein
        protein = Protein(string=line, dimensionality=3)

        # fold protein and set score with algorithm
        greedy = BestGreedy(protein, iterations=1000)
        greedy.run()
        folder = SimulatedAnnealing(greedy.protein, iterations=annealing3d_dic[i], mutations_per_iteration=3, temperature=2)

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
        iter_list.append(annealing3d_dic[i])


    with open('data/all_stats.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['algorithm', 'length', 'cysteine', 'dimensionality', 'iterations', 'scores', 'runtime'])
        writer.writerows(zip(algorithms, length, cysteine, dimensionalities, iter_list, scores, times))