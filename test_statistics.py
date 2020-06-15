from code import Protein, Amino, Timer
from code import BestGreedy, BestOfRandom, HillClimber
import matplotlib.pyplot as plt
import csv

if __name__ == "__main__":

    time = Timer()

    scores = []
    length = []
    cysteine = []
    times = []

    # open source file and read first line
    with open('data/all_proteins.txt', 'r') as infile:
        proteins = infile.readlines()

    for i, line in enumerate(proteins):
        
        # intialize protein
        protein = Protein(string=line)

        time.start()

        # fold protein and set score with algorithm
        folder = BestOfRandom(protein, iterations=1000)
        folder.run()
        
        time.stop()

        # get score
        score = folder.protein.get_score()

        length.append(len(line))
        cysteine.append('C' in line)
        scores.append(score)
        times.append(time.get_time())


    with open('data/random_stats.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(length, cysteine, scores, times))


    """

    string = 'PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP'

    iterations_list = list(range(1, 1000, 10))
    scores = []

    for iterations in iterations_list:
        
        # intialize protein
        protein = Protein(string=string)

        # fold protein and set score with algorithm
        folder = BestGreedy(protein, iterations=iterations)
        folder.run()

        # add score of protein to list
        score = folder.protein.get_score()
        print(iterations, score)
        scores.append(score)

    # save stats
    with open('data/greedy_stats.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(iterations_list, scores))
    
    # plot scores against number of iterations
    plt.figure()
    plt.plot(iterations_list, scores)
    plt.title('Greedy algorithm (protein length: 36)')
    plt.xlabel('iterations')
    plt.ylabel('score')
    plt.show() 

    """  