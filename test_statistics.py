from code import Protein, Amino, Timer
from code import BestGreedy
import matplotlib.pyplot as plt
import csv

if __name__ == "__main__":
    string = 'PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP'

    iterations_list = list(range(1, 20000, 100))
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