from code import Protein, Amino, Timer
from code.algorithms import random_protein_folder as rpf
import matplotlib.pyplot as plt

if __name__ == "__main__":
    string = 'PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP'

    iterations_list = list(range(1, 100000, 1000))
    scores = []

    for iterations in iterations_list:
        # intialize protein
        protein = Protein(string=string)

        # fold protein and set score with algorithm
        protein = rpf.fold(protein, iterations=iterations)

        # add score of protein to list
        scores.append(protein.get_score())

    # plot scores against number of iterations
    plt.figure()
    plt.plot(iterations_list, scores)
    plt.title('Random algorithm (protein length: 36)')
    plt.xlabel('iterations')
    plt.ylabel('score')
    plt.show()   