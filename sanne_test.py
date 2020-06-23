import matplotlib.pyplot as plt
from code import visualize
from code import Protein, Amino, Timer
from code import BestGreedy, BestOfRandom, HillClimber, SimulatedAnnealing
from code.algorithms.branch_and_bound import BranchBound

import profile, copy

if __name__ == "__main__":
    string = 'PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP'

    time = Timer()

    # intialize protein
    protein = Protein(string=string)

    bb = BranchBound(protein, p1=0.9, p2=0.6)

    time.start()

    bb.run()

    time.stop()

    print(sum([len(bb.all_scores[key]) for key in bb.all_scores]))

    print(f"Time: {time.get_time()} seconds")

    # print protein score
    if bb.protein.get_score():
        print(f"Score: {bb.protein.get_score()}")

    # plot protein
    visualize.visualize_protein(bb.protein)
    """
    # fold protein and set score with algorithm
    folder = HillClimber(greedy.protein, iterations=50000, mutations_per_iteration=4)
    folder.run()

    # print protein score
    if folder.protein.get_score():
        print(f"Score: {folder.protein.get_score()}")

    # plot protein
    visualize.visualize_protein(folder.protein)
    """