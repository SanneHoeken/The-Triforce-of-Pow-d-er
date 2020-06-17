import csv
from code import visualize
from code import Protein, Amino, Timer
from code import BestGreedy, BestOfRandom, HillClimber, SimulatedAnnealing

if __name__ == "__main__":

    time = Timer()

    protein_input = input("Protein string: ")
    algorithm_input = input("Algorithm: ")
    d_input = int(input("Dimensionality: "))

    protein = Protein(string=protein_input, dimensionality=d_input)

    """
    RANDOM
    """
    if algorithm_input == "random":
        
        iters = int(input("Iterations: "))

        print("Running algorithm ...")

        time.start() 

        folder = BestOfRandom(protein, iterations=iters)
        folder.run()

        time.stop() 

    """
    GREEDY
    """
    if algorithm_input == "greedy":
        
        iters = int(input("Iterations: "))

        print("Running algorithm ...")

        time.start() 

        folder = BestGreedy(protein, iterations=iters)
        folder.run()

        time.stop() 

    """
    HILL CLIMBER
    """

    if algorithm_input == "hillclimber":


        print("Generating starting state ...")
        
        greedy_folder = BestGreedy(protein, iterations=500)
        greedy_folder.run()

        print("Starting state is best out of 500 greedy folded proteins.")
        print(f"Score of starting state is: {greedy_folder.protein.get_score()}")
        
        iters = int(input("Iterations for Hillclimber: "))
        muts = int(input("Mutations per iteration: "))

        print("Running algorithm ...")

        time.start() 

        folder = HillClimber(greedy_folder.protein, iterations=iters, mutations_per_iteration=muts)
        folder.run()

        time.stop() 


    """
    SIMULATED ANNEALING
    """

    if algorithm_input == "simulated annealing":

        print("Generating starting state ...")
        
        greedy_folder = BestGreedy(protein, iterations=500)
        greedy_folder.run()

        print("Starting state is best out of 500 greedy folded proteins.")
        print(f"Score of starting state is: {greedy_folder.protein.get_score()}")

        iters = int(input("Iterations for Simulated Annealing: "))
        muts = int(input("Mutations per iteration: "))
        temp = float(input("Temperature: "))

        print("Running algorithm ...")

        time.start() 

        folder = SimulatedAnnealing(greedy_folder.protein, iterations=iters, mutations_per_iteration=muts, temperature=temp)
        folder.run()

        time.stop() 


    """
    RESULTS
    """
    
    print(f"Time: {time.get_time()} seconds")

    # print protein score
    if folder.protein.get_score():
        print(f"Score: {folder.protein.get_score()}")

    # write protein output to csv-file
    path = input("Path to csv output file: ")

    with open(path, 'w') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['amino', 'fold'])
        for amino in protein.get_aminos():
            csv_writer.writerow([f'{amino.type}', f'{amino.fold}'])
        csv_writer.writerow(['score', f'{protein.get_score()}'])

    print("Wrote output file. Plotting protein ...")

    # plot protein
    visualize.visualize_protein(folder.protein)
