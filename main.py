import csv
from code import visualize
from code import Protein, Amino, Timer
from code import BestGreedy, BestOfRandom, HillClimber, SimulatedAnnealing
from code import BFSPlus, BBProteinFolder

if __name__ == "__main__":

    time = Timer()

    protein_input = input("Type the protein string you want to fold: ")
    algorithm_input = int(input("""\nWhich algorithm do you want to use?
    Type '1' for Random.
    Type '2' for Greedy.
    Type '3' for Hillclimber.
    Type '4' for Simulated Annealing.
    Type '5' for Breadth First Search ++++
    Type '6' for Depth First Search with Branch and Bound\n"""))
    
    if algorithm_input in {1, 2, 3, 4}:
        d_input = int(input("""\nIn which dimensionality do you want to fold your protein? 
        Type '2' for 2-dimensional
        Type '3' for 3-dimensional.\n"""))
    else:
        d_input = 2

    protein = Protein(string=protein_input, dimensionality=d_input)

    ########
    #RANDOM#
    ########

    if algorithm_input == 1:
        
        iters = int(input("""\nHow many iterations do you want the Random algorithm to run?
        Type the amount.\n"""))

        print("\nRunning Random algorithm ...\n")

        time.start() 

        folder = BestOfRandom(protein, iterations=iters)
        folder.run()

        time.stop() 

    ########
    #GREEDY#
    ########

    elif algorithm_input == 2:
        
        iters = int(input("""\nHow many iterations do you want the Greedy algorithm to run?
        Type the amount.\n"""))

        print("\nRunning Greedy algorithm ...\n")

        time.start() 

        folder = BestGreedy(protein, iterations=iters)
        folder.run()

        time.stop() 

    ##############
    #HILL CLIMBER#
    ##############

    elif algorithm_input == 3:

        starting = int(input("""\nFor the Hillclimber algorithm, a starting state is required. 
        How do you want to generate a starting protein configuration?
        Type '1' for Random.
        Type '2' for Greedy.\n"""))

        if starting == 1:

            starting_iters = int(input("""\nHow many iterations do you want the Random algorithm to run?
            Type the amount.\n"""))
            
            print("\nGenerating starting state ...\n")

            starting_folder = BestOfRandom(protein, iterations=starting_iters)
            starting_folder.run()

        if starting == 2:

            starting_iters = int(input("""\nHow many iterations do you want the Greedy algorithm to run?
            Type the amount.\n"""))
            
            print("\nGenerating starting state ...\n")

            starting_folder = BestGreedy(protein, iterations=starting_iters)
            starting_folder.run()

        print(f"Score of starting state is: {starting_folder.protein.get_score()}\n")
        
        iters = int(input("""How many iterations do you want the Hillclimber algorithm to run?
        Type the amount.\n"""))
        muts = int(input("""\nHow many mutations do you want the algorithm to execute per iteration?
        Type the amount.\n"""))

        print("\nRunning Hillclimber algorithm ...\n")

        time.start() 

        folder = HillClimber(starting_folder.protein, iterations=iters, mutations_per_iteration=muts)
        folder.run()

        time.stop() 


    #####################
    #SIMULATED ANNEALING#
    #####################

    elif algorithm_input == 4:

        starting = int(input("""\nFor the Simulated Annealing algorithm, a starting state is required. 
        How do you want to generate a starting protein configuration?
        Type '1' for Random.
        Type '2' for Greedy.\n"""))

        if starting == 1:

            starting_iters = int(input("""\nHow many iterations do you want the Random algorithm to run?
            Type the amount.\n"""))
            
            print("\nGenerating starting state ...\n")

            starting_folder = BestOfRandom(protein, iterations=starting_iters)
            starting_folder.run()

        if starting == 2:

            starting_iters = int(input("""\nHow many iterations do you want the Greedy algorithm to run?
            Type the amount.\n"""))
            
            print("\nGenerating starting state ...\n")

            starting_folder = BestGreedy(protein, iterations=starting_iters)
            starting_folder.run()

        print(f"Score of starting state is: {starting_folder.protein.get_score()}\n")
        
        iters = int(input("""How many iterations do you want the Simulated Annealing algorithm to run?
        Type the amount.\n"""))
        muts = int(input("""\nHow many mutations do you want the algorithm to execute per iteration?
        Type the amount.\n"""))
        temp = float(input("""\nWhat starting temperature do you want the algorithm to use?
        Type the amount of degrees.\n"""))

        print("\nRunning Simulated Annealing algorithm ...\n")

        time.start() 

        folder = SimulatedAnnealing(starting_folder.protein, iterations=iters, mutations_per_iteration=muts, temperature=temp)
        folder.run()

        time.stop() 

    ###########################
    #BREADTH FIRST SEARCH ++++#
    ###########################

    elif algorithm_input == 5:

        default = int(input("""\nDo you want to use default parameters or advanced parameters?
        Type 1 for default.
        Type 2 for advanced.\n"""))

        if default == 2:
            pruning_depth = int(input("""\nFrom what depth do you want to start pruning? (Suggested: 8)\n"""))
            pruning_distance = int(input("""\nWhat is the value of the pruning distance factor? (Suggested: 4)\n"""))
            queue_size = int(input("""\nWhat maximal size for the queue? (Suggested: 2000)\n"""))

        print("\nRunning Breadth First Search ++++ algorithm ...\n")

        time.start() 

        if default == 2:
            folder = BFSPlus(protein, pruning_depth, pruning_distance, queue_size)
        else:
            folder = BFSPlus(protein)

        folder.fold()
        folder.set_score()

        time.stop() 

    ##########################################
    #DEPTH FIRST SEARCH WITH BRANCH AND BOUND#
    ##########################################

    elif algorithm_input == 6:

        print("\nRunning Depth First Search with Branch and Bound algorithm ...\n")

        time.start() 

        folder = BBProteinFolder(protein)
        folder.fold()
        folder.set_score()

        time.stop() 

    #########
    #RESULTS#
    #########
    
    print(f"Runtime: {time.get_time()} seconds")

    # print protein score
    if folder.finished_folded_protein.get_score():
        print(f"Score: {folder.finished_folded_protein.get_score()}")

    # write protein output to csv-file
    path = input("""\nType the path to a csv output file. For example: 'data/output.csv'\n""")

    with open(path, 'w') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['amino', 'fold'])
        for amino in folder.finished_folded_protein.get_aminos():
            csv_writer.writerow([f'{amino.type}', f'{amino.fold}'])
        csv_writer.writerow(['score', f'{folder.finished_folded_protein.get_score()}'])

    print("\nWrote output file. Plotting protein ...")

    # plot protein
    visualize.visualize_protein(folder.finished_folded_protein)
