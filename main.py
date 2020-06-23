import csv
from PIL import Image
from code import visualize
from code import Protein, Amino, Timer
from code import BestGreedy, BestOfRandom, HillClimber, SimulatedAnnealing, BranchBound, BFSplus

if __name__ == "__main__":

    time = Timer()

    protein_input = input("Type the protein string you want to fold: ")

    if len(protein_input) < 36:
        algorithm_input = int(input("""\nWhich algorithm do you want to use?
        Type '1' for Random.
        Type '2' for Greedy.
        Type '3' for Hillclimber.
        Type '4' for Simulated Annealing.
        Type '5' for Breadth First Search ++++
        Type '6' for Branch and Bound\n"""))
    
    else:
        algorithm_input = int(input("""\nWhich algorithm do you want to use?
        Type '1' for Random.
        Type '2' for Greedy.
        Type '3' for Hillclimber.
        Type '4' for Simulated Annealing.
        Type '5' for Breadth First Search ++++\n"""))

    if algorithm_input in {1, 2, 3, 4, 5}:
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
        
        image = Image.open('data/iterations_random.png')
        image.show()

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
        
        image = Image.open('data/iterations_greedy.png')
        image.show()

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

            image = Image.open('data/iterations_random.png')
            image.show()

            starting_iters = int(input("""\nHow many iterations do you want the Random algorithm to run?
            Type the amount.\n"""))
            
            print("\nGenerating starting state ...\n")

            starting_folder = BestOfRandom(protein, iterations=starting_iters)
            starting_folder.run()

        if starting == 2:

            image = Image.open('data/iterations_greedy.png')
            image.show()

            starting_iters = int(input("""\nHow many iterations do you want the Greedy algorithm to run?
            Type the amount.\n"""))
            
            print("\nGenerating starting state ...\n")

            starting_folder = BestGreedy(protein, iterations=starting_iters)
            starting_folder.run()

        print(f"Score of starting state is: {starting_folder.protein.get_score()}\n")
        
        image = Image.open('data/iterations_hill_annealing.png')
        image.show()

        iters = int(input("""How many iterations do you want the Hillclimber algorithm to run?
        Type the amount.\n"""))
        muts = int(input("""\nHow many mutations do you want the algorithm to execute per iteration? (Suggested: 4)
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

            image = Image.open('data/iterations_random.png')
            image.show()

            starting_iters = int(input("""\nHow many iterations do you want the Random algorithm to run?
            Type the amount.\n"""))
            
            print("\nGenerating starting state ...\n")

            starting_folder = BestOfRandom(protein, iterations=starting_iters)
            starting_folder.run()

        if starting == 2:

            image = Image.open('data/iterations_greedy.png')
            image.show()

            starting_iters = int(input("""\nHow many iterations do you want the Greedy algorithm to run?
            Type the amount.\n"""))
            
            print("\nGenerating starting state ...\n")

            starting_folder = BestGreedy(protein, iterations=starting_iters)
            starting_folder.run()

        print(f"Score of starting state is: {starting_folder.protein.get_score()}\n")

        image = Image.open('data/iterations_hill_annealing.png')
        image.show()

        iters = int(input("""How many iterations do you want the Simulated Annealing algorithm to run?
        Type the amount.\n"""))
        muts = int(input("""\nHow many mutations do you want the algorithm to execute per iteration? (Suggested: 4)
        Type the amount.\n"""))
        temp = float(input("""\nWhat starting temperature do you want the algorithm to use? (Suggested: 2)
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

        param_input = int(input("""\nDo you want to use default parameters or advanced parameters?
        Type '1' for default.
        Type '2' for advanced.\n"""))

        if param_input == 2:
            if d_input == 2:
                pruning_depth = int(input("""\nFrom what depth do you want to start pruning? (Suggested: 8)\n"""))
                pruning_distance = int(input("""\nWhat is the value of the pruning distance factor? (Suggested: 4)\n"""))
                queue_size = int(input("""\nWhat maximal size for the queue? (Suggested: 2000)\n"""))
            else:
                pruning_depth = int(input("""\nFrom what depth do you want to start pruning? (Suggested: 4)\n"""))
                pruning_distance = int(input("""\nWhat is the value of the pruning distance factor? (Suggested: 6)\n"""))
                queue_size = int(input("""\nWhat maximal size for the queue? (Suggested: 10000)\n"""))


        print("\nRunning Breadth First Search ++++ algorithm ...\n")

        time.start() 

        if param_input == 2:
            BFSPlus(protein, d_input, pruning_depth, pruning_distance, queue_size)
        else:
            folder = BFSPlus(protein, d_input)

        folder.fold()
        folder.set_score()

        time.stop() 

    ##################
    #BRANCH AND BOUND#
    ##################

    elif algorithm_input == 6:
        
        prob1 = float(input("""\nWith with probability do you want the algorithm to discard a protein configuration if the score of that configuration is worse than the mean score? (Suggested: 0.8) \n"""))
        prob2 = float(input("""\nWith with probability do you want the algorithm to discard a protein configuration if the score of that configuration is better than the mean score but worse than the best score? (Suggested: 0.5) \n"""))
        print("\nRunning Branch and Bound algorithm ...\n")

        time.start() 

        folder = BranchBound(protein, p1=prob1, p2=prob2)
        folder.run()

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
