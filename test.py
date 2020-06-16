import csv
import profile
from code import Protein, Amino, Timer, visualize
from code.algorithms import faster_random as rpf

if __name__ == "__main__":
    data_file = "data/proteins.txt"

    # intialize protein
    protein = Protein(file=data_file)

    time = Timer()
    time.start()

    # fold protein and set score with algorithm
    protein = rpf.fold(protein, iterations=5000)
    
    time.stop()

    print(f"Time: {time.get_time()} seconds")

    # write protein output to csv-file
    with open('data/output.csv', 'w') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['amino', 'fold'])
        for amino in protein.get_aminos():
            csv_writer.writerow([f'{amino.type}', f'{amino.fold}'])
        csv_writer.writerow(['score', f'{protein.get_score()}'])

    # print protein score
    if protein.get_score():
        print(f"Score: {protein.get_score()}")
        print(f"Finished protein: {protein.to_string_with_previous()}")


    # plot protein
    visualize.visualize_protein(protein)
