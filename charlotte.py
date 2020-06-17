import csv
import profile
import pickle
import os

from code.visualisations import visualize
from code.classes import protein as prt, amino, timer
from code.algorithms import charlotte_protein_folder as prt_folder

if __name__ == "__main__":
    data_file = "data/proteins.txt"

    pickled = {}
    pickle_path = "./charlotte_pickles/pickle" + str(len(os.listdir("./charlotte_pickles")))

    # intialize protein
    protein = prt.Protein(file=data_file)

    time = timer.Timer()
    time.start()

    # fold protein with algorithm
    folder = prt_folder.CharlotteProteinFolder(protein)
    # profile.run('folder.fold()')
    folder.fold()
    folder.set_score()
    
    time.stop()

    print(f"Time: {time.get_time()} seconds")

    pickled['protein'] = protein
    pickled['folder'] = folder
    pickled['score'] = protein.score
    pickled['time'] = time.get_time()
    pickle.dump(pickled, pickle_path)

    # write protein output to csv-file 
    with open('data/output.csv', 'w') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['amino', 'fold'])
        for amino in folder.finished_folded_protein.get_aminos():
            csv_writer.writerow([f'{amino.type}', f'{amino.fold}'])
        csv_writer.writerow(['score', f'{protein.get_score()}'])

    # print protein score
    if folder.finished_folded_protein.get_score():
        print(f"Score: {folder.finished_folded_protein.get_score()}")
        print(f"Finished protein: {folder.finished_folded_protein.to_string_with_previous()}")

    # plot protein
    visualize.visualize_protein(folder.finished_folded_protein)
