import csv
import profile

from code.visualisations import visualize
from code.classes import protein as prt, amino, timer
from code.algorithms import charlotte_protein_folder_3d as prt_folder
from code.algorithms.help_methods.pickling import *

if __name__ == "__main__":
    data_file = "data/proteins.txt"

    # intialize protein
    protein = prt.Protein(file=data_file, dimensionality=3)

    time = timer.Timer()
    time.start()

    # fold protein with algorithm
    folder = prt_folder.CharlotteProteinFolder3D(protein)
    # profile.run('folder.fold()')
    folder.fold()
    folder.set_score()
    
    time.stop()

    print(f"Time: {time.get_time()} seconds")

    # Pickling
    pickled = {}
    pickled['protein'] = protein.source_string
    pickled['folder'] = folder
    pickled['pruning_depth'] = folder.pruning_depth
    pickled['pruning_distance_heur'] = folder.pruning_distance_heur
    pickled['relevance_score_heur'] = folder.relevance_score_heur
    pickled['max_queue_size'] = folder.max_queue_size
    pickled['score'] = folder.finished_folded_protein.get_score()
    pickled['time'] = time.get_time()
    saved_file = store_pickle("charlotte_pickles", pickled)
    
    show_pickle(saved_file)

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
