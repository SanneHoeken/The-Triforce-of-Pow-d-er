import csv
import profile

from code.visualisations import visualize
from code.classes import protein as prt, amino, timer
from code.algorithms import thomas_protein_folder as prt_folder

if __name__ == "__main__":
    data_file = "data/proteins.txt"

    # intialize protein
    protein = prt.Protein(file=data_file)

    time = timer.Timer()
    time.start()

    # fold protein with algorithm
    folder = prt_folder.BBProteinFolder(protein)
    #profile.run('folder.fold()')
    folder.fold()
    folder.set_score()
    
    time.stop()

    print(f"Time: {time.get_time()} seconds")


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
