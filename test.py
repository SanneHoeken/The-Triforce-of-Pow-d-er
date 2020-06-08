import matplotlib.pyplot as plt
import csv

from code.classes import protein as prt, amino, random_protein_folder as prt_folder, timer

if __name__ == "__main__":
    data_file = "data/proteins.txt"

    # intialize protein
    protein = prt.Protein(file=data_file)

    time = timer.Timer()
    time.start()

    # fold protein and calculate score
    folder = prt_folder.RandomProteinFolder(protein)
    folder.fold()
    folder.calculate_score() # Note: calculate_score() could be implemented inside fold(), but this is probably more user-friendly
    
    time.stop()

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

    # visualize protein
    x = [amino.coordinate[0] for amino in protein.aminos]
    y = [amino.coordinate[1] for amino in protein.aminos]

    x_H = [amino.coordinate[0] for amino in protein.aminos if amino.type == 'H']
    y_H = [amino.coordinate[1] for amino in protein.aminos if amino.type == 'H']

    x_P = [amino.coordinate[0] for amino in protein.aminos if amino.type == 'P']
    y_P = [amino.coordinate[1] for amino in protein.aminos if amino.type == 'P']

    x_C = [amino.coordinate[0] for amino in protein.aminos if amino.type == 'C']
    y_C = [amino.coordinate[1] for amino in protein.aminos if amino.type == 'C']

    plt.figure()
    plt.xlim(-10, 10)
    plt.ylim(-10, 10)
    plt.xticks(range(-10, 10))
    plt.yticks(range(-10, 10))
    plt.grid(b=True)
    plt.plot(x, y, color='k')
    plt.scatter(x_H, y_H, color='r', marker='o')
    plt.scatter(x_P, y_P, color='b', marker='o')
    plt.scatter(x_C, y_C, color='g', marker='o')
    plt.show()