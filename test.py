import matplotlib.pyplot as plt

from code.classes import protein as prt, amino, protein_folder as prt_folder

if __name__ == "__main__":
    data_file = "data/proteins.txt"

    protein = prt.Protein(file=data_file)
    folder = prt_folder.ProteinFolder(protein)
    folder.fold()
    folder.calculate_score() # Note: calculate_score() could be implemented inside fold(), but this is probably more user-friendly
    
    if protein.get_score():
        print(protein.get_score())

    # visualize protein
    x = [amino.coordinate[0] for amino in protein.aminos]
    y = [amino.coordinate[1] for amino in protein.aminos]

    x_H = [amino.coordinate[0] for amino in protein.aminos if amino.type == 'H']
    y_H = [amino.coordinate[1] for amino in protein.aminos if amino.type == 'H']

    x_P = [amino.coordinate[0] for amino in protein.aminos if amino.type == 'P']
    y_P = [amino.coordinate[1] for amino in protein.aminos if amino.type == 'P']

    plt.figure()
    plt.xlim(-10, 10)
    plt.ylim(-10, 10)
    plt.xticks(range(-10, 10))
    plt.yticks(range(-10, 10))
    plt.grid(b=True)
    plt.plot(x, y, color='k')
    plt.scatter(x_H, y_H, color='r', marker='o')
    plt.scatter(x_P, y_P, color='b', marker='o')
    plt.show()