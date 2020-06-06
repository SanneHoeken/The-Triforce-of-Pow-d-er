import matplotlib.pyplot as plt

from code.classes import protein

if __name__ == "__main__":
    data_file = "data/proteins.txt"
    protein = protein.Protein(file=data_file)

    if protein.score:
        print(protein.score)

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