import matplotlib.pyplot as plt

from code.classes import protein

if __name__ == "__main__":
    data_file = "data/proteins.txt"
    protein = protein.Protein(data_file)

    print(protein.score)

    # visualize protein
    x = [amino.coordinate[0] for amino in protein.aminos]
    y = [amino.coordinate[1] for amino in protein.aminos]

    x_H = []
    y_H = []

    plt.figure()
    plt.xlim(-10, 10)
    plt.ylim(-10, 10)
    plt.xticks(range(-10, 10))
    plt.yticks(range(-10, 10))
    plt.grid(b=True)
    plt.plot(x, y)
    plt.scatter(x_H, y_H)
    plt.show()