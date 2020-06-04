import matplotlib.pyplot as plt

from code.classes import protein

if __name__ == "__main__":
    data_file = "data/proteins.txt"
    protein = protein.Protein(data_file)

    x = [amino.coordinate[0] for amino in protein.aminos]
    y = [amino.coordinate[1] for amino in protein.aminos]

    plt.figure()
    plt.xlim(-5, 5)
    plt.ylim(-5, 5)
    plt.xticks(range(-5, 5))
    plt.yticks(range(-5, 5))
    plt.grid(b=True)
    plt.plot(x, y)
    plt.show()