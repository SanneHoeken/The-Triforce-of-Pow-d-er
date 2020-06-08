import matplotlib.pyplot as plt

def visualize_protein(protein):
    
    x = [amino.coordinate[0] for amino in protein.aminos]
    y = [amino.coordinate[1] for amino in protein.aminos]

    x_H = [amino.coordinate[0] for amino in protein.aminos if amino.type == 'H']
    y_H = [amino.coordinate[1] for amino in protein.aminos if amino.type == 'H']

    x_P = [amino.coordinate[0] for amino in protein.aminos if amino.type == 'P']
    y_P = [amino.coordinate[1] for amino in protein.aminos if amino.type == 'P']

    x_C = [amino.coordinate[0] for amino in protein.aminos if amino.type == 'C']
    y_C = [amino.coordinate[1] for amino in protein.aminos if amino.type == 'C']

    plt.figure()
    plt.xlim(-4, 4)
    plt.ylim(-4, 4)
    plt.xticks(range(-4, 4))
    plt.yticks(range(-4, 4))
    plt.grid(b=True)
    plt.plot(x, y, color='k')
    plt.scatter(x_H, y_H, color='r', marker='o')
    plt.scatter(x_P, y_P, color='b', marker='o')
    plt.scatter(x_C, y_C, color='g', marker='o')
    plt.show()