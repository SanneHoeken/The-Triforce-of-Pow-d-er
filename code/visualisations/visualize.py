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

    x_range = round(0.5 * len(protein.aminos))
    y_range = round(0.5 * len(protein.aminos))

    plt.figure()
    plt.xlim(-x_range, x_range)
    plt.ylim(-y_range, y_range)
    plt.xticks(range(-x_range, x_range))
    plt.yticks(range(-y_range, y_range))
    plt.grid(b=True)
    plt.plot(x, y, color='k')
    plt.scatter(x_H, y_H, color='r', marker='o')
    plt.scatter(x_P, y_P, color='b', marker='o')
    plt.scatter(x_C, y_C, color='g', marker='o')
    plt.show()