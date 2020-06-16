import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def visualize_protein(protein):
    
    if protein.D == 2:
        x = [amino.coordinate[0] for amino in protein.aminos]
        y = [amino.coordinate[1] for amino in protein.aminos]

        x_H = [amino.coordinate[0] for amino in protein.aminos if amino.type == 'H']
        y_H = [amino.coordinate[1] for amino in protein.aminos if amino.type == 'H']

        x_P = [amino.coordinate[0] for amino in protein.aminos if amino.type == 'P']
        y_P = [amino.coordinate[1] for amino in protein.aminos if amino.type == 'P']

        x_C = [amino.coordinate[0] for amino in protein.aminos if amino.type == 'C']
        y_C = [amino.coordinate[1] for amino in protein.aminos if amino.type == 'C']

        x_max = max(x)
        x_min = min(x)
        y_max = max(y)
        y_min = min(y)

        plt.figure()
        plt.xlim(x_min - 1, x_max + 1)
        plt.ylim(y_min - 1, y_max + 1)
        plt.xticks(range(x_min - 1, x_max + 1))
        plt.yticks(range(y_min - 1, y_max + 1))
        plt.grid(b=True)
        plt.plot(x, y, color='k')
        plt.scatter(x_H, y_H, color='r', marker='o')
        plt.scatter(x_P, y_P, color='b', marker='o')
        plt.scatter(x_C, y_C, color='g', marker='o')
        plt.show()

    elif protein.D == 3:
        x = [amino.coordinate[0] for amino in protein.aminos]
        y = [amino.coordinate[1] for amino in protein.aminos]
        z = [amino.coordinate[2] for amino in protein.aminos]

        x_H = [amino.coordinate[0] for amino in protein.aminos if amino.type == 'H']
        y_H = [amino.coordinate[1] for amino in protein.aminos if amino.type == 'H']
        z_H = [amino.coordinate[2] for amino in protein.aminos if amino.type == 'H']

        x_P = [amino.coordinate[0] for amino in protein.aminos if amino.type == 'P']
        y_P = [amino.coordinate[1] for amino in protein.aminos if amino.type == 'P']
        z_P = [amino.coordinate[2] for amino in protein.aminos if amino.type == 'P']

        x_C = [amino.coordinate[0] for amino in protein.aminos if amino.type == 'C']
        y_C = [amino.coordinate[1] for amino in protein.aminos if amino.type == 'C']
        z_C = [amino.coordinate[2] for amino in protein.aminos if amino.type == 'C']

        x_max = max(x)
        x_min = min(x)
        y_max = max(y)
        y_min = min(y)
        z_max = max(z)
        z_min = min(z)

        plt.figure()
        ax = plt.axes(projection="3d")
        ax.set_xlim(x_min - 1, x_max + 1)
        ax.set_ylim(y_min - 1, y_max + 1)
        ax.set_zlim(z_min - 1, z_max + 1)
        ax.set_xticks(range(x_min - 1, x_max + 1))
        ax.set_yticks(range(y_min - 1, y_max + 1))
        ax.set_zticks(range(z_min - 1, z_max + 1))
        ax.grid(b=True)
        ax.plot(x, y, z, color='k')
        ax.scatter(x_H, y_H, z_H, color='r', marker='o')
        ax.scatter(x_P, y_P, z_P, color='b', marker='o')
        ax.scatter(x_C, y_C, z_C, color='g', marker='o')
        plt.show()