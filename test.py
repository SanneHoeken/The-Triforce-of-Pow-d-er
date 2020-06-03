import matplotlib.pyplot as plt

from protein import Protein

protein = Protein(['H', 'H', 'P', 'H', 'H', 'H', 'P', 'H', 'P', 'H', 'H', 'H', 'P', 'H'])

x = [coordinate[0] for coordinate in protein.coordinates]
y = [coordinate[1] for coordinate in protein.coordinates]
n = protein.length

plt.figure()
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.xticks(range(-5, 5))
plt.yticks(range(-5, 5))
plt.grid(b=True)
plt.plot(x, y)
plt.show()