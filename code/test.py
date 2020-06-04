import matplotlib.pyplot as plt

from protein import Protein

protein = Protein(['H', 'H', 'P', 'H', 'H', 'H', 'P', 'H'])

x = [amino.coordinate[0] for amino in protein.aminos]
y = [amino.coordinate[1] for amino in protein.aminos]

n = protein.length

plt.figure()
plt.xlim(-5, 5)
plt.ylim(-5, 5)
plt.xticks(range(-5, 5))
plt.yticks(range(-5, 5))
plt.grid(b=True)
plt.plot(x, y)
plt.show()