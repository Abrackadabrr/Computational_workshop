import numpy as np
import matplotlib.pyplot as plt

err = np.array([0.0300383,0.0272742, 0.0262324])
h = np.array([4.83, 2.35, 1.17])

res = np.array([[0.589887, 0.0249499],
[0.700703, 0.0262324],
[0.955483, 0.0251233],
[1.15115, 0.0245842],
[1.50582, 0.0256542],
[1.88466, 0.0262035],
[2.80281, 0.0300383]])

err = res[:, 1]
h = res[:, 0]

fig, ax = plt.subplots()

ax.plot(h, err)
ax.grid()
plt.show()

