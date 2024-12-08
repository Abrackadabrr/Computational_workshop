import numpy as np
import matplotlib.pyplot as plt

res_2d = np.array([[-17.1758626309616,-23.17576867647152,152.641756485608,],
                      [-15.33031317524716,-18.93769785522725,152.080306635044,],
                      [-15.42825788118584,-19.13692295433297,152.1995387398545,],
                      [-15.12362944936843,-19.77725473267999,152.3071700729725,],
                      [-15.26124454584735,-19.80009113992879,152.2881310841184,]])
order_2d = np.array([1, 2, 3, 4, 5])
res_1d = np.array([])
order_1d = np.array([1, 2, 3, 5, 7, 9])
res_3d = np.array([7225.917709660245,7222.691681864098,7222.668948560624,7222.6689145325])
order_3d = np.array([1, 2, 3])

res = res_2d
order = order_2d
diff = list(np.abs(res - res[-1]))
diff = [np.linalg.norm(i) for i in diff]
print(diff)

fig, ax = plt.subplots()
ax.plot(order, diff)
ax.grid()
ax.set_title('Объемный интеграл')
ax.set_xlabel('Порядок интегрирования')
ax.set_ylabel('Ошибка результата')
ax.set_yscale('log')

plt.show()