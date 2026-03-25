import numpy as np
import matplotlib.pyplot as plt

# Создаем сетку
x = np.linspace(-10, 10, 30)
y = np.linspace(-10, 10, 30)
X, Y = np.meshgrid(x, y)

# Создаем график
fig, ax = plt.subplots(figsize=(10, 8))
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Фазовый портрет системы')
ax.grid()

# Задаем ЛОСДУ
dx = X + Y
dy = -7*X + 3*Y

# Строим поле направлений
ax.streamplot(X, Y, dx, dy, color='blue', density=1.5, arrowsize=1.5)

# Отображаем график
plt.show()