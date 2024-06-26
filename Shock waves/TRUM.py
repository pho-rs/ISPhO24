import matplotlib.pyplot as plt
import numpy as np
import math as m


M_values = [1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
theta = 10
# make data
for M in M_values:
    b_values = []
    t_values = []
    for i in range(1, 10000):
        b = m.pi * i / 20000

        t = m.atan(2.0 * m.tan(b)**(-1) * (M**2 * m.sin(b)**2 -1)/(M**2 * (1.4 + m.cos(2*b)) + 2))
        if t > 0:
            b_values.append(b * 360 / 2 / m.pi)
            t_values.append(t * 360 / 2 / m.pi)

    plt.plot(t_values, b_values, label=f'M={M}')
    plt.legend()

x = [theta] * len(b_values)
plt.plot(x, b_values, color='silver')
plt.xlabel(r'$\theta, ^\circ$')
plt.ylabel(r'$\beta, ^\circ$')
plt.tight_layout()
plt.show()