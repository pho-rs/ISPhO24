import matplotlib.pyplot as plt
import numpy as np
import math
import tqdm

fig, axs = plt.subplots(1)

def single_sim(a_x, q_x, dxi=0.1, ximax = 150):
    x=0.0
    dx=1e-1

    xi=0.0

    _x = []
    _dx = []
    _xi = []

    while xi < ximax:
        ddx = -(a_x-2*q_x*math.cos(2*xi))*x

        dx += ddx * dxi
        x += dx * dxi

        _x.append(x)
        _dx.append(dx)
        _xi.append(xi)
        xi += dxi
            
    return [np.array(_xi), np.array(_x)]

a_min = -3.0
a_max = 7.0
a_iterations = 500

q_min = 0
q_max = 5.0
q_iterations = 500

_a = np.linspace(a_min, a_max, a_iterations)
_q = np.linspace(q_min, q_max, q_iterations)

data_x = np.ndarray(shape = (a_iterations, q_iterations), dtype=float)

for i in tqdm.tqdm(range(a_iterations)):
    for j in range(q_iterations):
        ax = _a[i]
        qx = _q[j]
        
        sim = single_sim(ax, qx)
        
        _x = np.abs(sim[1])

        
        x1 = np.max(_x[:_x.shape[0] // 2])
        x2 = np.max(_x)
        
       
        data_x[i][j] = x2 / x1

bool_data_x = np.ndarray(shape = (a_iterations, q_iterations), dtype = int)

for i in range(a_iterations):
    for j in range(q_iterations):
        if data_x[i][j] < 1.1:
            bool_data_x[i][j] = 1
        else:
            bool_data_x[i][j] = 0

c = axs.pcolormesh(_q, _a, bool_data_x)
axs.set_xlabel(r'$q_x$')
axs.set_ylabel(r'$a_x$')
plt.show()

