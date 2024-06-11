import matplotlib.pyplot as plt
import numpy as np
import cmath
import math

a_x = -0.1
q_x = 0.8

dxi=0.01
ximax = 150.0

x=0.0
dx=4.0

fig, axs = plt.subplots(2)

xi=0.0

_x = []
_dx = []
_xi = []

while xi < ximax:
    ddx = -(a_x-2*q_x*math.cos(2*xi))*x

    x += dx * dxi
    dx += ddx * dxi
    
    _x.append(x)
    _dx.append(dx)
    _xi.append(xi)
    xi += dxi


axs[0].plot(_xi,_x,label=r'$a_x$={:.2f}, $q_x$={:.2f}'.format(a_x,q_x))
axs[0].set_xlabel(r'$\xi$')
axs[0].set_ylabel(r'$x$')
axs[1].plot(_x,_dx,label=r'$a_x$={:.2f}, $q_x$={:.2f}'.format(a_x,q_x))
axs[1].set_xlabel(r'$x$')
axs[1].set_ylabel(r'd$x$/d$\xi$')

axs[0].legend()
axs[1].legend()

plt.tight_layout()
plt.show()

