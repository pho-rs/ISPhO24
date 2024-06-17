import matplotlib.pyplot as plt
import numpy as np
import math
import cmath
import tqdm

fig, axs = plt.subplots(1)


eps1 = 1.0
eps2 = 4.0 

n1 = cmath.sqrt(eps1)
n2 = cmath.sqrt(eps2)

theta_min = 0.0
theta_max = 90.0
n_theta = 100

_theta = np.linspace(theta_min, theta_max, n_theta)

_res = np.ndarray(shape = (n_theta), dtype=float)

for i in range(len(_theta)):
    th1 = _theta[i]/180*math.pi
    th2 = cmath.asin(cmath.sin(th1)*n1/n2)

    rs = (n1*cmath.cos(th1) - n2 *cmath.cos(th2))/(n1*cmath.cos(th1) + n2 *cmath.cos(th2))
    rp = (n1*cmath.cos(th2) - n2 *cmath.cos(th1))/(n1*cmath.cos(th2) + n2 *cmath.cos(th1))

    res = (abs(rp) / abs(rs))**2

    _res[i]=res

c = axs.plot(_theta, _res)
axs.set_xlabel(r'$\theta, ^\circ$')
axs.set_ylabel(r'$R_p/R_s$')

plt.show()

