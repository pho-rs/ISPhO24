import matplotlib.pyplot as plt
import numpy as np
import math
import cmath
import tqdm

fig, axs = plt.subplots(1)


Phi = 60

n_min = 1.0
n_max = 2.0
n_n = 300

_n = np.linspace(n_min, n_max, n_n)

_res = np.ndarray(shape = (n_n), dtype=float)

n1 = 1.0

for i in range(len(_n)):
    n2 = _n[i]

    th1 = Phi/180.0*math.pi
    th2 = cmath.asin(cmath.sin(th1)*n1/n2)

    rs = (n1*cmath.cos(th1) - n2 *cmath.cos(th2))/(n1*cmath.cos(th1) + n2 *cmath.cos(th2))
    rp = (n1*cmath.cos(th2) - n2 *cmath.cos(th1))/(n1*cmath.cos(th2) + n2 *cmath.cos(th1))

    res = (abs(rp) / abs(rs))**2

    _res[i]=res

c = axs.plot(_n, _res)
axs.set_xlabel(r'$n$')
axs.set_ylabel(r'$R_p/R_s$')

plt.show()

