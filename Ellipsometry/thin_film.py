import matplotlib.pyplot as plt
import numpy as np
import math
import cmath
import tqdm

fig, axs = plt.subplots(1)

I = complex(0,1)


d_min = 0.0
d_max = 200.0
d_n = 500

wl = 500.0

_d = np.linspace(d_min, d_max, d_n)


_res = np.ndarray(shape = (d_n), dtype=float)

n1 = 1.65
n2 = complex(4.32,-0.073)

for i in range(len(_d)):
    d = _d[i]
    
    r01 = (n1 - 1.0) / (1 + n1)
    r12 = (n2 - n1) / (n2 + n1)

    phi = n1 * d * 2 * math.pi / wl

    rs = (r12 + r01 * cmath.exp(2*I*phi))/(r01*r12 + cmath.exp(2*I*phi))
    
    res = (abs(rs))**2

    _res[i]=res

c = axs.plot(_d, _res)

axs.set_xlabel(r'$d, $nm')
axs.set_ylabel(r'$R$')

plt.tight_layout()
plt.show()

