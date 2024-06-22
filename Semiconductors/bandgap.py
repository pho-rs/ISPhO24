import matplotlib.pyplot as plt
import numpy as np
import math

filename='table.txt'

k_0 = np.loadtxt(fname=filename,usecols=(0),skiprows=1,ndmin=1)
E = np.loadtxt(fname=filename,usecols=(1),skiprows=1,ndmin=1)

fig, ax = plt.subplots(1)

if len(k_0) > 1:
    k_0, E = zip( *sorted( zip(k_0, E) ) )

ax.plot(k_0, E,'.-')
ax.set_xlabel(r'$k_0$')
ax.set_ylabel(r'$E_g$, eV')
ax.set_xlim(left=0)
ax.set_ylim(bottom=0)
plt.show()
