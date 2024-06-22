import matplotlib.pyplot as plt
import numpy as np
import cmath
import math
import scipy

xmin=0.0
xmax=1.0
dx_val = 1e-6

def Psi_sim(Np, Nn, kA=0.001, region='p'):
    x = xmin

    psi = 0.0
    dpsi= 0.0
    Ndop = 0.0

    dx = dx_val
    if region == 'p':
        psi = math.log(Np) - math.log(Np) * kA
        dpsi = - math.log(Np) * kA * math.sqrt(Np)
        Ndop = -Np
    elif region == 'n':
        psi = -math.log(Nn) + math.log(Nn) * kA
        dpsi = -math.log(Nn) * kA * math.sqrt(Nn)
        Ndop = Nn
        dx = dx_val * -1.0

    _x = []
    _psi = []
    _dpsi = []

    i=0
    while abs(x) < xmax and psi > -math.log(Nn) and psi < math.log(Np):
        if i % 20 == 0:
            _psi.append(psi)
            _dpsi.append(dpsi)
            _x.append(x)
        i+=1

        ddpsi = 2 * math.sinh(psi) + Ndop

        psi += dpsi * dx + ddpsi * dx**2 / 2
        dpsi += ddpsi * dx
        x += dx


    return np.array(_x), np.array(_psi), np.array(_dpsi)

def cross(_psi1, _dpsi1, _x1, _psi2, _dpsi2, _x2):
    _dpsi2_ = np.interp(_psi1.copy(), _psi2, _dpsi2)
    cross_dpsi = np.interp(0.0,_dpsi2_ - _dpsi1,_dpsi1)
    cross_psi = np.interp(0.0,_dpsi2_ - _dpsi1, _psi1)
    cross_x1 = np.interp(-cross_psi,-_psi1, _x1)
    cross_x2 = np.interp(cross_psi,_psi2, _x2)

    return cross_psi, cross_dpsi, cross_x1, cross_x2

fig, axs = plt.subplots()

def main(Np, Nn, J):
    _x1, _psi1, _dpsi1 = Psi_sim(Np, Nn, region = 'p')
    _x2, _psi2, _dpsi2 = Psi_sim(Np, Nn, region = 'n')
    
#    axs[0,0].plot(_x1,_psi1,label='p')
#    axs[0,0].plot(_x2,_psi2,label='n')
#    axs[0,0].legend()

    cross_psi, cross_dpsi, cross_x1, cross_x2 = cross(_psi1,_dpsi1,_x1,_psi2,_dpsi2,_x2)

 #   axs[0,0].scatter(cross_x1,cross_psi)
 #   axs[0,0].scatter(cross_x2,cross_psi)

#    axs[0,1].plot(_psi1,_dpsi1,label='p')
#    axs[0,1].plot(_psi2,_dpsi2,label='n')
#    axs[0,1].scatter(cross_psi,cross_dpsi)
#    axs[0,1].legend()

    x1 = []
    psi1 = []
    dpsi1 = []
    x2 = []
    psi2 = []
    dpsi2 = []

    for i in range(len(_x1)):
        if _x1[i]<cross_x1:
            x1.append(_x1[i]-cross_x1)
            psi1.append(_psi1[i])
            dpsi1.append(_dpsi1[i])
    for i in range(len(_x2)-1,0,-1):
        if _x2[i]>cross_x2:
            x2.append(_x2[i]-cross_x2)
            psi2.append(_psi2[i])
            dpsi2.append(_dpsi2[i])

    x = np.array(x1+x2)
    psi = np.array(psi1+psi2)
    dpsi = np.array(dpsi1+dpsi2)

    axs.plot(x,psi)
    axs.set_xlabel(r'$x/L$')
    axs.set_ylabel(r'$\Psi$')

Np=2.0e5
Nn=1.0e5
J=0.0
main(Np, Nn, J)

plt.tight_layout()
plt.show()
