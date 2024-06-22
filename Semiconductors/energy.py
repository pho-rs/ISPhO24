import matplotlib.pyplot as plt
import numpy as np
import cmath
import math

I = 1j

def lambda_1(alpha, beta):
    return -I*alpha + cmath.sqrt(-beta)
def lambda_2(alpha, beta):
    return -I*alpha - cmath.sqrt(-beta)

def calc_err(k, k_0, E, b):
    l_21 = lambda_1(k, k_0**2 * (E+1/b))
    l_22 = lambda_2(k, k_0**2 * (E+1/b))

    l_11 = lambda_1(k, k_0**2 * E)
    l_12 = lambda_2(k, k_0**2 * E)

    M = np.array([[1.0,-1.0,1.0,-1.0],
                [l_11,-l_21,l_12,-l_22],
                [cmath.exp(-l_11*b),-cmath.exp(l_21*(1-b)),cmath.exp(-l_12*b),-cmath.exp(l_22*(1-b))],
                [l_11*cmath.exp(-l_11*b),-l_21*cmath.exp(l_21*(1-b)),l_12*cmath.exp(-l_12*b),-l_22*cmath.exp(l_22*(1-b))]])

    M1 = np.array([[1.0,-1.0,1.0,],
                [l_11,-l_21,l_12,],
                [cmath.exp(-l_11*b),-cmath.exp(l_21*(1-b)),cmath.exp(-l_12*b)]])

    err = abs(np.linalg.det(M))/abs(np.linalg.det(M1))

    return float(err)


def main():
    fig, axs = plt.subplots(1)

    k_0 = 4.0
    b = 0.1
    YSCALE=2.8

    E_size = 10000
    k_size=25

    E_F = 4/k_0**2*math.pi**2

    __E = np.linspace(-1.0/b, -1.0/b+E_F*YSCALE, int(E_size))
    __k = np.linspace(0,math.pi, k_size)
    
    k_res = []
    E_res = []

    k_err = []
    E_err = []
    err_err = []

    i=0
    k=math.pi

    for k in __k:
        print(f'{i}/{k_size}',end='\r')
        i+=1

        _err = []
        _E = []
        for E in __E:
            err = calc_err(k, k_0, E, b)

            _err.append(err)
            _E.append(E)

            k_err.append(k)
            E_err.append(E)
            err_err.append(math.log(err+1e-25))

            if len(_err) < 2:
                continue

            n = len(_err)
            if _err[n-2] < _err[n-3] and _err[n-2] < _err[n-1]:

                if _E[n-3] > 0 and _E[n-1]<0:
                    continue

                E1 = _E[n-3]
                E2 = _E[n-1]

                err1 = _err[n-3]
                err2 = _err[n-1]

                n_err = err1
                n_E = E1
                dE = E2 - E1

                mindE = 1.0e-10
                min_err = 0.001
                while n_err > min_err and dE>0.0:
                    
                    n_E = E1 + dE/2.0
                    n_err = calc_err(k,k_0, n_E,b)

                    if n_err > err1 or n_err > err2:
                        break
                    if dE < mindE and n_err > 1.0:
                        break

                    derr1 = err1 - n_err
                    derr2 = err2 - n_err

                    if derr1 > derr2:
                        E1 += dE * err1/(err1 + err2)/2.0
                        err1 = calc_err(k,k_0, E1,b)
                        
                    else:
                        E2 -= dE * err2/(err1 + err2)/2.0
                        err2 = calc_err(k,k_0, E2,b)

                    dE = E2 - E1

                if n_err <= min_err:
                    k_res.append(k)
                    E_res.append(n_E)
                    
    
    axs.plot(k_res, E_res,'.',label=r'$k_0={:.2f}, b/a={:.2f}$'.format(k_0, b))

    E_gr = np.min(E_res)

    E_F2 = 1/k_0**2*math.pi**2
    E_F6 = 9/k_0**2*math.pi**2

    k_con = []
    E_con = []

    k_val = []
    E_val = []

    E_val2 = []
    k_val2 = []

    E_con2 = []
    k_con2 = []

    for i in range(len(k_res)):
        if E_res[i] <= E_gr + E_F and E_res[i] >= E_gr + E_F2 and k_res[i] < math.pi/4:
            E_val.append(E_res[i])
            k_val.append(k_res[i])

        if E_res[i] >= E_gr + E_F and E_res[i] <= E_gr + E_F6 and k_res[i] < math.pi/4:
            E_con.append(E_res[i])
            k_con.append(k_res[i])

        if E_res[i] <= E_gr + E_F2:
            E_val2.append(E_res[i])
            k_val2.append(k_res[i])

        if E_res[i] >= E_gr + E_F6:
            E_con2.append(E_res[i])
            k_con2.append(k_res[i])


    k_con = np.array(k_con)
    E_con = np.array(E_con)
    k_val = np.array(k_val)
    E_val = np.array(E_val)


    k_con, E_con = zip( *sorted( zip(k_con, E_con) ) )
    k_val, E_val = zip( *sorted( zip(k_val, E_val) ) )

    con = np.polyfit(np.power(k_con,2), E_con,1)
    val = np.polyfit(np.power(k_val,2), E_val,1)

    print(f'e_c={con[1]}\t~alpha={con[0]}')
    print(f'e_v={val[1]}\t~beta={-val[0]}')

    #print(f'{(E_con[1]-E_val[1])*k_0**2*0.151} eV')

    #axs.fill_between([0, math.pi], [np.max(E_val2),np.max(E_val2)],[np.min(E_val), np.min(E_val)],color='red',alpha=0.3)
    #axs.fill_between([0, math.pi], [np.max(E_val),np.max(E_val)],[np.min(E_con), np.min(E_con)],color='red',alpha=0.3)
    #axs.fill_between([0, math.pi], [np.max(E_con),np.max(E_con)],[np.min(E_con2), np.min(E_con2)],color='red',alpha=0.3)

    #axs.plot(k_val2, E_val2,'.-', color='tab:blue',label=r'$k_0={:.2f}, b/a={:.2f}$'.format(k_0, b))
    #axs.plot(k_val, E_val,'.-', color='tab:blue')
    #axs.plot(k_con, E_con,'.-', color='tab:blue')
    #axs.plot(k_con2, E_con2,'.-', color='tab:blue')

    axs.plot([0,cmath.pi],[E_gr+E_F, E_gr+E_F],label='Estimated Fermi level')

    axs.plot(k_con,np.power(k_con,2)*con[0]+con[1],label=r'$\epsilon_c+\tilde{\alpha}k^2$')
    axs.plot(k_val,np.power(k_val,2)*val[0]+val[1],label=r'$\epsilon_v-\tilde{\beta}k^2$')

    axs.set_ylabel('\u03B5')
    axs.set_xlabel('k')

    axs.legend()
    
    
    plt.show()

    return

    fig, axs = plt.subplots(1)

    k_0 = 10
    b = 0.1

    E_size = 10000
    k_size = 50

    __E = np.linspace(-1.01,1.0, E_size) 
    __k = np.linspace(0,math.pi, k_size)
    _err = []
    _E = []
    _k = []

    i=0
    for k in __k:
        print(f'{i}/{k_size}')
        i+=1
        for E in __E:
            l_11 = lambda_1(k, k**2 - k_0**2 * E)
            l_12 = lambda_2(k, k**2 - k_0**2 * E)

            l_21 = lambda_1(k, k**2 + k_0**2 * (-1.0-E))
            l_22 = lambda_2(k, k**2 + k_0**2 * (-1.0-E))

            M = np.array([[1.0,-cmath.exp(l_21),1.0,-cmath.exp(-l_22)],
                        [l_11,-l_21*cmath.exp(l_21),l_12,-l_22*cmath.exp(l_22)],
                        [cmath.exp(l_11*b),-cmath.exp(l_21*b),cmath.exp(l_12*b),-cmath.exp(l_22*b)],
                        [l_11*cmath.exp(l_11*b),-l_21*cmath.exp(l_21*b),l_12*cmath.exp(l_12*b),-l_22*cmath.exp(l_22*b)]])
            
            err = abs(np.linalg.det(M))

            if err < 1e-2:
                print(f'k={k}\tE={E}\nl_11={l_11}\tl_12={l_12}\nl_21={l_21}\tl_22={l_22}\n')

            _err.append(math.log(err+1e-25))
            _E.append(E)
            _k.append(k)

    _err = np.array(_err).reshape((k_size, E_size))
    _E = np.array(_E).reshape((k_size, E_size))
    _k = np.array(_k).reshape((k_size, E_size))
    
    c0=axs.pcolor(_k,_E, _err)
    fig.colorbar(c0, ax=axs)

    plt.show()


    return

    j=0
    for q_x in __q_x:
        print(f'{j}/{_q_x_size}')
        for a_x in __a_x:
            x=0.01
            dx=0.00
            t=0.0
            x_max = x
            x_half_max = x
        
            while t < tmax:
                ddx = -(a_x-2*q_x*math.cos(2*t))*x

                dx += ddx * dt
                x += dx * dt

                if abs(x) > x_max:
                    x_max = abs(x)

                    if t < tmax/2:
                        x_half_max = x_max

                t += dt
            
            if x_max/x_half_max < 1.1:
                _st.append(1.0)
            else:
                _st.append(0.0)

            _x_max.append(math.log(x_max/x_half_max))
            _a_x.append(a_x)
            _q_x.append(q_x)

        j+=1

    with open('out.txt', 'w') as f_out:
        for i in range(len(_x_max)):
            print(f'{_a_x[i]}\t{_q_x[i]}\t{_x_max[i]}\t{_st[i]}',file=f_out)


    _x_max = np.array(_x_max).reshape((_q_x_size, _a_x_size))
    _a_x = np.array(_a_x).reshape((_q_x_size, _a_x_size))
    _q_x = np.array(_q_x).reshape((_q_x_size, _a_x_size))
    _st = np.array(_st).reshape((_q_x_size, _a_x_size))

    c0=axs[0].pcolor(_q_x,_a_x, _x_max)
    fig.colorbar(c0, ax=axs[0])
    c1=axs[1].pcolor(_q_x,_a_x, _st)
    fig.colorbar(c1, ax=axs[1])

    axs[0].set_xlabel('q_x')
    axs[0].set_ylabel('a_x')
    axs[1].set_xlabel('q_x')
    axs[1].set_ylabel('a_x')

    plt.show()


main()