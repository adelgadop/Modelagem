from __future__ import division
from tty import CFLAG
from scipy.sparse.linalg import spsolve
import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import os

def wave(t,n):
    om = 2*np.pi/1800
    return np.sin(om*t[n])

def euler_back(C, n, CFL):
    C[1:, n] = C[1:, n-1] - CFL*(C[1:, n-1]- C[:-1, n-1])
    return C

def leap2(C, n, CFL, f_space=False, f_time = False, alfa=0.105):
    if n == 1:
    # Euler 2d order scheme for the first time step
        C = euler_back(C, n, CFL)   # rev ok

    elif n > 1:
        C[1:-1, n] = C[1:-1, n-2] - CFL*(C[2:, n-1] - C[:-2, n-1])   # rev ok
        # radiacional
        C[-1, n] = C[-1, n-1] - CFL*(C[-1, n-1]- C[-2, n-1])
    
    # Filtro Asselin desde a mitade da grade
    if f_space == True:
        C[1:50,n] = C[1:50, n-1] + alfa*(C[0:49, n] - 2*C[1:50, n-1] + C[2:51, n-1])
    
    if f_time == True:
        C[:50,2:-1] = C[:50, 1:-2] + alfa*(C[:50, 1:-2] - 2*C[:50, 1:-2] + C[:50, 2:-1])
    else:
        pass
    return C

def ordem4(C,n, CFL, f_space=False, f_time = False, alfa=0.105):
    # Aprox leapfrog 4a ordem
    # --------------
    if n == 1:
        # Euler 2d order scheme for the first time step
        C = euler_back(C, n, CFL)  # rev ok

    elif n > 1:
        C[2:-2, n] = C[2:-2,n-2] - CFL/6*(C[:-4,n-1] - 8*C[1:-3,n-1] + 8*C[3:-1,n-1] - C[4:,n-1])
        # radiacional
        C[-2, n] = C[-2, n-1] - CFL*(C[-2, n-1] - C[-3, n-1])
        C[-1, n] = C[-1, n-1] - CFL*(C[-1, n-1] - C[-2, n-1])
    
    # Filtro Asselin desde a mitade da grade
    if f_space == True:
        C[1:50,n] = C[1:50, n-1] + alfa*(C[0:49, n] - 2*C[1:50, n-1] + C[2:51, n-1])
    
    if f_time == True:
        C[:50,2:-1] = C[:50, 1:-2] + alfa*(C[:50, 1:-2] - 2*C[:50, 1:-2] + C[:50, 2:-1])
    else:
        pass
    return C

def crank_matrix(CFL, x, uc=0.5):
    import scipy.sparse as sp
    uns   = np.ones(len(x))
    r     = uns*CFL/2
    diags = (-1, 0, 1)  # -1 low diagonal, 0 main diagonal, 1 upper diagonal
    A = sp.spdiags( [-uc*r, uns, uc*r], diags, len(x), len(x) )
    A = (sp.lil_matrix(A)).tocsr()
    B = sp.spdiags( [(1-uc)*r, uns, -(1-uc)*r], diags, len(x), len(x) )
    B = (sp.lil_matrix(B)).tocsr()

    return A, B

def crank(A, B, C, n, CFL, f_space=False, f_time = False, alfa=0.105):
    C[:, n] = spsolve(A, B*C[:, n-1])
    C[-1, n] = C[-1, n-1] - CFL*(C[-1, n-1] - C[-2, n-1])

    # Filtro Asselin desde a mitade da grade
    if f_space == True:
        C[1:50,n] = C[1:50, n-1] + alfa*(C[0:49, n] - 2*C[1:50, n-1] + C[2:51, n-1])
    
    if f_time == True:
        C[:50,2:-1] = C[:50, 1:-2] + alfa*(C[:50, 1:-2] - 2*C[:50, 1:-2] + C[:50, 2:-1])
    else:
        pass
    
    return C

def hovm(X,T,C_s, ylabel,name, n, dt, t):
    fig, ax = plt.subplots(1,2,figsize=(12,4), gridspec_kw={'wspace':.15})
    levels = [0,0.2,0.4,0.6,0.8,0.9,1.]
    im = ax[0].contourf(X/1000, T/3600, C_s, levels, origin='lower', extend='both',
                        colors=['azure', 'bisque','orange', 'r', 'brown', 'darkred'])
    im.cmap.set_under('w')
    im.cmap.set_over('k')
    ax[0].set_ylabel("Tempo em horas")
    ax[0].set_xlabel("km")
    ax[0].set_title(f"Diagrama Hovmoller", loc='left')
    ax[1].set_ylabel(ylabel)
    ax[1].set_xlabel("Tempo (h)")
    ax[1].set_title(f"Ponto j=50, t[{n}] = {round((n*dt/3600),1)} horas, $\Delta t$={dt} s", loc='left')
    fig.colorbar(im, ax=ax[0],orientation="vertical",shrink=0.9) #fraction=0.04, pad=0.08
    ax[1].plot(t/3600, C_s[50,:])
    fig.savefig("fig/"+name+".png", 
            dpi = 400, bbox_inches='tight', facecolor='w')

