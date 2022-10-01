from __future__ import division
from tty import CFLAG
from scipy.sparse.linalg import spsolve
import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import os

def wave(t, n):
    om = 2*np.pi/1800
    res = np.sin(om*t[n])
    fonte = lambda res: res if res >= 0 else 0

    return fonte(res)

def euler_back(C, n, CFL, dt, F):
    C[1:, n] = C[1:, n-1] + F[1:, n-1]*dt - CFL*(C[1:, n-1]- C[:-1, n-1])
    return C

def filtro(C, alfa, f_space=False, f_time=False):
    """Filtro Robert-Asselin para remover o modo computacional

    Args:
        C (array): Campos nulos com a emissão na metade
        n (int): passo de tempo
        alfa (float, optional): Coeficiente.
        f_space (bool, optional): Filtro no espaço. Defaults to False.
        f_time (bool, optional): Filtro no tempo. Defaults to False.

    Returns:
        array: valores de C filtrados
    """
    # Filtro Asselin desde a mitade da grade
    if f_space == True:
        C[1:100,:] = C[1:100, :] + alfa*(C[:99, :] - 2*C[1:100, :] + C[2:101, :])
    
    if f_time == True:
        C[:100,1:-1] = C[:100, 1:-1] + alfa*(C[:100, :-2] - 2*C[:100, 1:-1] + C[:100, 2:])
        
    else:
        pass
    
    return C

def leap2(C, n, CFL, dt, F, f_space=False, f_time = False, alfa=0.105):
    if n == 1:
    # Euler
        C = euler_back(C, n, CFL, dt, F)   

    elif n > 1:
        C[1:-1, n] = C[1:-1, n-2] + 2*dt*F[1:-1, n-1]- CFL*(C[2:, n-1] - C[:-2, n-1])   
        # radiacional
        C[-1, n] = C[-1, n-1] + F[-1,n-1]*dt - CFL*(C[-1, n-1]- C[-2, n-1])
    
    # Filtro Robert-Asselin
    C = filtro(C, alfa, f_space, f_time)
   
    return C

def ordem4(C,n, CFL, dt, F, f_space=False, f_time = False, alfa=0.105):
    # Aprox leapfrog 4a ordem
    # --------------
    if n == 1:
        # Euler
        C = euler_back(C, n, CFL, dt, F)  

    elif n > 1:
        C[2:-2, n] = C[2:-2,n-2] + F[2:-2,n-1]*2*dt - CFL/6*(C[:-4,n-1] - 8*C[1:-3,n-1] + 8*C[3:-1,n-1] - C[4:,n-1])
        # radiacional
        C[-2, n] = C[-2, n-1] + F[-2,n-1]*dt - CFL*(C[-2, n-1] - C[-3, n-1])
        C[-1, n] = C[-1, n-1] + F[-1,n-1]*dt - CFL*(C[-1, n-1] - C[-2, n-1])
    
    # Filtro Robert-Asselin
    C = filtro(C, alfa, f_space, f_time)
    
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

    # Filtro Robert-Asselin
    C = filtro(C, alfa, f_space, f_time)
    
    return C

def hovm(X,T,C_s, ylabel,name, n, dt, t, fonte=True, levels=[-1,-.1,0,.5,1,10,50, 150, 200]):
    colores = ['c','snow','w','bisque', 'peachpuff','orange', 'r', 'brown','k']

    if fonte == True:
        fig, ax = plt.subplots(1,2,figsize=(12,4), gridspec_kw={'wspace':.15})
        im = ax[0].contourf(X/1000, T/3600, C_s, levels, origin='lower', extend='both',
                            colors=colores)
        #im.cmap.set_under('w')
        #im.cmap.set_over('k')
        ax[0].set_ylabel("Tempo em horas")
        ax[0].set_xlabel("km")
        ax[0].set_title(f"Diagrama Hovmoller", loc='left')
        ax[1].set_ylabel(ylabel)
        ax[1].set_xlabel("Tempo (h)")
        ax[1].set_title(f"Ponto j=100, t[{n}] = {round((n*dt/3600),1)} horas, $\Delta t$={dt} s", loc='left')
        fig.colorbar(im, ax=ax[0],orientation="vertical",shrink=0.9) #fraction=0.04, pad=0.08
        ax[1].plot(t/3600, C_s[100,:])
    else:
        fig, ax = plt.subplots(1,figsize=(6,4))
        im = ax.contourf(X/1000, T/3600, C_s, levels, origin='lower', extend='both',
                            colors=colores)
        #im.cmap.set_under('w')
        #im.cmap.set_over('k')
        ax.set_ylabel("Tempo em horas")
        ax.set_xlabel("km")
        ax.set_title(f"Diagrama Hovmoller", loc='left')
        fig.colorbar(im, ax=ax,orientation="vertical",shrink=0.9) #fraction=0.04, pad=0.08

    fig.savefig("fig/"+name+".png", 
                dpi = 400, bbox_inches='tight', facecolor='w')

