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

def euler_back(C, Nt, F, t, dt, CFL, Nx):
    for n in range(1, Nt-1):
        F[101, n] = 1/2*(wave(t, n)+wave(t, n-1))
        # advecção ordem 1

        for i in range(1, Nx-1):
            C[i, n+1] = C[i, n] + F[i, n]*dt - CFL*(C[i, n]- C[i-1, n])
        # radiacional
        C[-1, n+1] = C[-1, n] - CFL*(C[-1, n]- C[-2, n])
    return C

def filter(C, gamma, alfa, n, Nx, f_space, f_time):
    """Filtro Robert-Asselin-Williams para remover o modo computacional
    gamma = 0.1 ou 0.01
    alfa = 0.53
    """
    for i in range(1, Nx-1):
        if f_time == True:
            C[i,  n] = C[i, n  ] +     gamma*alfa/2*(C[i, n+1] - 2*C[i, n] + C[i, n-1])
            C[i,n+1] = C[i, n+1] - gamma*(1-alfa)/2*(C[i, n+1] - 2*C[i, n] + C[i, n-1])
        
        elif f_space == True:
            C[i ,  n] = C[i,   n] +     gamma*alfa/2*(C[i+1, n] - 2*C[i, n] + C[i-1, n])
            C[i+1, n] = C[i+1, n] - gamma*(1-alfa)/2*(C[i+1, n] - 2*C[i, n] + C[i-1, n])
    return C

def leap2(C, F, t, dt, CFL, Nt, Nx, gamma, alfa, f_space, f_time):
    for n in range(1, Nt-1):
        F[101, n] = 1/2*(wave(t, n)+wave(t, n-1))
               
        if n == 1: # Euler
            for i in range(1, Nx-1):
                C[i, n+1] = C[i, n] + F[i, n]*dt - CFL*(C[i, n]- C[i-1, n])
                
        if n > 1: # leapfrog
            for i in range(1, Nx-1):
                C[i, n+1] = C[i, n-1] + F[i, n]*2*dt - CFL*(C[i+1, n]- C[i-1, n])
            
        # radiacional
        C[-1, n+1] = C[-1, n] - CFL*(C[-1, n]- C[-2, n])  
            
        # Robert-Asselin-Williams filter
        for i in range(1, Nx-1):
            if f_time == True:
                C[i,  n] = C[i, n  ] +     gamma*alfa/2*(C[i, n+1] - 2*C[i, n] + C[i, n-1])
                C[i,n+1] = C[i, n+1] - gamma*(1-alfa)/2*(C[i, n+1] - 2*C[i, n] + C[i, n-1])
            
            if f_space == True:
                C[i ,  n] = C[i,   n] +     gamma*alfa/2*(C[i+1, n] - 2*C[i, n] + C[i-1, n])
                C[i+1, n] = C[i+1, n] - gamma*(1-alfa)/2*(C[i+1, n] - 2*C[i, n] + C[i-1, n])

    return C

def ordem4(C,  F, t, dt, CFL, Nt, Nx,  gamma, alfa, f_space, f_time):
    for n in range(1, Nt-1):
        # Aprox leapfrog 4a ordem
        # --------------
        F[101, n] = 1/2*(wave(t, n)+wave(t, n-1))
        if n == 1:
            # Euler
            C = euler_back(C, Nt, F, t, dt,CFL, Nx)  

        elif n > 1:
            
            C[2:-2, n+1] = C[2:-2,n-1] + F[2:-2,n]*2*dt - CFL/6*(C[:-4,n] - 8*C[1:-3,n] + 8*C[3:-1,n] - C[4:,n])
                                
            # radiacional
            C[-2, n+1] = C[-2, n] + F[-2,n]*dt - CFL*(C[-2, n] - C[-3, n])
            C[-1, n+1] = C[-1, n] + F[-1,n]*dt - CFL*(C[-1, n] - C[-2, n])
            
        # Robert-Asselin-Williams filter
        for i in range(1, Nx-1):
            if f_time == True:
                C[i,  n] = C[i, n  ] +     gamma*alfa/2*(C[i, n+1] - 2*C[i, n] + C[i, n-1])
                C[i,n+1] = C[i, n+1] - gamma*(1-alfa)/2*(C[i, n+1] - 2*C[i, n] + C[i, n-1])
            
            if f_space == True:
                C[i ,  n] = C[i,   n] +     gamma*alfa/2*(C[i+1, n] - 2*C[i, n] + C[i-1, n])
                C[i+1, n] = C[i+1, n] - gamma*(1-alfa)/2*(C[i+1, n] - 2*C[i, n] + C[i-1, n])
      
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

def crank(A, B, C, F, t, dt, n, Nx, CFL, gamma, alfa, f_space=False, f_time = False):
    F[101, n] =  1/2*(wave(t, n) + wave(t, n-1))
    C[:, n+1] = spsolve(A, B*C[:, n]) + F[:,n]*dt
    C[-1, n+1] = C[-1, n] + F[-1,n]*dt - CFL*(C[-1, n] - C[-2, n])
    
    # Robert-Asselin-Williams filter
    for i in range(1, Nx-1):
        if f_time == True:
            C[i,   n] = C[i, n  ] +     gamma*alfa/2*(C[i, n+1] - 2*C[i, n] + C[i, n-1])
            C[i, n+1] = C[i, n+1] - gamma*(1-alfa)/2*(C[i, n+1] - 2*C[i, n] + C[i, n-1])
            
        if f_space == True:
            C[i ,  n] = C[i,   n] +     gamma*alfa/2*(C[i+1, n] - 2*C[i, n] + C[i-1, n])
            C[i+1, n] = C[i+1, n] - gamma*(1-alfa)/2*(C[i+1, n] - 2*C[i, n] + C[i-1, n])
    return C

def hovm(X,T,C_s, F, ylabel,name, n, dt, t, alfa, gamma, filtro=False, fonte=True, levels=[-1,-.1,0,.5,1,10,50, 150, 200]):
    colores = ['c','snow','w','bisque', 'peachpuff','orange', 'r', 'brown','k']

    if fonte == True:
        fig, ax = plt.subplots(1,2,figsize=(12,4), gridspec_kw={'wspace':.15})
        im = ax[0].contourf(X/1000, T/3600, C_s, levels, origin='lower', # extend='both',
                            colors=colores)
        #im.cmap.set_under('w')
        #im.cmap.set_over('k')
        ax[0].set_ylabel("Tempo em horas")
        ax[0].set_xlabel("km")
        if filtro == True:
            ax.set_title(f"Diagrama Hovmoller, filtro (gamma={gamma}, alfa={alfa})", loc='left')
        ax[0].set_title(f"Diagrama Hovmoller", loc='left')
        ax[1].set_ylabel(ylabel[1])
        ax[1].set_xlabel("Tempo (h)")
        ax[1].set_title(f"Ponto j=101, t[{n}] = {round((n*dt/3600),1)} horas, $\Delta t$={dt} s", loc='left')
        cbar = fig.colorbar(im, ax=ax[0],orientation="vertical") #fraction=0.04, pad=0.08 ,shrink=0.8
        cbar.ax.set_title(ylabel[0], fontsize=8)
        ax[1].plot(t/3600, F[101,:])
    else:
        fig, ax = plt.subplots(1,figsize=(6,4))
        im = ax.contourf(X/1000, T/3600, C_s, levels, origin='lower', #extend='both',
                            colors=colores)
        #im.cmap.set_under('w')
        #im.cmap.set_over('k')
        ax.set_ylabel("Tempo em horas")
        ax.set_xlabel("km")
        if filtro == True:
            ax.set_title(f"Diagrama Hovmoller, \n filtro (gamma={gamma}, alfa={alfa})", loc='left')
        else:
            ax.set_title(f"Diagrama Hovmoller", loc='left')
        cbar = fig.colorbar(im, ax=ax,orientation="vertical") #fraction=0.04, pad=0.08
        cbar.ax.set_title(ylabel[0], fontsize=8)

    fig.savefig("fig/"+name+".png", 
                dpi = 400, bbox_inches='tight', facecolor='w')

