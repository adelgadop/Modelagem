from __future__ import division
from scipy.sparse.linalg import spsolve
import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import os

def Conc(x, nr):
    """
    Função de concentração com variação Gaussiana
    ---------------------------------------------
    i:      index
    nr:     número de pontos a partir do centro da gaussiana
    """
    return 10*np.exp(-(x - 51*5000)**2/(nr*5000)**2)

def Rect(x):
    """
    Rectangular function
    --------------------
    Centered in 50, 51, 52
    """
    return np.where((x >= (50*5_000)) & (x <= (52*5_000)), 5, 0 )

def crank_matrix(x, CFL, uc):
    """
    Return matrices A and B for advection equations
    -----------------------------------------------
    x      : linear space with Nx points
    r      : CFL/2
    uc     : is a Crank-Nicolson parameter equals to 1/2
    """
    import scipy.sparse as sp
    uns   = np.ones(len(x))
    r     = uns*CFL/2
    diags = (-1, 0, 1)  # -1 low diagonal, 0 main diagonal, 1 upper diagonal
    A = sp.spdiags( [-uc*r, uns, uc*r], diags, len(x), len(x) )
    A = sp.lil_matrix(A)
    B = sp.spdiags( [(1-uc)*r, uns, -(1-uc)*r], diags, len(x), len(x) )
    B = sp.lil_matrix(B)
    
    # periodic
    A[0,-1] = -CFL/4
    A[-1,0] = CFL/4
    B[0,-1] = CFL/4
    B[-1,0] = -CFL/4
    return A.tocsr(), B.tocsr()  # A.toarray()

def deriv(c,U,dx):
    """
    Derivada RK4  Rev: ok
    """
    dc = np.zeros(len(c))
    p = U/(2*dx)
    dc[2:-2] = -(U/(12*dx))*(8*(c[3:-1]-c[1:-3]) - (c[4:] - c[:-4]))
    dc[1]  = -p*(c[2] - c[-1])
    dc[-2] = -p*(c[-1]- c[-3])
    dc[0]  = -p*(c[1] - c[-2])
    dc[-1] = -p*(c[1] - c[-2])
    return dc

def sol_analytical(fun, x, U, Nx, Nt, nr, dt, t):
    RE = np.zeros((Nx, Nt))   # Matriz
    # Space: the last is Nx-1 (e.g., not use Nx)
    # Time: the last in Nt-1 (e.g., not use n+1)

    # Condições iniciais
    # -------------------
    if fun == "Gaussiana":
        RE[:,0] = Conc(x, nr)

    elif fun == "Retângulo":
        RE[:,0] = Rect(x)

    if U > 0:
        for n in range(1, Nt):
            if dt == 500:
                RE[1:, n] = RE[:-1, n-1]  # sol analitica U*dt/dx = 1
                RE[0 , n] = RE[ -1, n-1]
            else:
                if fun == "Gaussiana":
                    ab = 10*np.exp(-(x - U*t[n] - (Nx-50)*5000)**2/(nr*5000)**2)
                    ba = 10*np.exp(-(x - U*t[n] + (Nx-51)*5000)**2/(nr*5000)**2)
                    RE[:,n] = np.maximum(ab, ba)
                    #RE[0, n] = RE[-1, n-1]

                elif fun == "Retângulo":
                    ab = np.where((x >= (50*5_000 + U*t[n])) & (x <= (52*5_000 + U*t[n])), 5, 0 )
                    ba = np.where((x >= (-50*5_000 + U*t[n])) & (x <= (-48*5_000 + U*t[n])), 5, 0 )
                    RE[:,n] = np.maximum(ab,ba)
                    #RE[0, n] = RE[-1, n-1]
                
    else:
        print("precisa de codigo")
   
    return RE
    
def sol_num(aprox, cond_front, fun, CFL, nr, Nx, Nt, x, dx, U, t, dt):
    """
    aprox       : Tipo de aproximação: "ordem1", "leapfrog", "ordem4", "Matsuno", "implicito", "RK4" (Runge Kutta 4).
    cond_front  : Condição de fronteira, tipo "fixa", "periódica" e "radiacional".
    fun         : Função tipo "Gaussiana" e "Retângulo"
    hora        : Tempo em horas (int)
    """
    C = np.zeros((Nx, Nt))   # Matriz
    # Space: the last is Nx-1 (e.g., not use Nx)
    # Time: the last in Nt-1 (e.g., not use n+1)
    u  = C.copy()
    
    A, B = crank_matrix(x, CFL, 0.5)  # 0.5 because 1 + n/2

    # Condições iniciais
    # -------------------
    if fun == "Gaussiana":
        C[:,0] = Conc(x, nr)
        u[:,0] = Conc(x, nr)
        c = C[:,0].copy()

    elif fun == "Retângulo":
        C[:,0] = Rect(x)
        u[:,0] = Rect(x)
        c = C[:,0].copy()
        
    for n in range(1, Nt):

        # Aproximação numérica:
        # ---------------------
        if aprox == "ordem1": # ok
            C[1:, n] = C[1:, n-1] - CFL*(C[1:, n-1]- C[:-1, n-1])

            # Condição de fronteira
            # ---------------------
            if cond_front == 'fixa':
                C[[0,-1], n] = [0,0]

            elif cond_front == 'periódica':
                C[0, n]  = C[-1, n-1]
                
            elif cond_front == 'radiacional':
                C[0, n] = C[ -1, n-1] - CFL*(C[-1, n-1]- C[-2, n-1])

        # -------------------------------------------------
        elif aprox == "leapfrog":  # ok
            if n == 1:
            # Euler 2d order scheme for the first time step
                C[1:-1,n] = C[1:-1,n-1] - CFL/2*(C[2:,n-1] - C[:-2,n-1])   # rev ok

            elif n > 1:
                C[1:-1, n] = C[1:-1, n-2] - CFL*(C[2:, n-1] - C[:-2, n-1])   # rev ok
            
            # Condição de fronteira
            # ---------------------
            if cond_front == 'fixa': # ok
                C[[0,-1], n] = [0,0]
            
            elif cond_front == 'periódica': # ok
                C[-1, n] = C[ 99, n-1] 
                C[ 0, n] = C[-1, n-1] 
      
            elif cond_front == 'radiacional': # ok
                C[-1, n] = C[-1, n-1] - CFL*(C[-1, n-1]- C[-2, n-1])
                 
       # -------------------------------------------------
       # Fourth order on space leapfrog on time
        elif aprox == "ordem4": 
            if n == 1:
                # Euler 2d order scheme for the first time step
                C[1:-1,n] = C[1:-1,n-1] - CFL/2*(C[2:,n-1] - C[:-2,n-1])
                C[0,n] = C[-1,n-1]
            elif n > 1:
                C[2:-2, n] = C[2:-2,n-2] - CFL/6*(8*(C[3:-1,n-1]-C[1:-3,n-1])-(C[4:,n-1]-C[:-4,n-1]))
                        
            #elif cond_front == 'periódica': # (-1, 0, 1) # second order
                C[-2, n] = C[-2, n-2] - CFL*(C[-1, n-1] - C[-3, n-1])
                C[1,n]   = C[1, n-1] - CFL*(C[1, n-1]- C[0, n-1])
                C[-1, n] = C[99, n-2] - CFL/6*(8*(C[100,n-1]-C[98,n-1])-(C[1,n]-C[97,n-1]))
                C[ 0, n] = C[100, n-2] - CFL*(C[1, n] - C[99, n-1])
 

        # -------------------------------------------------       
        elif aprox == "Matsuno": # somente periódica
            # Start with Euler forward
            u[1:-1,n] = C[1:-1, n-1] - CFL*(C[1:-1, n-1] - C[0:-2, n-1])
            C[1:-1,n] = C[1:-1, n-1] - CFL*(u[1:-1, n-1] - u[0:-2, n-1])

            if cond_front == 'periódica':
                u[-1, n] = C[-1, n-1] - CFL*(C[-1, n-1]- C[-2, n-1])
                C[-1, n] = C[-1, n-1] - CFL*(u[-1, n-1] - u[-2,n-1])
                u[0, n] = u[-1, n-1] #- CFL*(C[-1, n-1]- C[-2, n-1])
                C[0, n] = C[-1, n-1] #- CFL*(u[-1, n-1] - u[-2,n-1])

        elif aprox == "Crank": 
            # spsolve: solve the sparse linear system Ax=B
            # A*C[:,n] = B*C[:, n-1]  
            C[:, n] = spsolve(A, B*C[:, n-1])
         
        # -------------------------------------------------
        elif aprox == "RK4":
            """
            Runge-Kutta scheme
            ------------------
            """
            k1 = deriv(c            , U, dx)
            k2 = deriv(c + dt*0.5*k1, U, dx)
            k3 = deriv(c + dt*0.5*k2, U, dx)
            k4 = deriv(c + dt*k3    , U, dx)

            C[:,n] = c = c + (k1 + 2*k2 + 2*k3 + k4)*dt/6
       
        else:
            print("Error")
         
    return C
    
def plot_sol_num(C, fun, aprox, cond_front, ylabel, dP, Nt, CFL, U, dx, dt):
    
    fig, ax = plt.subplots(1,1, figsize=(8, 5))
    ax.plot(C[:,0], color = 'b', lw=3, label="Condição inicial")

    for n in range(1, Nt):  
        if n % dP == 0:
            ax.plot(C[:,n], color='m', linestyle='dashed', label=f"PT {n}")
    hora = round((n*dt/3600),1)
    ax.plot(C[:,-1], color='g', lw=3, label = f"Final {hora} horas ")
    ax.set_ylabel(ylabel)
    ax.set_title(f"$\Delta$t = {dt} segundos, $\Delta$x = {dx} metros, CFL = {round(CFL,2)}.",loc='left')
    ax.legend(fontsize=8, ncol=2)
    ax.set_xlabel("Pontos da grade")
    ax.text(2,3, r"$\vec U$ = "+ f"{U} m/s.", fontsize=12)
    #ax.text(2,4, f"CFL = {CFL.round(2)}", fontsize=12)
    fig.savefig("fig/" + aprox + "_" + fun[:3] + "_" + cond_front[:3] +"_" + str(hora) +".png", 
                dpi = 300, bbox_inches='tight', facecolor='w')

def fig2gif(Nt, dP, dt, CFL, C_ref, c, ylabel, aprox, nr, cond_front, fun):
    # Fazemos gifs
    filenames = []
    for n in range(0,int(Nt),dP):
        # plot the line chart
        fig, ax = plt.subplots(figsize=[8,6])
        ax.set_title(f"Hora: {round((n*dt/3600),1)}, CFL = {round(CFL,1)}, " + r"$\Delta$t"+ f"= {round(dt,2)}")
        ax.plot(C_ref[:,n], color = "b", label="Sol. Analítica")
        ax.plot(c[:,n], color = 'r', linestyle = 'dashed', label= aprox+f' nr: {nr}')
        ax.set_ylim(-4,15)
        ax.legend()
        ax.set_ylabel(ylabel)
        
        # create file name and append it to a list
        filename = f'fig/gifs/{n}.png'
        filenames.append(filename)
        
        # save frame
        fig.savefig(filename, dpi=300)
        plt.close() # build gif
        
    with imageio.get_writer('gifs/'+fun[:3]+'_' + aprox+'_'+ cond_front[:3] +'.gif', mode='I', duration = 1) as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)
            
    # Remove files
    for filename in set(filenames):
        os.remove(filename)