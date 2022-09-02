import numpy as np
import matplotlib.pyplot as plt

def sol_analytical(cond_front, fun, X, U, T, Nx):
    
    if cond_front == "radiacional":
        if fun == "Gaussiana":
            C = {nr:10*np.exp(-(X - U*T - 51*5000)**2/(nr*5000)**2) for nr in [10, 2]}

        if fun == "Retângulo":
            rect = np.linspace(-10, 10, 101)
            C = {np.where(abs(np.linspace(-10, 10, 101)) <= 0.5, 1, 0)*nr for nr in [10, 2]}
       
    elif cond_front == "periódica":
        C = {nr:np.maximum(10*np.exp(-(X - U*T - 51*5000)**2/(nr*5000)**2),
             10*np.exp(-(X - U*T - (51-(Nx-2))*5000)**2/(nr*5000)**2))
             for nr in [10, 2]}

    return C
    
def Conc(i, nr):
    """
    Função de concentração com variação Gaussiana
    ---------------------------------------------
    i:      index
    nr:     número de pontos a partir do centro da gaussiana
    """
    return 10*np.exp(-(i*5000 - 51*5000)**2/(nr*5000)**2)
    
def sol_num(aprox, cond_front, fun, CFL, nr, Nx, Nt, x, dx, U, t):
    """
    aprox       : Tipo de aproximação: "ordem1", "leapfrog", "ordem4", "Matsuno", "implicito", "RK4" (Runge Kutta 4).
    cond_front  : Condição de fronteira, tipo "fixa", "periódica" e "radiacional".
    fun         : Função tipo "Gaussiana" e "Retângulo"
    hora        : Tempo em horas (int)
    """
    C = np.zeros((Nx, Nt))   # Matriz

    # Condições iniciais
    # -------------------
    if fun == "Gaussiana":
        C[:,0] = Conc(x/dx, nr)

    elif fun == "Retângulo":
        C[[int(51 - 2/2), 51, int(51 + 2/2)], 0] = 1/2*nr
            
    for n in range(1, Nt): 
        # Aproximação numérica:
        # ---------------------
        if aprox == "ordem1":
            C[1:Nx, n] = C[1:Nx, n-1] - CFL*(C[1:Nx, n-1]- C[0:Nx-1, n-1])

        elif aprox == "leapfrog":
            # Primeiro passo: Forward Time Central Space
            C[1:Nx-1, 1] = C[1:Nx-1, 0] - CFL/2*(C[2:Nx, 0] - C[0:Nx-2, 0])
            C[0, 1]      = Conc(x[0]/dx - U*t[1], nr)
            Q = (1 - CFL)/(1 + CFL)
            C[1:Nx-1, n] = C[1:Nx-1, n - 2] - CFL*(C[2:Nx, n-1] - C[0:Nx-2, n-1])

            # Boundary condition
            C[0, n]      = Conc((x[0] - U*t[n])/dx, nr)    
            C[Nx-1, n]   = C[Nx - 2, n-1] - Q*C[Nx - 2, n] + Q*C[Nx-1, n-1]

        elif aprox == "Matsuno":
            pass

                    
        # Condição de fronteira
        # ---------------------
        if cond_front == 'fixa':
            C[[0,-1], n-1] = [0,0]

        elif cond_front == 'periódica':
            C[[0, 1], n] = C[[-1,-2],n-1]
                
        elif cond_front == 'radiacional':
            C[ 0, n] = C[ 0, n-1] - CFL*(C[ 1, n-1]- C[ 0, n-1])
            C[-1, n] = C[-1, n-1] - CFL*(C[-1, n-1]- C[-2, n-1])
            
    
    return C
    
def plot_sol_num(fun, aprox, cond_front, ylabel, dP, hora, Nx, Nt, CFL, U, x, dx, t, dt):
    if fun == "Retângulo":
        fig, ax = plt.subplots(1,1, figsize=(8, 5))
        nr = 10
        C = sol_num(aprox, cond_front, fun, CFL, nr, Nx, Nt, x, dx, U, t)

        ax.plot(C[:,0], color = 'b', lw=3, label="Condição inicial")

        for n in range(1, Nt):  
            if n % dP == 0:
                ax.plot(C[:,n], color='m', linestyle='dashed', label=f"PT {n}")

        ax.plot(C[:,-1], color='g', lw=3, label = f"Final {hora} horas ")
        ax.set_ylabel(ylabel)
        ax.set_title("Concentração"+f" para $\Delta$t = {dt.round(0)} segundos, $\Delta$x = {dx} metros.",loc='left')
        ax.legend(fontsize=8, ncol=2)
        ax.set_xlabel("Pontos da grade")
        ax.text(2,3, r"$\vec U$ = "+ f"{U} m/s.", fontsize=12)
        ax.text(2,4, f"CFL = {CFL.round(2)}", fontsize=12)
        fig.savefig("fig/" + aprox + "_" + fun[:3] + "_" + cond_front[:3] +"_" + str(hora) +".png", 
                    dpi = 300, bbox_inches='tight', facecolor='w')

    elif fun == "Gaussiana":
        fig, ax = plt.subplots(2,1, figsize=(10, 8))
        
        for j, nr in enumerate([10, 2]):
            
            C = sol_num(aprox, cond_front, fun, CFL, nr, Nx, Nt, x, dx, U, t)
            ax[j].plot(C[:,0], color = 'b', lw=3, label="Condição inicial")

            for n in range(1, Nt):  
                if n % dP == 0:
                    ax[j].plot(C[:,n], color='m', linestyle='dashed', label=f"PT {n}")

            ax[j].plot(C[:,-1], color='g', lw=3, label = f"Final {hora} horas ")
            ax[j].set_ylabel(ylabel)
            ax[j].set_title("Concentração"+f" para $\Delta$t = {dt.round(0)} segundos, $\Delta$x = {dx} metros e nr = {nr}.",loc='left')
            ax[0].legend(fontsize=8, ncol=2)
            ax[1].set_xlabel("Pontos da grade")
            ax[1].text(2,3, r"$\vec U$ = "+ f"{U} m/s.", fontsize=12)
            ax[1].text(2,4, f"CFL = {CFL.round(2)}", fontsize=12)
            
        fig.savefig("fig/" + aprox + "_" + fun[:3] + "_" + cond_front[:3] +"_" + str(hora) +".png", 
                    dpi = 300, bbox_inches='tight', facecolor='w')