""" 
Exercício 4
-----------
Shallow water 2D
- f = 0 (scenario 1)
- f = cte (scenario 2)
- f = varying Coriolis force (scenario 3)
- Radiational boundary conditions 
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import imageio.v2 as imageio
import os

def plot_gridC():
    """
    Make Grid-C Arakawa
    Idea original: Daniel Melo Costa. Obrigado pelas dicas.
    """
    dx = 100*1e3
    dy = 100*1e3

    x2 = np.arange(-4000*1e3, 4000*1e3+dx, dx)
    y2 = np.arange(-4000*1e3, 4000*1e3+dy, dy)

    x = (x2[1:] + x2[ :-1])/2
    y = (y2[1:] + y2[ :-1])/2

    Xu, Yu =np.meshgrid(x2, y )
    Xv, Yv =np.meshgrid( x, y2)
    Xh, Yh =np.meshgrid( x, y )
    
    """Plot"""
    fig, ax = plt.subplots(1, figsize=(5.5,5), gridspec_kw={'hspace':0.05})
    ax.set_title("Grade C (Arakawa), $\Delta x = \Delta y$ = {:.0f} km, \nnº pontos: u ".format(dx/1000)+ f"{Xu.shape}, v {Xv.shape}, h {Xh.shape}\n\n", 
                 fontsize= 10, loc='left')
    
    ax.scatter(Xu/1000, Yu/1000, color='blue', marker='.', s=200,label='u')
    ax.scatter(Xv/1000, Yv/1000, color='red' , marker='.', s=200,label='v')
    ax.scatter(Xh/1000, Yh/1000, color='green', marker='.', s=200,label='h')
    ax.set_xlabel('x [km]',fontsize=10)
    ax.set_ylabel('y [km]',fontsize=10)
    ax.set_xlim(-120, 120)
    ax.set_ylim(-120, 120)
    ax.legend(ncol=3, bbox_to_anchor=(0.51, 1.1), loc='upper right')
    fig.savefig("fig/GridC_ex4.png", dpi=300, bbox_inches='tight', facecolor='w')
        
    return Xu, Yu, Xv, Yv, Xh, Yh, dx, dy

def animation(scen, Xh, Yh, out, H, test, f_space, f_time, pts= 6,
              lv1=1, lv2=1, lv3=1, arr1=0.001, arr2=0.001, arr3=0.001):
            
            lvs1 = list(np.linspace(-lv1, lv1, 16*2+1))
            lvs2 = list(np.linspace(-lv2, lv2, 16*2+1))
            lvs3 = list(np.linspace(-lv3, lv3, 16*2+1))
            lvs = {'scen1': lvs1, 'scen2':lvs2, 'scen3':lvs3}
            arrows = {'scen1': arr1, 'scen2': arr2, 'scen3':arr3}
            staged_u = np.zeros(Xh.shape)
            staged_v = np.zeros(Xh.shape)

            filenames = []

            for n in range(0, len(out['h'])):
                fig, ax = plt.subplots()
                
                if test == True:
                    arrows = {'scen1': 0.001, 'scen2':0.001, 'scen3':0.001}
                    im = ax.contourf(Xh/1000, Yh/1000, out['h'][int(n)], 
                                    #cmap = plt.cm.RdBu_r, 
                                    levels = lvs[scen], extend='both')
                else:
                    arrows = {'scen1': 0.001, 'scen2':0.001, 'scen3':0.001}
                    im = ax.contourf(Xh/1000, Yh/1000, out['h'][int(n)], 
                                    cmap = plt.cm.RdBu_r, 
                                    levels = lvs[scen], extend='both')
                    
                fig.colorbar(im, ax=ax, orientation="vertical")
                ax.set_title("t={:.0f} dias, passo de tempo = {:.0f}".format(out['t'][int(n)]/(24*3600), n), 
                            fontsize=10, loc='left')
                ax.set_xlim(-4000, 4000)
                ax.set_ylim(-4000, 4000)
                
                u, v = out['u'][n], out['v'][n]
                staged_u[ :, :] = (u[ :, 1:] + u[ :, :-1]) * 0.5
                staged_v[ :, :] = (v[ 1:, :] + v[:-1, : ]) * 0.5
                
                Q = ax.quiver(Xh[::pts, ::pts]/1000, Yh[::pts, ::pts]/1000,
                              staged_u[::pts, ::pts], staged_v[::pts, ::pts],
                              units = "xy", scale = arrows[scen], width = 30)

                # create file name and append it to a list
                filename = f'fig/gifs/{n}.png'
                filenames.append(filename)
                plt.close() # build gif

                # save frame
                fig.savefig(filename, dpi=200)
                

            if test == True:
                if (f_space == True) & (f_time == False):
                    ngif = 'gifs/lf_t_'+str(H)+'m_'+scen+'_Fh_fs.gif'

                elif (f_space == False) & (f_time == True):
                    ngif = 'gifs/lf_t_'+str(H)+'m_'+scen+'_Fh_ft.gif'
                    
                elif (f_space == True) & (f_time == True):
                     ngif = 'gifs/lf_t'+str(H)+'m_'+scen+'_Fh_fst.gif'
                    
                else:
                    ngif = 'gifs/lf_t_'+str(H)+'m_'+scen+'_Fh.gif'

                with imageio.get_writer(ngif, mode='I', duration = .10) as writer:
                    for filename in filenames:
                        image = imageio.imread(filename)
                        writer.append_data(image)

            else:
                if (f_space == True) & (f_time == False):
                    ngif = 'gifs/lf_'+str(H)+'m_'+scen+'_Fu_fs.gif'

                elif (f_space == False) & (f_time == True):
                    ngif = 'gifs/lf_'+str(H)+'m_'+scen+'_Fu_ft.gif'
                    
                elif (f_space == True) & (f_time == True):
                     ngif = 'gifs/lf_'+str(H)+'m_'+scen+'_Fu_fst.gif'
                    
                else:
                    ngif = 'gifs/lf_'+str(H)+'m_'+scen+'_Fu.gif'
                    
                with imageio.get_writer(ngif, mode='I', duration = .10) as writer:
                    for filename in filenames:
                        image = imageio.imread(filename)
                        writer.append_data(image)
                        
            # Remove files
            for filename in set(filenames):
                os.remove(filename) 
            print(f"Execution of animation done for {scen}.")
            

def swe_lf_lw(scen, Xu, Yu, Xv, Yv, Xh, Yh, dx, dy, LF=True, LW=False, Nrx=10, Nry=4, H=250,
                 dt = 120, Nt=2160, lat=-20, test=True, rad=True, amp=1, pts=6, splot=True,
                 gamma=0.01, alfa=0.53, f_time=True, f_space=True, sample_interval=60):
    """Leap-frog scheme applied to shallow-water equations.

    Args:
        scen (string): Coriolis scenario
        Xu (array): u arrange in x axis
        Yu (array): u arrange in y axis
        Xv (array): v arrange in x axis
        Yv (array): v arrange in y axis
        Xh (array): h arrange in x axis
        Yh (array): h arrange in y axis
        dx (number): x resolution in meters
        dy (number): y resolution in meters
        H (int, optional): Height. Defaults to 250.
        dt (int): Time step [s].
        Nt (number, optional): _description_. Defaults to 2160 = 3*24*30.
        lat (int, optional): _description_. Defaults to -20.
        test (bool, optional): H in time=0 to test radiational conditions. Defaults to True.
        rad (bool, optional): Radiational. Defaults to True.
        gamma (float, optional): Robert-Asselin factor. Defaults to 0.009.
        f_time (bool, optional): Filter in time. Defaults to True.
        f_space (bool, optional): Filter in space. Defaults to True.
        sample_interval (int): Time to save information. Defaults = 20 s (each 20 seconds)
    """
    Nx, Ny = Xh.shape
    q = 2*np.pi/86400      # angular velocity (2pi/(24*3600 s)) [1/s]
    a = 6371000            # Earth radius [m]
    g = 9.8                # m/s^2 (gravity acceleration)
    c = np.sqrt(g*H)       # Phase velocity [m/s]
    
    if scen == 'scen1':
        f = 0
        fu = f*Yu
        fv = f*Yv
        
    elif scen == 'scen2':
        f = 2*q*np.sin(np.deg2rad(lat))
        fu = np.ones(Yu.shape)*f
        fv = np.ones(Yv.shape)*f
        
    elif scen == 'scen3':
        beta = 2*q*np.cos(np.deg2rad(0))/a
        fu = beta*Yu
        fv = beta*Yv
        
    # We define initial conditions and variables
    """ time definitions
    <var>nm1  (n-1)
    <var>n    ( n )
    <var>nm   ( )
    <var>np1  (n+1) 
    """
    unm, vnm, hnm = np.zeros(Xu.shape), np.zeros(Xv.shape), np.zeros(Xh.shape)
    un, vn, hn = np.zeros(Xu.shape), np.zeros(Xv.shape), np.zeros(Xh.shape)
    unp1, vnp1, hnp1 = un.copy(), vn.copy(), hn.copy()
    
    # Plot source -----------------------------------------------------------#
    staged_u = np.zeros(Xh.shape)
    staged_v = np.zeros(Xh.shape)
    
    if test == True:
        Fu = 0*np.exp(-Xu**2/(Nrx*dx)**2 - Yu**2/(Nry*dy)**2)/(24*3600)
        Fh = amp*np.exp(-(Xh - 15*dx)**2/(4*dx)**2 - Yh**2/(10*dy)**2)
        hn += Fh
        
        if splot == True:
            fig, ax = plt.subplots(figsize=(6.5,5))
            im = ax.contourf(Xh/1000, Yh/1000, hn, cmap='bwr',
                            levels=list(np.linspace(.001,amp,50))) #cmap='bwr'
            staged_u[ :, :] = (un[ :, 1:] + un[ :, :-1]) * 0.5
            staged_v[ :, :] = (vn[ 1:, :] + vn[:-1, : ]) * 0.5        
            ax.quiver(Xh[::pts, ::pts]/1000, Yh[::pts, ::pts]/1000, 
                    staged_u[::pts, ::pts], staged_v[::pts, ::pts],
                    units = "xy", scale = 0.0001, width = 50)
            ax.set_title('(a) Fonte de massa para o teste (Ex. 4)')
            fig.colorbar(im, ax=ax, orientation='vertical') # pad, shrink
            fig.savefig("fig/Fh_test.png", dpi = 250, bbox_inches='tight', facecolor='w')
        else:
            pass

    else:
        Fu = amp*np.exp(-Xu**2/(Nrx*dx)**2 - Yu**2/(Nry*dy)**2)/(24*3600)
        Fh = 0*np.exp(-Xh**2/(3*dx)**2 - Yh**2/(3*dy)**2)
        un += Fu
        
        if splot == True:
            fig, ax = plt.subplots(figsize=(6.5,5))
            im = ax.contourf(Xu/1000, Yu/1000, un*86400, levels=list(np.linspace(.001,amp,50))) #cmap='bwr'
            staged_u[ :, :] = (un[ :, 1:] + un[ :, :-1]) * 0.5
            staged_v[ :, :] = (vn[ 1:, :] + vn[:-1, : ]) * 0.5        
            ax.quiver(Xh[::pts, ::pts]/1000, Yh[::pts, ::pts]/1000, 
                    staged_u[::pts, ::pts], staged_v[::pts, ::pts],
                    units = "xy", scale = 0.001, width = 50)
            ax.set_title('(b) Fonte zonal constante (Ex. 4)')
            fig.colorbar(im, ax=ax, orientation='vertical') # pad, shrink
            fig.savefig("fig/Fu.png", dpi = 250, bbox_inches='tight', facecolor='w')
        
        else:
            pass
        
    # Loop simulation Euler FTBS at n = 1 and thereafter Leapfrog ---------------------------------------#
    # ---------------------------------------------------------------------------------------------------#
    time_step = 1
    t_0 = time.perf_counter() # timing computation loop
    
    # Sampling variables each sample_interval
    h_lst = list(); u_lst = list(); v_lst = list()      # h, u, v sampling
    div_lst = list() ; vor_lst = list()                 # divergence and vorticity
    t_lst  = [0]                                        # sampling time
    mass   = np.array([])                               # mass conservation
    en_k   = np.array([])                               # kinetic energy
    en_p   = np.array([])                               # potential energy
    enst   = np.array([])                               # enstrophy energy

    h_lst.append(hn)
    u_lst.append(un)
    v_lst.append(vn)
        
    # Calculation:
    """ time definitions for <var> as u, v, h
    <var>nm1  (n-1)
    <var>n    ( n )
    <var>np1  (n+1)
    
    Python array structure
        (i, j) -> i as rows    (variation in Y)
                -> j as columns (variation in X)
    """
    # Iterations ----------------------------------------------------------------------
    if c*dt/dx <= 0.35:    
        # 0. Euler-forward scheme at n = 1
        # --------------------------------
        if scen == 'scen1':
                    
            unp1[:, 1:-1] = un[:, 1:-1] + dt*(-g *(hn[: , 1: ] - hn[:  , :-1])/dx +
                                              Fu[:, 1:-1])
            
            vnp1[1:-1, :] = vn[1:-1, :] + dt*(-g *(hn[1:, :  ] - hn[:-1, :  ])/dy)
            
        elif scen == 'scen2':
        
            unp1[:, 1:-1] = un[:, 1:-1] + dt*(-g *(hn[:   , 1: ] - hn[: , :-1])/dx +
                                              f/4*(vn[1:  , :-1] + vn[1: ,  1:] + 
                                                   vn[:-1 , 1: ] + vn[:-1, :-1]) + 
                                              Fu[:, 1:-1])
            vnp1[1:-1, :] = vn[1:-1, :] + dt*(-g *(hn[1:, :] - hn[:-1, :])/dy -
                                              f/4*(un[:-1, 1:] + un[1:, 1:] +
                                                   un[1:, :-1] + un[:-1, :-1]))
            
        elif scen == 'scen3':
            
            unp1[:, 1:-1] = un[:, 1:-1] + dt*(-g *(hn[:  , 1:] - hn[: , :-1])/dx +
                                              1/4*((vn*fv)[1: , :-1] + (vn*fv)[1:  , 1: ] + 
                                                   (vn*fv)[:-1,  1:] + (vn*fv)[:-1 , :-1]) +
                                              Fu[:, 1:-1])
            vnp1[1:-1, :] = vn[1:-1, :] + dt*(-g *(hn[1:, :] - hn[:-1, :])/dy -
                                              1/4*((un*fu)[:-1, 1:] + (un*fu)[1:, 1:] +
                                                   (un*fu)[1:, :-1] + (un*fu)[:-1, :-1]))
            
        hnp1[:, :] = hn[:, :] - dt*H*((un[:, 1:] - un[:  , :-1])/dx +
                                      (vn[1:, :] - vn[:-1, :  ])/dy)
                
        unm1, vnm1, hnm1 = un.copy()  , vn.copy()  , hn.copy()
        un  , vn  , hn   = unp1.copy(), vnp1.copy(), hnp1.copy()

        while (time_step < Nt):
            # 1. Leap-frog scheme at n > 1
            # ----------------------------
            # vector
            if scen == 'scen1':
                unp1[:, 1:-1] = (unm1[:, 1:-1] + 2*dt*(-g * (hn[:, 1:] - hn[:, :-1])/dx +
                                                       Fu[:, 1:-1]))
                
                vnp1[1:-1, :] = (vnm1[1:-1, :] + 2*dt*(-g * (hn[1:, :] - hn[:-1, :])/dy))
                
            elif scen == 'scen2':
                unp1[:, 1:-1] = (unm1[:, 1:-1] + 2*dt*(-g * (hn[:, 1:] - hn[:, :-1])/dx +
                                                       f/4 * (vn[1:, :-1] + vn[1:, 1:] +
                                                              vn[:-1, 1:] + vn[:-1, :-1]) +
                                                        Fu[:, 1:-1]))
                
                vnp1[1:-1, :] = (vnm1[1:-1, :] + 2*dt*(-g*(hn[1:, :] - hn[:-1, :])/dy -
                                                       f/4*(un[:-1, 1:] + un[1: , 1: ] +
                                                            un[1:, :-1] + un[:-1, :-1])))
            
            elif scen == 'scen3':
                unp1[:, 1:-1] = (unm1[:, 1:-1] + 2*dt*(-g*(hn[:, 1:] - hn[:, :-1])/dx +
                                                        1/4*((vn*fv)[1:, :-1] + (vn*fv)[1:, 1:] +
                                                             (vn*fv)[:-1, 1:] + (vn*fv)[:-1, :-1]) +
                                                        Fu[:, 1:-1]))
                
                vnp1[1:-1, :] = (vnm1[1:-1, :] + 2*dt*(-g*(hn[1:, :] - hn[:-1, :])/dy -
                                                        1/4*((un*fu)[:-1, 1:] + (un*fu)[1:, 1:] +
                                                            (un*fu)[1:, :-1] + (un*fu)[:-1, :-1])))
                
            hnp1[:, :] = (hnm1[:, :] - 2*dt*H*((un[:, 1:] - un[:, :-1])/dx +
                                                (vn[1:, :] - vn[:-1, :])/dy ))           
                       
            # 2. Boundary condition
            # ---------------------
            if rad == False:
                unp1[: ,  0] = 0  # west
                unp1[: , -1] = 0  # east
                vnp1[-1, : ] = 0  # north
                vnp1[ 0, : ] = 0  # south
            
            else:
                """ Radiational conditions & rigid at east"""
                
                # East (right) -> du/dt - fv + c*du/dx = 0
                unp1[: , -1] = 0   # Rigid conditions

                if scen == 'scen1':
                    # radiational conditions
                    # ----------------------
                    # west (left) -> du/dt - c*du/dx = 0
                    unp1[: ,  0] = (un[:, 0 ] + c*dt/dx * (un[:, 1] - un[:, 0]))

                    # North -> dv/dt + c*dv/dy = 0
                    vnp1[-1, : ] = (vn[-1, :] - c*dt/dy * (vn[-1, :] - vn[-2, :]))
                                
                    # South -> dv/dt - c*dv/dy = 0
                    vnp1[ 0, : ] = (vn[0, : ] + c*dt/dy * (vn[1, :] - vn[0, :]))
                    
                elif scen == 'scen2':
                    # radiational conditions
                    # ----------------------
                    # west (left) -> du/dt - fv - c*du/dx = 0
                    unp1[: ,  0] = (un[:, 0] + dt*f*(vn[1: , 0] + vn[:-1, 0])/2 +
                                    c*dt/dx * (un[:, 1] - un[:, 0]))

                    # North -> dv/dt + fu + c*dv/dy = 0
                    vnp1[-1, : ] = (vn[-1, :] - dt*f*(un[-1, 1: ] + un[-1, :-1])/2 - 
                                    c*dt/dy * (vn[-1, :] - vn[-2, :]))
                                
                    # South -> dv/dt + fu - c*dv/dy = 0
                    vnp1[ 0, : ] = (vn[0, :] - dt*f*(un[0, 1: ] + un[0, :-1])/2 +
                                    c*dt/dy * (vn[1, :] - vn[0, :]))
            
                elif scen == 'scen3':
                    # radiational conditions
                    # ----------------------
                    # west (left) -> du/dt - fv - c*du/dx = 0
                    unp1[: ,  0] = (un[:, 0] + dt*((vn*fv)[1: , 0] + (vn*fv)[:-1, 0])/2 +
                                    c*dt/dx * (un[:, 1] - un[:, 0]))

                    # North -> dv/dt + fu + c*dv/dy = 0
                    vnp1[-1, : ] = (vn[-1, :] - dt*((un*fu)[-1, 1: ] + (un*fu)[-1, :-1])/2 - 
                                    c*dt/dy * (vn[-1, :] - vn[-2, :]))
                                
                    # South -> dv/dt + fu - c*dv/dy = 0
                    vnp1[ 0, : ] = (vn[0, :] - dt*((un*fu)[0, 1: ] + (un*fu)[0, :-1])/2 +
                                    c*dt/dy * (vn[1, :] - vn[0, :]))

            # 3. Apply a Robert-Asselin filter to suppress the computational mode
            # -------------------------------------------------------------------
            """Filtro Robert-Asselin-Williams para remover o modo computacional
            gamma = 0.1 ou 0.01
            alfa = 0.53
            """
            if f_time == True:
                un = un + gamma*alfa/2*(unm1 - 2*un + unp1)
                unp1 = unp1 - gamma*(1-alfa)/2*(unm1 - 2*un + unp1)

                vn = vn + gamma*alfa/2*gamma*(vnm1 - 2*vn + vnp1)
                vnp1 = vnp1 - gamma*(1-alfa)/2*(vnm1 - 2*vn + vnp1)
            
            elif f_space == True:  # only index not vector
                for i in range(1, Nx):
                    for j in range(1, Ny-1):    
                        un[j, i] = un[j, i] + gamma*alfa/2 * (un[j+1, i+1] - 2*un[j, i] + un[j-1, i-1])
                        un[j+1, i+1] = un[j+1, i+1] - gamma*(1-alfa)/2*(un[j+1, i+1] - 2*un[j, i] + un[j-1, i-1])
                        
                for i in range(1, Nx-1):
                    for j in range(1, Ny):                            
                        vn[j, i] = vn[j, i] + gamma*alfa/2*(vn[j+1, i+1] - 2*vn[j, i] + hn[j-1, i-1])
                        vn[j+1, i+1] = vn[j+1, i+1] - gamma*(1-alfa)/2*(vn[j+1, i+1] - 2*vn[j, i] + vn[j-1, i-1])
                        
                for i in range(1, Nx-1):
                    for j in range(1, Ny-1):                            
                        hn[j, i] = hn[j, i] + gamma*alfa/2*(hn[j+1, i+1] - 2*hn[j, i] + hn[j-1, i-1])
                        hn[j+1, i+1] = hn[j+1, i+1] - gamma*(1-alfa)/2*(hn[j+1, i+1] - 2*hn[j, i] + hn[j-1, i-1]) 
            else:
                pass
            
            # 4. Store the resulting fields at regular time
            # -----------------------------------------------------------------------------                   
            # 4.1 Divergence (du/dx + dv/dy)
            # ------------------------------
            div = ((un[:, 1:] - un[:, :-1])/dx +
                   (vn[1:, :] - vn[:-1, :])/dy)
            
            # 4.2 Vorticity (dv/dx - du/dy)
            # -----------------------------
            uplot = (un[:, 1:] + un[:, :-1])/2
            vplot = (vn[1:, :] + vn[:-1, :])/2
            vor = np.gradient(vplot, axis=1)/dx - np.gradient(uplot, axis=0)/dy
            # vor = ((vn[1:-1, 1:] - vn[1:-1, :-1])/dx +   # 81, 79
            #        (un[1:, 1:-1] - un[:-1, 1:-1])/dy)    # 79, 81
            
            # 4.3 Conservation mass and energy
            # --------------------------------
            mass = np.concatenate((mass, [np.nansum(hn * dx*dy)]))
            en_p = np.concatenate((en_p, [g/2 * np.nansum(hn**2)* dx*dy]))
            en_k = np.concatenate( (en_k, [H/2 * (np.nansum((un**2) * dx*dy ) + 
                                                  np.nansum((vn**2) * dx*dy ))]))
                                        
            # 4.4 Samples each sample interval
            # ----------------------------
            if (time_step % sample_interval == 0):
                print("Time: \t{:.2f} hours".format(time_step*dt/3600))
                print("Step: \t{} / {}".format(time_step, Nt))
                print("Mass: \t{:.3f}\n".format(np.sum(hn)))
                h_lst.append(hn)
                u_lst.append(un)
                v_lst.append(vn)
                div_lst.append(div)
                vor_lst.append(vor)
                t_lst.append(time_step*dt)        # seconds
                
            # 5. Switch the time-step results
            # ---------------------------------  
            unm1, vnm1, hnm1 = un.copy(), vn.copy(), hn.copy()       # n   -> n-1
            un  , vn  , hn   = unp1.copy(), vnp1.copy(), hnp1.copy() # n+1 -> n       
                
            time_step += 1 
        
        # End the time loop of the model and the entire model code --------------------------------------                
        out = {'u':u_lst, 'v':v_lst, 'h':h_lst, 't':t_lst, 
               'div':div_lst, 'vor':vor_lst, 
               'm': mass, 'ek':en_k, 'ep':en_p}
  
        print(f"Main computation loop done for {scen}!."+"\nExecution time: {:.2f} s".format(time.perf_counter() - t_0)) 
        print("Courant number {:.2f}".format(c*dt/dx))
                
        return out
    
    else:
        return print("Ops CFL > 0.35, please reduce dt")
    
    
def plot_swe(cenarios, data, Xh, Yh, Xu, Yu, H, nidx, test, 
             plot_time = True, plot_div_vor = True,
             width=30, pts=3, lv1=.5, lv2=.5, lv3=.5, lv_div = 5*1e-7, lv_vor = 1e-5,
             arr1=0.001, arr2=0.001, arr3=0.001, arr_div_amp = 1, arr_vor_amp = 1):
    
    if plot_time == True:
        # Arrows
        staged_u = np.zeros(Xh.shape)
        staged_v = np.zeros(Xh.shape)

        lvs1 = list(np.linspace(-lv1, lv1, 16*2+1))
        lvs2 = list(np.linspace(-lv2, lv2, 16*2+1))
        lvs3 = list(np.linspace(-lv3, lv3, 16*2+1))
        lvs = {'scen1': lvs1, 'scen2':lvs2, 'scen3':lvs3}

        arrows = {'scen1': arr1, 'scen2':arr2, 'scen3':arr3}
        titles = ['f = 0', 'f em -20°S', 'Beta equatorial']
        nidxs = [int(nidx/2), int(nidx/2), int(nidx/2), nidx, nidx, nidx ]
            
        fig, axes = plt.subplots(3,2, figsize=(5.5,12), sharex=True, sharey=True, gridspec_kw={'wspace':0.15, 'hspace':0.15} )
        for ax, scen, tit, nidx in zip(axes.flatten(), cenarios*2, titles*2, nidxs):
            
            im = ax.contourf(Xh/1000, Yh/1000, data[scen]['h'][nidx], levels=lvs[scen], extend='both') #cmap=plt.cm.RdBu_r
            cbar = fig.colorbar(im, ax=ax, orientation="vertical")
            cbar.ax.set_ylabel('h [m]', fontsize=8)
            u, v = data[scen]['u'][nidx], data[scen]['v'][nidx]
            staged_u[ :, :] = (u[ :, 1:] + u[ :, :-1]) * 0.5
            staged_v[ :, :] = (v[ 1:, :] + v[:-1, :]) * 0.5
            Q = ax.quiver(Xh[::pts, ::pts]/1000, Yh[::pts, ::pts]/1000, staged_u[::pts, ::pts], staged_v[::pts, ::pts],
                units = "xy", scale = arrows[scen], width = width)
            #ax.quiver(-3000,3000,1,0,scale=20,width=width)
            ax.set_xlim(-4000, 4000)
            ax.set_ylim(-4000, 4000)
            #qk = ax.quiverkey(Q, 0.9, 0.8, 1000, "0.1 m/s", labelpos = "E", coordinates = "figure")
            ax.set_title(tit+", dia {:.0f}, n = {:.0f}, H = {:.0f}".format(data[scen]['t'][nidx]/(3600*24), nidx, H), fontsize=10, loc='left')
        #cax = fig.add_axes([.05,.07,0.9,.02]) # left, bottom, width, height
        #cbar = fig.colorbar(im, cax=cax, orientation='horizontal') # pad, shrink
        #cbar.ax.tick_params(labelsize=9)
        if test == True:
            fig.savefig("fig/leapfrog_test_"+str(H)+"m_cenarios.png", dpi = 300, bbox_inches='tight', facecolor='w')
            
        else:
            fig.savefig("fig/leapfrog_"+str(H)+"m_cenarios.png", dpi = 300, bbox_inches='tight', facecolor='w')
            
    else:
        pass
            
    if plot_div_vor == True:
        # Divergence & Vorticity
        # ----------------------

        lvs1 = list(np.linspace(-lv_div,lv_div,2*16+1))
        lvs2 = list(np.linspace(-lv_vor, lv_vor,2*16+1))
        titles = ['f = 0', 'f em -20°S', 'Beta equatorial']
            
        fig, axes = plt.subplots(3,2, figsize=(8,12), sharex=True, sharey=True, gridspec_kw={'wspace':0.15, 'hspace':0.15} )
        # Divergencia
        # -----------
        for ax, scen in zip([axes[0,0], axes[1,0], axes[2,0]], cenarios):
            im = ax.contourf(Xh/1000,Yh/1000, data[scen]['div'][nidx], levels=lvs1, cmap=plt.cm.RdBu_r,  extend='both') #,
            #cbar = fig.colorbar(im, ax=ax, orientation="vertical")
            # u, v = data[scen]['u'], data[scen]['v']
            # staged_u[:, :, :] = (u[:, 1:, :] + u[:, :-1, :]) * 0.5
            # staged_v[:, :, :] = (v[:, :, 1:] + v[:, :, :-1]) * 0.5
            Q = ax.quiver(Xu[::pts, ::pts]/1000, Yu[::pts, ::pts]/1000, 
                        data[scen]['u'][nidx][::pts, ::pts], data[scen]['v'][nidx][::pts, ::pts],
                        units = "xy", scale = arrows[scen]*arr_div_amp)
            #qk = ax.quiverkey(Q, 0.9, 0.8, 1000, "0.1 m/s", labelpos = "E", coordinates = "figure")
            ax.set_xlim(-4000, 4000)
            ax.set_ylim(-4000, 4000)
            ax.set_title("Divergencia: dia {:.0f}, H = {:.0f}, ".format(data[scen]['t'][nidx]/(3600*24),  H)+scen, fontsize=10, loc='left')
            cax = fig.add_axes([.125,.07,0.35,.015]) # left, bottom, width, height
            cbar = fig.colorbar(im, cax=cax, orientation='horizontal') # pad, shrink
            cbar.ax.tick_params(labelsize=8)
            cbar.ax.set_xlabel('$\\nabla \\cdot \\vec{V}$'+' [s$^{-1}$]')
            
        # Vorticidade
        # -----------
        for ax, scen in zip([axes[0,1], axes[1,1], axes[2,1]], cenarios):
            im = ax.contourf(Xh/1000,Yh/1000, data[scen]['vor'][nidx], levels=lvs2, cmap=plt.cm.RdBu_r,  extend='both') #levels=lvs,
            #cbar = fig.colorbar(im, ax=ax, orientation="vertical")
            Q = ax.quiver(Xu[::pts, ::pts]/1000, Yu[::pts, ::pts]/1000, 
                        data[scen]['u'][nidx][::pts, ::pts], data[scen]['v'][nidx][::pts, ::pts],
                        units = "xy", scale = arrows[scen]*arr_vor_amp)
            #qk = ax.quiverkey(Q, 0.9, 0.8, 1000, "0.1 m/s", labelpos = "E", coordinates = "figure")
            ax.set_xlim(-4000, 4000)
            ax.set_ylim(-4000, 4000)
            ax.set_title("Vorticidade: dia {:.0f}, H = {:.0f}, ".format(data[scen]['t'][nidx]/(3600*24),  H)+scen, fontsize=10, loc='left')
            cax = fig.add_axes([.55,.07,0.35,.015]) # left, bottom, width, height
            cbar = fig.colorbar(im, cax=cax, orientation='horizontal') # pad, shrink
            cbar.ax.tick_params(labelsize=8)
            cbar.ax.set_xlabel('$\\nabla \\times \\vec{V}$'+' [s$^{-1}$]')

        if test == True:
            fig.savefig("fig/leapfrog_test_div_vor"+str(H)+"m_"+scen+".png", dpi = 300, bbox_inches='tight', facecolor='w')
        else:
            fig.savefig("fig/leapfrog_div_vor"+str(H)+"m_"+scen+".png", dpi = 300, bbox_inches='tight', facecolor='w')
    
    else:
        pass
        
    
def lf_swe(scen, amp=1, H=10, days=100, dt=60, tstep=24):
    """ Shallow water equations not linearized, a simple energy-conserving model
        using the Leapfrog scheme
        
    Args:    
        amp:    1 (Fu activated), 0 (Fu = 0)
        H:      Height for initial conditions
        days:   integration period
        dt:     temporal resolution [s]
        tsetp:  time step [hours]
    """

    dx = 100*1e3   # m
    dy = 100*1e3   # m

    #Grade C
    x2 = np.arange(-4000*1e3, 4000*1e3 + dx, dx)
    y2 = np.arange(-4000*1e3, 4000*1e3 + dy, dy)
    x = (x2[1:] + x2[:-1])/2
    y = (y2[1:] + y2[ :-1])/2

    Xu, Yu = np.meshgrid(x2, y)
    Xv, Yv = np.meshgrid(x, y2)
    Xh, Yh = np.meshgrid(x,  y)
    Xz, Yz = np.meshgrid(x2,y2)
    
    # Define constants
    # ------------------------------------------

    Tf = days*24*3600                         # final time
    t = np.arange(0, Tf, dt)
    q = 2*np.pi/86400                         # angular velocity (2pi/(24*3600 s)) [1/s]
    a = 6371000                               # Earth radius [m]
    g = 9.8                                   # m/s^2 (gravity acceleration)
    c = np.sqrt(g*H)                          # Phase velocity [m/s]
    if scen == 'scen1':
        f = 0
        fu = f*Yu
        fv = f*Yv
    elif scen == 'scen2':
        f = 2*q*np.sin(np.deg2rad(-20))
        fu = np.ones(Yu.shape)*f
        fv = np.ones(Yv.shape)*f
    else:
        beta = 2*q*np.cos(np.deg2rad(0))/a        # Plano Beta equatorial     
        fu = beta*Yu
        fv = beta*Yv

    # We define initial conditions and variables
    """ time definitions
    <var>nm1  (n-1)
    <var>n    ( n )
    <var>np1  (n+1) 
    """

    # time (n-1)
    unm1, vnm1, hnm1  = np.zeros(Xu.shape), np.zeros(Xv.shape), np.zeros(Xh.shape)
    # time ( n )
    un  , vn  , hn    = np.zeros(Xu.shape), np.zeros(Xv.shape), np.zeros(Xh.shape)
    # time (n+1)
    unp1, vnp1, hnp1  = np.zeros(Xu.shape), np.zeros(Xv.shape), np.zeros(Xh.shape)

    # Source conditional for zonal wind and height
    # ------------------------------------------------------------------------------
    alfa = 0.02
    Fu = amp*np.exp(-(Yu-0*dy)**2/(4*dy)**2)*np.exp(-(Xu-15*dx)**2/(10*dx)**2)/(24*3600)
    Fh = 10*np.exp(-Xh**2/(4*dx)**2 - Yh**2/(10*dy)**2)

    dec = (alfa**6) * t[::-1]/(3600*24)**2 * np.exp(-alfa*(t/(3600*24)))
    dec = dec/np.nanmax(dec)
    un   += Fu*dec[0]
    hn   += H + Fh
    hnm1 += H + Fh
    hnp1 += H + Fh

    # Sampling variables each sample_interval
    # ---------------------------------------
    h_lst = list()                                     # h sampling
    u_lst = list(); v_lst = list()                     # u, v sampling
    div_lst = list(); vor_lst = list()                 # divergence and vorticity
    t_lst  = [0]
    sample_interval = 60*tstep                         # Time to save information [min]
    h_lst.append(hn)
    u_lst.append(un)
    v_lst.append(vn)

    U, V = un.copy(), vn.copy()
    mass   = np.array([])                              # mass conservation
    en_k   = np.array([])                              # kinetic energy
    en_p   = np.array([])                              # potential energy
    enst   = np.array([])                              # enstrophy energy

    time_step = 0

    for n in range(1, t.shape[0]):
        """Shallow water equations linearized
        U: Mass flux U
        V: Mass flux V
        B: Bernoulli defined at the locations where h is defined
        Zeta: absolute potential vorticity is redefined to the corners of th C-grid
        """
        # Mass fluxes U and V at inside points related to h
        # -------------------------------------------------
        # Mass flux in U
        U = un.copy()*np.nan
        U[:,1:-1] = un[:,1:-1] * 1/2*(hn[:, :-1] + hn[:, 1:])
        U[:,   0] = un[: ,0] * hn[:,0]
        U[:,  -1] = un[:,-1] * hn[:,-1]

        # Mass flux in V
        V = vn.copy()*np.nan
        V[1:-1, :] = vn[1:-1, :] * 1/2*(hn[:-1, :]+hn[1:, :])
        V[0   , :] = vn[  0, :] * hn[0 , :]
        V[-1  , :] = vn[ -1, :] * hn[-1, :]

        # Bernoulli (ok)
        # ----------------------------------------------------------------------------------------
        B = (g*hn + 1/2 * (1/2*(un[: , 1:]**2 + un[:, :-1]**2) +
                           1/2*(vn[1:,  :]**2 + vn[:-1, :]**2)))

        # ZETA - absolute potential vorticity related with fv
        # ----------------------------------------------------------------------------------------
        ZETA = np.ones((Xz.shape[0], Yz.shape[1]))*np.nan
        
        # Inside points -> E = (f + dv/dx - du/dy)/h
        ZETA[1:-1, 1:-1] = ((fv[1:-1, :-1] + (vn[1:-1, 1:] - vn[1:-1, :-1])/dx - (un[1:, 1:-1] - un[:-1, 1:-1])/dy) / 
                            ((hn[:-1, :-1] + hn[:-1, 1:] + hn[1:, :-1] + hn[1:, 1:])/4))
        
        # South - North [0, -1] in v (Y axis)
        ZETA[ 0, 1:-1] = (fv[ 0, :-1] + (vn[ 0, 1:] - vn[0 , :-1])/dx) / ((hn[ 0, :-1] + hn[ 0, 1:])/2) 
        ZETA[-1, 1:-1] = (fv[-1, :-1] + (vn[-1, 1:] - vn[-1, :-1])/dx) / ((hn[-1, :-1] + hn[-1, 1:])/2)
        
        # East - West [0, -1] in u (X axis)
        ZETA[1:-1, 0] = (fv[1:-1, 0] - (un[1:, 0] - un[:-1 , 0])/dy) / ((hn[:-1, 0] + hn[1:, 0])/2) 
        ZETA[1:-1,-1] = (fv[1:-1,-1] - (un[1:,-1] - un[:-1 ,-1])/dy) / ((hn[:-1,-1] + hn[1:,-1])/2) 

        # Iterations as leap-frog scheme
        # -----------------------------------------------------------------------------------------
        unp1[:   , 1:-1] = (unm1[:,1:-1] + 2*dt*(1/2*(ZETA[1: , 1:-1]* 1/2 *(V[1: , :-1] + V[1: , 1:]) +
                                                      ZETA[:-1, 1:-1]* 1/2 *(V[:-1, :-1] + V[:-1, 1:])) - 
                                                 (B[:, 1:] -B[:, :-1])/dx + Fu[:, 1:-1] * dec[n]))
        
        vnp1[1:-1,:    ] = (vnm1[1:-1,:] - 2*dt*(1/2*(ZETA[1:-1,  1:]* 1/2 *(U[:-1,  1:] + U[1:, 1:]) +
                                                      ZETA[1:-1, :-1]* 1/2 *(U[:-1, :-1] + U[1:,:-1])) + 
                                                 (B[1:, :] - B[:-1, :])/dy))
        
        hnp1[:, :]       = (hnm1[:,:] + 2*dt*(-(U[:, 1:] - U[:, :-1])/dx - (V[1:, :] - V[:-1, :])/dy))
        
        # Boundary conditions (radiational and rigid at east)
        # ------------------------------------------------------------------------------------------
        # East -> du/dt - fv + c*du/dx = 0
        unp1[:,-1] = 0
        
        # West -> du/dt - fv - c*du/dx = 0    (radiational)
        unp1[: , 0] = un[: , 0] + dt*((vn*fv)[1: , 0] + (vn*fv)[:-1, 0])/2 + c*dt/dx * (un[:, 1] - un[:, 0])

        # South -> dv/dt + fu - c*dv/dy = 0   (radiational)
        vnp1[0 , :] = vn[0 , :] - dt*((un*fu)[0,  1:] + (un*fu)[0, :-1])/2 + c*dt/dx * (vn[1, :] - vn[0, :])

        # North -> dv/dt + fu + c*dv/dy = 0   (radiational)
        vnp1[-1, :] = vn[-1, :] - dt*((un*fu)[0,  1:] + (un*fu)[0, :-1])/2 - c*dt/dx * (vn[-1,:] - vn[-2, :])

        # Loop
        # ----------------------------------------------------------------------
        # n   -> n-1
        unm1 = un.copy()
        vnm1 = vn.copy()
        hnm1 = hn.copy()
        
        # n+1 -> n
        un = unp1.copy()
        vn = vnp1.copy()
        hn = hnp1.copy()

        # Divergence and vorticity
        # ----------------------------------------------------------------------    
        uplot = (un[: ,1:] + un[:  , :-1])/2
        vplot = (vn[1:, :] + vn[:-1,   :])/2
        
        div = (un[:,1:] - un[:, :-1])/dx + (vn[1:, :] - vn[:-1, :])/dy
        vor = np.gradient(vplot, axis=1)/dx - np.gradient(uplot, axis=0)/dy

        # Mass & energy conservation
        # ----------------------------------------------------------------------
        mass = np.concatenate((mass, [np.nansum(hn * dx*dy)]))
        en_p = np.concatenate((en_p, [g/2 * np.nansum((hn**2) * dx*dy)]))
        en_k = np.concatenate((en_k, [H/2 *(np.nansum((un**2) * dx*dy) + 
                                            np.nansum((vn**2) * dx*dy))]))
        
        enst = np.concatenate((enst, [np.nansum(vor**2)]))
            
        time_step += 1
        
        if (time_step % sample_interval == 0):
            print("Time: \t{:.2f} hours".format(time_step*dt/3600))
            print("Sample interval (days): \t{} / {:.0f}".format(time_step/sample_interval+1, Tf/(24*60*dt)))
            print("Mass: \t{:.2f}\n".format(np.sum(hn)))
            h_lst.append(hn)
            u_lst.append(un)
            v_lst.append(vn)
            div_lst.append(div)
            vor_lst.append(vor)
            t_lst.append(time_step*dt)        # seconds
    
    param = {'dx':dx,'dy':dy,'Xh':Xh,'Yh':Yh, 'Xu':Xu, 'Yu':Yu, 'Xv':Xv, 'Yv':Yv, 'Xz':Xz, 'Yz':Yz, 'H':H, 'Fh':Fh, 'Fu':Fu}
    
    return {'t':t_lst, 'h': h_lst, 'u':u_lst, 'v': v_lst, 'div':div_lst, 'vor':vor_lst, 'm':mass, 'ep':en_p, 'ek':en_k, 'ens': enst}, param
            
    

        
    
    
    
    