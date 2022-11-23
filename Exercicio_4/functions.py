""" 
Exercicio 4
-----------
Shallow water 2D
- f = 0 (scenario 1)
- f = cte (scenario 2)
- f = varying Coriolis force (scenario 3)
- Radiational boundary conditions 
"""

import numpy as np
import matplotlib.pyplot as plt

def plot_gridC(Nx, Ny):
    """
    Make Grade-C Arakawa
    Idea original: Daniel Melo Costa. Obrigado pelas dicas.
    """
    x2, dx = np.linspace(-4050*1e3, 4050*1e3, Nx+2, retstep = True)  # Array with x-points
    y2, dy = np.linspace(-4050*1e3, 4050*1e3, Ny+2, retstep = True)
    x = (x2[1:]+x2[0:-1])/2
    y = (y2[1:]+y2[0:-1])/2

    Yu, Xu =np.meshgrid(y ,x2)
    Yv, Xv =np.meshgrid(y2, x)
    yv = Yv[0,:]
    Yh, Xh =np.meshgrid(y, x)
    
    """Plot"""
    fig, (ax, ax2) = plt.subplots(2,1, figsize=(8,16), gridspec_kw={'hspace':0.05})
    ax.scatter(Xu/1000, Yu/1000, color='blue', marker='.', s=4,label='u')
    ax.scatter(Xv/1000, Yv/1000, color='red' , marker='.', s=4,label='v')
    ax.scatter(Xh/1000, Yh/1000, color='green', marker='.', s=4,label='h')
    ax.set_ylabel('Y [km]',fontsize=10)
    ax.set_xlim(-4600, 4600)
    ax.set_ylim(-4600, 4600)
    ax.legend(ncol=3, bbox_to_anchor=(0.6, .94))
    ax.set_title("Grade C (Arakawa), $\Delta x = \Delta y$ = {:.0f} km, nº pontos: u ".format(dx/1000)+ f"{Xu.shape}, v {Xv.shape}, h {Xh.shape}", 
                 fontsize= 10, loc='left')
    
    ax2.scatter(Xu/1000, Yu/1000, color='blue', marker='.', s=40,label='u')
    ax2.scatter(Xv/1000, Yv/1000, color='red' , marker='.', s=40,label='v')
    ax2.scatter(Xh/1000, Yh/1000, color='green', marker='.', s=40,label='h')
    ax2.set_xlabel('X [km]',fontsize=10)
    ax2.set_ylabel('Y [km]',fontsize=10)
    ax2.set_xlim(-490, 490)
    ax2.set_ylim(-490, 490)
    fig.savefig("fig/GridC.png", dpi=300, bbox_inches='tight', facecolor='w')
        
    return y, yv, Xu, Yu, Xv, Yv, Xh, Yh, dx, dy
    
def conditions(scen, Xu, Yu, Xh, Yh, yv, lat, Nx, Ny, dx, dy, Nt, Nrx, Nry, test=False):
    """Coriolis parameter
    Args:
        scen (character): scenario name (scen1 = 0, scen2 = f cte, scen3 = beta variation)
        y (array): grid in km
        Nrx = 10  Length wave in x axis
        Nry = 4   Length wave in y axis
    Returns:
        Array: Coriolis value in [1/s].
    """
    beta = np.zeros((Ny+2))
    div  = np.zeros((Nt+1,Nx+1,Ny+1))
    vor  = np.zeros((Nt+1,Nx+1,Ny+1))
     
    # Allocating arrays and initial conditions as Arakawa C-grid
    # ----------------------------------------------------------
    u  = np.zeros((Nt+1, Nx+2, Ny+1))  # zonal wind
    v  = np.zeros((Nt+1, Nx+1, Ny+2))  # meridional wind
    h  = np.zeros((Nt+1, Nx+1, Ny+1))  # height variation by time
    Fu = np.zeros((Nt+1, Nx+2, Ny+1))  # Zonal force like ENSO
    
    if test == False:
        Fu[:,:,:] = -np.exp(-Xu**2/(Nrx*dx)**2 - Yu**2/(Nry*dy)**2)/(24*3600)
    
        # Initial conditions for u, v, & h
        # --------------------------------
        u[0, :, :] = Fu[0,:,:]
        v[0, :, :] = 0
        h[0, :, :] = 0 #40*np.exp(-Xh**2/(Nrx*dx)**2 - Yh**2/(Nry*dy)**2)
        
    else:
        # Initial conditions for u, v, & h
        # --------------------------------
        u[0, :, :] = 0
        v[0, :, :] = 0
        h[0, :, :] = 10*np.exp(-Xh**2/(3*dx)**2 - Yh**2/(3*dy)**2)

    q = 2*np.pi/86400      # angular velocity (2pi/(24*3600 s)) [1/s]
    f = np.zeros((Ny+2))
    a = 6371000            # Earth radius [m]
        
    if scen == 'scen1':
        f = 0
        
    elif scen == 'scen2':
        f = 2*q*np.sin(np.deg2rad(lat))
        
    elif scen == 'scen3':
        lat_deg = 0 #yv/111139  # meters to degree
        beta = 2*q*np.cos(np.deg2rad(lat_deg))/a # plano beta equatorial
            
    return u, v, h, Fu, f, beta, div, vor

def leapfrog(scen, Nt, u, v, h, g, dx, dy, dt, H, Fu, f, beta, y, yv, c, div, vor, rad):
    """Leapfrog scheme to resolve shallow water equation 2D.
    Args:
        Nx (int): n° points in x axes
        Ny (int): n° points in y axes
        Nt (int): n° time steps
        u (array): zonal wind vectors
        v (array): meridional wind vectors
        y (array): y-axis dimension for u and h
        yv (array): y-axis dimension for v
        h (array): height 
        g (float): gravitational acceleration
        dx (int): spatial horizontal resolution at x
        dy (int): spatial horizontal resolution at y 
        dt (int): temporal resolution in seconds
        H (int): Mean height
        Fu (array): Constant source at zonal wind.
        f (array): coriolis
        c (float): sqrt(gH)
        div (array): Divergence
        vor (array): Vorticity
        rad (boolean): Radiational if true, else boundary fixed = 0
    Returns:
        arrays: u, v, h, div, vor
    """
    for n in range(0, Nt):
        if (scen == 'scen1') | (scen == 'scen2'):
                           
            if n == 1: # Euler-forward scheme
                
                u[1, 1:-1, :] = (u[0, 1:-1, :] +                               # u(i,j)            (80, 81) pontos
                                 dt*(-g*(h[0, 1:, :] - h[0, :-1, :])/dx +      # h(i+1,j) - h(i,j) 
                                     f/4*(v[0, :-1, 1:] +                      # v(i  , j  )           
                                          v[0, 1: , 1:] +                      # v(i+1, j  )          
                                          v[0, 1: ,:-1] +                      # v(i+1, j-1)        
                                          v[0, :-1,:-1]) +                     # v(i  , j-1)        
                                     Fu[0,1:-1,:])) # ok
                
                v[1, :, 1:-1] = (v[0, :, 1:-1] +                               # v(i,j)            (81, 80) pontos
                                 dt*(-g*(h[0, :, 1:] - h[0, :, :-1])/dy -      # h(i,j+1) - h(i,j) 
                                     f/4*(u[0, 1:, :-1] +                      # u(i, j)           
                                          u[0, 1:,  1:] +                      # u(i,j+1)          
                                          u[0,:-1,  1:] +                      # u(i-1,j+1)        
                                          u[0,:-1, :-1]))                      # u(i-1,j)          
                                   ) # ok
                
                h[1, :, :] = (h[0, :, :] -                                     # h(i, j)           (81, 81) pontos
                                dt*H*((u[0, 1:, :] -                           # u(i, j)       
                                       u[0,:-1, :])/dx +                       # u(i-1, j)   
                                      (v[0, :, 1:] -                           # v(i, j)     
                                       v[0, :,:-1])/dy )) #ok                  # v(i, j-1)   
                        
            elif n > 1: # Leap-frog scheme
                u[n+1, 1:-1, :] = (u[n-1, 1:-1, :] +                           # u(i, j)           (80, 81) pontos
                                   2*dt*(-g*(h[n, 1:, :] - h[n, :-1, :])/dx +  # h(i+1,j) - h(i,j) 
                                         f/4*(v[n,:-1, 1:] +                   # v(i  , j  )          
                                              v[n, 1:, 1:] +                   # v(i+1, j  )          
                                              v[n, 1:,:-1] +                   # v(i+1, j-1)        
                                              v[n,:-1,:-1]) +                  # v(i  , j-1)        
                                         Fu[n,1:-1,:])) # ok
                
                v[n+1, :, 1:-1] = (v[n-1, :, 1:-1] +                           # v(i,j)            (81, 80) pontos
                                   2*dt*(-g*(h[n, :, 1:] - h[n, :, :-1])/dy -  # h(i,j+1) - h(i,j) 
                                         f/4*(u[n,  1:, :-1] +                 # u(i, j)           
                                              u[n,  1:,  1:] +                 # u(i,j+1)          
                                              u[n, :-1,  1:] +                 # u(i-1,j+1)        
                                              u[n, :-1, :-1]))                 # u(i-1,j)          
                                   ) # ok
                
                h[n+1, :, :] = (h[n-1, :, :] -                                 # h(i  , j  )       (81, 81) pontos
                                2*dt*H*((u[n, 1:, :] -                         # u(i  , j  )       
                                         u[n,:-1, :])/dx +                     # u(i-1, j  )   
                                        (v[n, :, 1:] -                         # v(i  , j  )     
                                         v[n, :,:-1])/dy ))                    # v(i  , j-1)   
            
            # Radiational as boundary condition (west, east, north, south)
            # -----------------------------------------------------------
            if rad == True:
                # west (left) -> du/dt - fv - c*du/dx = 0
                u[n+1, 0, :] = (u[n, 0, :] + 
                                dt*f*(v[n, 0, 1: ] + 
                                      v[n, 0, :-1])/2 + 
                                c*dt/dx * (u[n, 1, :] - u[n, 0, :]))
                
                # East (right) -> du/dt - fv + c*du/dx = 0
                u[n+1, -1, :]  = 0
                
                # North -> dv/dt + fu + c*dv/dy = 0
                v[n+1, :, -1] = (v[n, :, -1] - 
                                 dt*f*(u[n, 1:, -1] + 
                                       u[n, :-1, -1])/2 - 
                                 c*dt/dy * (v[n, :, -1] - v[n, :, -2]))

                # South -> dv/dt + fu - c*dv/dy = 0
                v[n+1, :, 0] = (v[n, :, 0] - 
                                dt*f*(u[n, 1:, 0] + 
                                      u[n, :-1, 0])/2 + 
                                c*dt/dy * (v[n, :, 1] - v[n, :, 0]))
            
            else:
                u[n+1, 0, :] = 0  # west
                v[n+1, :,-1] = 0  # north
                v[n+1, :, 0] = 0  # south
                u[n+1,-1, :] = 0  # east                    
                                                     
        elif scen == 'scen3': # Beta with two "y" (y and yv for C-grid Arakawa)         
            
            if n == 1: # Euler-forward scheme
                u[1, 1:-1, :] = (u[0, 1:-1, :] +                               # u(i, j)             (80, 81) pontos
                                 dt*(-g*(h[0, 1:, :] - h[0, :-1, :])/dx +      # h(i+1,j) - h(i,j)    
                                     beta/4*(yv[1:]*(v[0,:-1, 1:] +            # yv(j)*[v(i,j) +      
                                                     v[0, 1:, 1:]) +           # v(i+1,j)]      
                                             yv[:-1]*(v[0, 1:,:-1] +           # yv[j-1]*[v[i+1,j-1]+ 
                                                      v[0,:-1,:-1])) +         # v[i,j-1]             
                                     Fu[0,1:-1,:])) # ok
                
                v[1, :, 1:-1] = (v[0, :, 1:-1] +                               # v(i,j)              (81, 80) pontos
                                 dt*(-g*(h[0, :, 1:] - h[0, :, :-1])/dy -      # h(i,j+1) - h(i,j)   
                                     beta/4*(y[:-1]*(u[0, 1:, :-1] +           # y[j]*[u[i,j] +       
                                                     u[0,:-1, :-1]) +          # u(i-1,j)   
                                             y[1: ]*(u[0, 1:,  1:] +           # y[1:]*[u(i,j+1) + 
                                                     u[0, :-1, 1:]))           # u(i-1,j+1)        
                                     ))
                
                h[1, :, :] = (h[0, :, :] -                                     # h(i, j)             (81, 81) pontos
                              dt*H*((u[0, 1:, :] -                             # u(i, j)       
                                     u[0,:-1, :])/dx +                         # u(i-1, j)   
                                    (v[0, :, 1:] -                             # v(i, j)     
                                     v[0, :,:-1])/dy )) #ok                    # v(i, j-1)   
                        
            elif n > 1: # Leap-frog scheme               
                u[n+1, 1:-1, :] = (u[n-1, 1:-1, :] +                           # u(i, j)             (80, 81) pontos
                                   2*dt*(-g*(h[n, 1:, :] - h[n, :-1, :])/dx +  # h(i+1,j) - h(i,j)    
                                         beta/4*(yv[1:]*(v[n,:-1, 1:] +        # yv(j)*[v(i,j) +      
                                                         v[n, 1:, 1:]) +       # v(i+1,j)]    
                                                 yv[:-1]*(v[n, 1:,:-1] +       # yv[j-1]*[v[i+1,j-1]+ 
                                                          v[n,:-1,:-1])) +     # v[i,j-1] 
                                         Fu[n,1:-1,:])) # ok
                
                v[n+1, :, 1:-1] = (v[n-1, :, 1:-1] +                           # v(i,j)              (81, 80) pontos
                                   2*dt*(-g*(h[n, :, 1:] - h[n, :, :-1])/dy -  # h(i,j+1) - h(i,j) 
                                         beta/4*(y[:-1]*(u[n, 1:, :-1] +       # y[j]*[u[i,j] + 
                                                         u[n,:-1, :-1]) +      # u(i-1,j)          
                                                 y[1: ]*(u[n, 1:,  1:] +       # y[1:]*[u(i,j+1) + 
                                                         u[n, :-1, 1:]))       # u(i-1,j+1)       
                                         )) #ok
                
                h[n+1, :, :] = (h[n-1, :, :] -                                 # h(i, j)             (81, 81) pontos
                                2*dt*H*((u[n, 1:, :] -                         # u(i, j)       
                                         u[n,:-1, :])/dx +                     # u(i-1, j)   
                                        (v[n, :, 1:] -                         # v(i, j)     
                                         v[n, :,:-1])/dy ))                    # v(i, j-1)    
            
            # Radiational as boundary condition (west, east, north, south)
            # -----------------------------------------------------------
            if rad == True:
                # west (left) -> du/dt - fv - c*du/dx = 0
                u[n+1, 0, :] = (u[n, 0, :] +
                                dt*beta*(yv[1: ] * v[n, 0, 1: ] + 
                                         yv[:-1] * v[n, 0, :-1])/2 +
                                c*dt/dx * (u[n, 1, :] - u[n, 0, :]))
                
                # East (right) -> du/dt - fv + c*du/dx = 0
                u[n+1, -1, :]  = 0
                
                # North -> dv/dt + fu + c*dv/dy = 0
                v[n+1, :, -1] = (v[n, :, -1] -
                                 dt*beta*y[-1]*(u[n, 1: , -1] + 
                                                u[n, :-1, -1])/2 - 
                                 c*dt/dy * (v[n, :, -1] - v[n, :, -2]))

                # South -> dv/dt + fu - c*dv/dy = 0
                v[n+1, :,  0] = (v[n, :,  0] - 
                                 dt*beta*y[0 ]*(u[n, 1: , 0] + 
                                                u[n, :-1, 0])/2 + 
                                 c*dt/dy * (v[n, :, 1] - v[n, :, 0]))
            
            else:
                u[n+1,  0, :] = 0  # west
                v[n+1,  :,-1] = 0  # north
                v[n+1,  :, 0] = 0  # south
                u[n+1, -1, :] = 0  # east
       
        # Divergence (du/dx + dv/dy)
        # -------------------------
        div[n, :, :] = ((u[n, 1:, :] - u[n, :-1, :])/dx +
                        (v[n, :, 1:] - v[n, :, :-1])/dy)
    
        # Vorticity (dv/dx - du/dy)
        # -------------------------
        uplot = (u[n, 1:, :] + u[n, :-1, :])/2
        vplot = (v[n, :, 1:] + v[n, :, :-1])/2
        
        # vor[n, :1, :1] = ( (vplot[:, 1:] - vplot[:, :-1])/dx -
        #                    (uplot[1:, :] - uplot[:-1, :])/dy )
        
        vor[n, :, :] = np.gradient(vplot, axis=0)/dx - np.gradient(uplot, axis=1)/dy
    
    # Mass, energy and enstrophy conservation
    # ---------------------------------------
    mass   = np.zeros((Nt+1))
    Ep     = np.zeros((Nt+1))
    Ek     = np.zeros((Nt+1))
    #enstro = np.zeros((Nt+1))
    
    for n in range(0, Nt):
        mass[n]   = np.nansum(h[n,:,:] * dx*dy)
        
        Ep[n]     = g/2* np.nansum( (h[n,:,:]**2) * dx*dy)
        Ek[n]     = H/2*(np.nansum( (u[n,:,:]**2) * dx*dy) + 
                         np.nansum( (v[n,:,:]**2) * dx*dy))
        
        #enstro[n] = np.nansum(vor[n]**2)
        
    return u, v, h, div, vor, mass, Ep, Ek