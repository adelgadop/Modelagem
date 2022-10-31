import time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import imageio.v2 as imageio
import os

def conditions(scen, lat, Nx, Ny, dx, dy, Nt, Nrx, Nry):
    """Coriolis parameter
    Args:
        scen (character): scenario name (scen1 = 0, scen2 = f cte, scen3 = beta variation)
        y (array): grid in km
        Nrx = 10  Length wave in x axis
        Nry = 4   Length wave in y axis

    Returns:
        Array: Coriolis value in [1/s].
    """
    x, _ = np.linspace(-(Nx+1)*dx, (Nx+1)*dx, Nx+1, retstep = True)  # Array with x-points
    y, _ = np.linspace(-(Ny+1)*dy, (Ny+1)*dy, Ny+1, retstep = True)  # Array with y-points
    X, Y = np.meshgrid(x, y)                              # Meshgrid for plotting
    X = np.transpose(X)                                   # To get plots right
    Y = np.transpose(Y)                                   # To get plots right

    beta = np.zeros((Ny+1, 1))
    div = np.zeros((Nt+1,Nx+1,Ny+1))
    vor = np.zeros((Nt+1,Nx+1,Ny+1))
     
    # Allocating arrays and initial conditions as Arakawa C-grid
    # ----------------------------------------------------------
    u  = np.zeros((Nt+1, Nx+2, Ny+1))  # zonal wind
    v  = np.zeros((Nt+1, Nx+1, Ny+2))  # meridional wind
    h  = np.zeros((Nt+1, Nx+1, Ny+1))  # height variation by time
    Fu = np.zeros((Nt+1, Nx+1, Ny+1))  # Zonal force like ENSO
    
    Fu[:,:,:] = np.exp(-X**2/(Nrx*4*dx)**2 - Y**2/(Nry*2*dy)**2) 

    q = 2*np.pi/86400      # angular velocity (2pi/24h) [1/s]
    f = np.zeros((Ny+1))
    a = 6371000            # Earth radius [m]
    if scen == 'scen1':
        f[:] = 0
    elif scen == 'scen2':
        f[:] = 2*q*np.sin(np.deg2rad(lat))
    elif scen == 'scen3':
        lat_deg = y/110    # y in km to degrees
        beta[:] = 2*q*np.cos(np.deg2rad(lat_deg[:]))/a
        f[:]=beta[:]*Y[0,:]
            
    return scen, x, y, X, Y, u, v, h, Fu, f, div, vor

def leapfrog(Nx, Ny, Nt, u, v, h, g, dx, dy, dt, H, Fu, f, c, div, vor):
    """_summary_

    Args:
        Nx (int): _description_
        Ny (int): _description_
        Nt (int): _description_
        u (array): _description_
        v (array): _description_
        h (array): _description_
        g (float): _description_
        dx (int): _description_
        dy (int): _description_
        dt (int): _description_
        H (int): Mean height
        Fu (array): Constant source at zonal wind.
        f (array): _description_
        c (float): _description_
        div (array): Divergence
        vor (array): Vorticity

    Returns:
        _type_: _description_
    """
    for n in range(1, Nt):

        if n == 1: # Euler-forward scheme
            for i in range(1, Nx-1):
                for j in range(1, Ny-1):
                    u[n,i+1,j]= u[n-1,i+1,j] + dt*(-g*(h[n-1,i+2,j] - h[n-1,i+1,j])/dx + \
                                f[j]/4*(v[n-1,i,j]+v[n-1,i+1,j]+v[n-1,i+1,j-1]+v[n-1,i,j-1]) + Fu[n,i,j]) 
                    v[n,i,j+1]= v[n-1,i,j+1] + dt*(-g*(h[n-1,i,j+2] - h[n-1,i,j+1])/dy - \
                                f[j]/4*(u[n-1,i,j]+u[n-1,i,j+1]+u[n-1,i-1,j+1]+u[n-1,i-1,j]))
                    h[n,i,j]= h[n-1,i,j] - dt*H*((u[n-1,i,j]-u[n-1,i-1,j])/dx + (v[n-1,i,j]-v[n-1,i,j-1])/dy)
                    
        else: # Leap-frog scheme
            for i in range(1, Nx-1):
                for j in range(1, Ny-1):
                    u[n+1,i+1,j]= u[n-1,i+1,j] + 2*dt*(-g*(h[n,i+2,j] - h[n,i+1,j])/dx + \
                                f[j]/4*(v[n,i,j]+v[n,i+1,j]+v[n,i+1,j-1]+v[n,i,j-1]) + Fu[n,i,j])
                    v[n+1,i,j+1]= v[n-1,i,j+1] + 2*dt*(-g*(h[n,i,j+2] - h[n,i,j+1])/dy - \
                                f[j]/4*(u[n,i,j]+u[n,i,j+1]+u[n,i-1,j+1]+u[n,i-1,j]))
                    h[n+1,i,j]= h[n-1,i,j] - 2*dt*H*((u[n,i,j]-u[n,i-1,j])/dx + (v[n,i,j]-v[n,i,j-1])/dy)
        
        # Radiational as boundary condition (west, north, south)
        # ------------------------------------------------------
        # west
        u[n+1,0,:]= u[n,0,:] + dt*(f[:]*(v[n,0,1:]+v[n,0,:-1])/2 + c*(u[n,1,:]-u[n,0,:])/dx)
        # north
        v[n+1,:,-1]= v[n,:,-1] - dt*(f[-1]*(u[n,1:,-1] + u[n,:-1,-1])/2 + c*(v[n,:,-2]-v[n,:,-1]/dy))
        # south
        v[n+1,:,0]= v[n,:,0] - dt*(f[-1]*(u[n,1:,0]+ u[n,:-1,0])/2 + c*(v[n,:,1]-v[n,:,0]/dy))
        
        # Fixed boundary condition at east
        u[n,-2::,:] =  0
        h[n,-2::,:] =  5
    
    # Divergence (du/dx + dv/dy)
    # -------------------------
    div[:,:,:]=(u[:,1:, :] - u[:,:-1, :])/dy + (v[:,:,1:] - v[:,:,:-1])/dx
    
    # Vorticity (dv/dx - du/dy)
    # ---------
    vor[:,:,:] = (v[:,:,1:] - v[:,:,:-1])/dx - (u[:,1:,:] - u[:,:-1,:])/dy
        
    return u, v, h, div, vor
