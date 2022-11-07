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

    beta = np.zeros((Ny+1))
    div  = np.zeros((Nt+1,Nx+1,Ny+1))
    vor  = np.zeros((Nt+1,Nx+1,Ny+1))
     
    # Allocating arrays and initial conditions as Arakawa C-grid
    # ----------------------------------------------------------
    u  = np.zeros((Nt+1, Nx+2, Ny+1))  # zonal wind
    v  = np.zeros((Nt+1, Nx+1, Ny+2))  # meridional wind
    h  = np.zeros((Nt+1, Nx+1, Ny+1))  # height variation by time
    Fu = np.zeros((Nt+1, Nx+1, Ny+1))  # Zonal force like ENSO
    
    fonte_u = -np.exp(-X**2/(Nrx*dx)**2 - Y**2/(Nry*dy)**2)
    Fu[:,:,:] =  0 #fonte_u
    
    # Initial conditions for u, v, & h
    # --------------------------------
    u[0, :, :] = 0.
    v[0, :, :] = 0.
    h[0, 20, 20] = 30 #fonte_u*-20

    q = 2*np.pi/86400      # angular velocity (2pi/(24*3600 s)) [1/s]
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
            
    return scen, x, y, X, Y, u, v, h, Fu, f, beta, div, vor

def leapfrog(scen, Nx, Ny, Nt, u, v, h, g, dx, dy, dt, H, Fu, f, beta, y, c, div, vor):
    """Leapfrog scheme to resolve shallow water equation 2D.

    Args:
        Nx (int): n° points in x axes
        Ny (int): n° points in y axes
        Nt (int): n° time steps
        u (array): zonal wind vectors
        v (array): meridional wind vectors
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

    Returns:
        arrays: u, v, h, div, vor
    """
    for n in range(1, Nt):
        # Radiational as boundary condition (west, north, south, east)
        # ------------------------------------------------------
        # west (esquerdo) -> du/dt - fv - c*du/dx = 0
        #u[n+1,0,:] = u[n,0,:] + dt*(f[:]*(v[n,0,1:]+v[n,0,:-1])/2 + c*(u[n,1,:]-u[n,0,:])/dx)
        u[n+1,0,:] = 0
        
        # north -> dv/dt + fu + c*dv/dy = 0
        #v[n+1,:,-1] = v[n,:,-1] - dt*(f[-1]*(u[n,1:,-1] + u[n,:-1,-1])/2 + c*(v[n,:,-1]-v[n,:,-2])/dy)
        v[n+1,:,-1] = 0
        
        # south -> dv/dt + fu - c*dv/dy = 0
        #v[n+1,:,0] = v[n,:,0] - dt*(f[-1]*(u[n,1:,0]+ u[n,:-1,0])/2 - c*(v[n,:,1]- v[n,:,0])/dy)
        v[n+1,:,0] = 0
        
        # East -> du/dt - fv + c*du/dx = 0
        #u[n+1,-1,:] = u[n,-1,:] + dt*(f[:]*(v[n,-1,1:]+v[n,-1,:-1])/2 - c*(u[n,-1,:]-u[n,-2,:])/dx)
        u[n+1,-1,:] =  0
        
        if (scen == 'scen1') | (scen == 'scen2'):
            if n == 1: # Euler-forward scheme
                for i in range(1, Nx):
                    for j in range(1, Ny):
                        u[1,i,j]= u[0,i,j] + dt*(-g*(h[0,i+1,j] - h[0,i,j])/dx + \
                                    f[j]/4*(v[0,i,j] + v[0,i+1,j] + v[0,i+1,j-1] + v[0,i,j-1]) + Fu[n,i,j]) # ok
                        
                        v[1,i,j]= v[0,i,j] + dt*(-g*(h[0,i,j+1] - h[0,i,j])/dy - \
                                    f[j]/4*(u[0,i,j] + u[0,i,j+1] + u[0,i-1,j+1] + u[0,i-1,j])) # ok
                        
                        h[1,i,j]= h[0,i,j] - dt*H*( (u[0,i,j] - u[0,i-1,j])/dx + (v[0,i,j] - v[0,i,j-1])/dy ) #ok
                        
            else: # Leap-frog scheme
                for i in range(1, Nx):
                    for j in range(1, Ny):
                        u[n+1,i,j] = u[n-1,i,j] + 2*dt*(-g*(h[n,i+1,j] - h[n,i,j])/dx + \
                                     f[j]/4*(v[n,i,j] + v[n,i+1,j] + v[n,i+1,j-1] + v[n,i,j-1]) + Fu[n,i,j]) #ok
                        v[n+1,i,j] = v[n-1,i,j] + 2*dt*(-g*(h[n,i,j+1] - h[n,i,j])/dy - \
                                     f[j]/4*(u[n,i,j] + u[n,i,j+1] + u[n,i-1,j+1] + u[n,i-1,j]))  # ok
                        h[n+1,i,j] = h[n-1,i,j] - 2*dt*H*((u[n,i,j] - u[n,i-1,j])/dx + (v[n,i,j] - v[n,i,j-1])/dy) # ok
                               
        elif scen == 'scen3': # Beta 
            if n == 1: # Euler-forward scheme
                for i in range(1, Nx):
                    for j in range(1, Ny):
                        u[n,i,j]= u[n-1,i,j] + dt*(-g*(h[n-1,i+1,j] - h[n-1,i,j])/dx + \
                                    beta[j]/4*(y[j]*v[n-1,i,j]*(v[n-1,i,j]+v[n-1,i+1,j]) + \
                                        y[j-1]*v[n-1,i+1,j-1]*(v[n-1,i+1,j-1]+v[n-1,i,j-1])) + Fu[n,i,j]) 
                        v[n,i,j]= v[n-1,i,j] + dt*(-g*(h[n-1,i,j+1] - h[n-1,i,j])/dy - \
                                    beta[j]/4*(y[j]*u[n-1,i,j]*h[n-1,i,j]*(u[n-1,i,j]+u[n-1,i-1,j]) + \
                                        y[j+1]*u[n-1,i,j+1]*h[n-1,i,j+1]*(u[n-1,i,j+1]+u[n-1,i-1,j+1])))
                        h[n,i,j]= h[n-1,i,j] - dt*H*((u[n-1,i,j]-u[n-1,i-1,j])/dx + (v[n-1,i,j]-v[n-1,i,j-1])/dy) 
                        
            else: # Leap-frog scheme
                for i in range(1, Nx):
                    for j in range(1, Ny):
                        u[n+1,i,j]= u[n-1,i,j] + 2*dt*(-g*(h[n,i+1,j] - h[n,i,j])/dx + \
                                      beta[j]/4*(y[j]*v[n,i,j]*(v[n,i,j]+v[n,i+1,j]) + \
                                        y[j-1]*v[n,i+1,j-1]*(v[n,i+1,j-1]+v[n,i,j-1])) + Fu[n,i,j])
                        v[n+1,i,j]= v[n-1,i,j] + 2*dt*(-g*(h[n,i,j+1] - h[n,i,j])/dy - \
                                    beta[j]/4*(y[j]*u[n,i,j]*h[n,i,j]*(u[n,i,j]+u[n,i-1,j]) + \
                                        y[j+1]*u[n,i,j+1]*h[n,i,j+1]*(u[n,i,j+1]+u[n,i-1,j+1])))
                        h[n+1,i,j]= h[n-1,i,j] - 2*dt*H*((u[n,i,j]-u[n,i-1,j])/dx + (v[n,i,j]-v[n,i,j-1])/dy)   
        
    # Divergence (du/dx + dv/dy)
    # -------------------------
    div[:,:,:]=(u[:,1:, :] - u[:,:-1, :])/dx + (v[:,:,1:] - v[:,:,:-1])/dy
    
    # Vorticity (dv/dx - du/dy)
    # ---------
    vor[:,:,:] = (v[:,:,1:] - v[:,:,:-1])/dx - (u[:,1:,:] - u[:,:-1,:])/dy
        
    return u, v, h, div, vor
