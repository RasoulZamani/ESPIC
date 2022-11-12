import numpy as np
from constants import *


def density (parts, cells_num, dx):
    """
    wheiting density on grids from particle position data

    Args:
        parts (np.array): particle info array
        cells_num (int) : number of cells
        dx      (float) : cell length

    Returns:
        rho (np.array): dencity on grids
    """
    # initialize rho: filling by zero 
    rho = [0.0 for i in range(cells_num + 1)]
    n_parts = len(parts)
   
    # first-order CIC weighting: qj=qc*(1-(xi-Xj)/dx)
    for k in range(n_parts):
        xi = (parts[k].x / dx)
        Xj = np.floor(xi)
        d = xi - Xj

        cur = int(Xj)
        nxt = (cur + 1)

        #print ("\n",cur,k)
        rho[cur] += parts[k].q * (1.0 - d)
        rho[nxt] += parts[k].q * d


    rho[cells_num + 1 - 1] += rho[0]
    rho[0] = rho[cells_num + 1 - 1]

    # convert list to array
    rho = np.array(rho)

    # for chek stuf:
    #print("rho is:\n", rho)
    
    #normalization to dx
    # rho /= dX
    return rho


def sor_solver (rho, cells_num, dx):
    """
    Solving Poisson Eq by Successive Over-Relaxation (SOR) Method  

    Args:
        rho (np.array): density of particles on grids
        cells_num(int): number of cells
        dx     (float): cell length

    Returns:
        phi (np.array): potential on grids (answer of SOR solver)
    """
    #  relaxation factor must be between 0,2
    omega = 2.0 / (1 + 2 * np.pi / (cells_num + 1))

    # initial phi values
    phi = np.array([0.0 for i in range(cells_num + 1)])

    # right hand side
    rhs = -np.copy(rho) * dx**2 / EPS0

    # solver
    for k in range(SOR_MAX_ITR):
        
        if k == SOR_MAX_ITR - 1:
            print( "\n Ops! SOR (likely) diverges :( \n")
        
        phinew = np.copy(phi)
        for i in range(cells_num ):
            nxt = (i + 1) if (i < cells_num - 1 ) else 0
            prv = (i - 1) if (i > 0) else (cells_num -1 )

            phinew[i] = (1 - omega) * phi[i] + (omega / -2.0) * (rhs[i] - phinew[prv] - phi[nxt])

        # checking convergence
        if k % 32 == 0:
            err = np.max(np.abs(phinew - phi))
            if err < SOR_ERR:
                phi = phinew
                break
        phi = np.copy(phinew)

    phi[cells_num + 1 - 1] = phi[0]
    return phi


def fieldOnNodes (phi,cells_num, dx):
    """
    Calculating field from potential on grids

    Args:
        phi (np.array)  : potential on grids
        cells_num (int) : number of cells
        dx      (float) : cell length

    Returns:
        field (np.array): E field on grids
    """
    efield = np.array([0.0 for i in range(cells_num + 1)])
    
    for i in range(cells_num + 1):
        nxt = (i + 1) if (i < cells_num + 1 - 1) else 0
        prv = (i - 1) if (i > 0) else (cells_num + 1 - 1)

        efield[i] = (phi[prv] - phi[nxt]) / (2 * dx)

    return efield


def fieldOnParticles (field, parts, cells_num, dx):
    """
    Interpolating field from grids to particles

    Args:
        field (np.array): E field on grids
        parts (np.array): particles on grids
        cells_num (int) : number of cells
        dx      (float) : cell length

    Returns:
        efield (np.array): E field on partilce
    """
    n_parts = len(parts)
    efield = np.array([0.0 for i in range(n_parts)])
    
    # first-order CIC backwards-weighting:
    # E(xi) = Ej * (Xj+1 - xi)/dx + Ej+1 * (xi-Xj)/dx  
    for k in range(n_parts):
        if parts[k].mv:
            xi = (parts[k].x / dx)
            j = int(np.floor(xi))

            nxt = (j + 1) if (j + 1) < (cells_num + 1) else 0

            efield[k] = (nxt - xi) * field[j] + (xi - j) * field[nxt]
    return efield


def rewind (direction, field, parts, dt):
    """
    First time rewind velocity by dT/2 forward or backward

    Args:
        direction (np.array): rewind direction. +1: forward, -1: backward
        field (np.array): E field on particles E(xi)
        parts (np.array): particles info array
        dt     (float)  : time step 

    Returns:
        parts (np.array): particles info array with updated velocity
    """
    n_parts = len(parts)
    for k in range(n_parts):
        # updating velocity
        if parts[k].mv:
            parts[k].v += direction * field[k] * parts[k].qm * dt / 2.0
    return parts


def moveParticles (field, parts, size, dt):
    """
    motion eq. integration for particles

    Args:
        field (np.array): E field on particles E(xi)
        parts (np.array): particles
        size  (float)   : size of environment 
        dt    (float)   : time step

    Returns:
        psrts (np.array): particles with updated position and velocity
    """
    n_parts = len(parts)
    for k in range(n_parts):
        if parts[k].mv:
            # updating velocity
            parts[k].v += field[k] * parts[k].qm * dt

            # updating position
            parts[k].x += parts[k].v * dt

            while parts[k].x < 0:
                parts[k].x += size
            while parts[k].x >= size:
                parts[k].x -= size
    
    return parts