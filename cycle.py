import numpy as np
from parameters import *


def density (parts):
    """
    wheiting density on grids from particle position data

    Args:
        parts (np.array): particle position array

    Returns:
        rho (np.array): dencity on grids
    """
    # initialize rho: filling by zero 
    rho = [0.0 for i in range(NODES)]
    n_parts = len(parts)
   
    # first-order CIC weighting: qj=qc*(1-(xi-Xj)/dx)
    for k in range(n_parts):
        xi = (parts[k].x / dX)
        Xj = np.floor(xi)
        d = xi - Xj

        cur = int(Xj)
        nxt = (cur + 1)

        #print ("\n",cur,k)
        rho[cur] += parts[k].q * (1.0 - d)
        rho[nxt] += parts[k].q * d


    rho[NODES - 1] += rho[0]
    rho[0] = rho[NODES - 1]

    # convert list to array
    rho = np.array(rho)

    #normalization to dx
    # rho /= dX
    if VERBOSE: 
        print("\n wheiting dencity on grids from particle position (rho) by dencity funcion done! \n")
    return rho


def sor_solver (rho):
    """
    Solving Poisson Eq by Successive Over-Relaxation (SOR) Method  

    Args:
        rho (np.array): density of particles on grids

    Returns:
        phi (np.array): potential on grids (answer of SOR solver)
    """
    #  relaxation factor
    omega = 2.0 / (1 + 2 * np.pi / NODES)

    # initial phi values
    phi = np.array([0.0 for i in range(NODES)])

    # right hand side
    rhs = -np.copy(rho) * dX**2 / eps0

    # solver
    for k in range(SOR_MAX_ITR):
        
        if k == SOR_MAX_ITR - 1:
            print( "\n Ops! SOR (likely) diverges :( \n")
        
        phinew = np.copy(phi)
        for i in range(NODES - 1):
            nxt = (i + 1) if (i < NODES - 2) else 0
            prv = (i - 1) if (i > 0) else (NODES - 2)

            phinew[i] = (1 - omega) * phi[i] + (omega / -2.0) * (rhs[i] - phinew[prv] - phi[nxt])

        # checking convergence
        if k % 32 == 0:
            err = np.max(np.abs(phinew - phi))
            if err < SOR_ERR:
                phi = phinew
                if VERBOSE:
                    print(f"Poisson solved by SOR in {k} steps by {err} error ")
                break
        phi = np.copy(phinew)

    phi[NODES - 1] = phi[0]
    return phi


def fieldOnNodes (phi):
    """
    Calculating field from potential on grids

    Args:
        phi (np.array): potential on grids

    Returns:
        field (np.array): E field on grids
    """
    efield = np.array([0.0 for i in range(NODES)])
    
    for i in range(NODES):
        nxt = (i + 1) if (i < NODES - 1) else 0
        prv = (i - 1) if (i > 0) else (NODES - 1)

        efield[i] = (phi[prv] - phi[nxt]) / (2 * dX)

    if VERBOSE:
        print("\n Calculating field from potential on grids (efild) done!\n")
    return efield


def fieldOnParticles (field, parts):
    """
    Interpolating field from grids to particles

    Args:
        field (np.array): E field on grids
        parts (np.array): particles on grids

    Returns:
        efield (np.array): E field on partilce
    """
    n_parts = len(parts)
    efield = np.array([0.0 for i in range(n_parts)])
    
    # first-order CIC backwards-weighting:
    # E(xi) = Ej * (Xj+1 - xi)/dx + Ej+1 * (xi-Xj)/dx  
    for k in range(n_parts):
        if parts[k].mv:
            xi = (parts[k].x / dX)
            j = int(np.floor(xi))

            nxt = (j + 1) if (j + 1) < NODES else 0

            efield[k] = (nxt - xi) * field[j] + (xi - j) * field[nxt]
    if VERBOSE:
        print("\n calculating E field on particles done! \n")
    return efield


def rewind (direction, field, parts):
    """
    First time rewind velocity by dT/2 forward or backward

    Args:
        direction (np.array): rewind direction. +1: forward, -1: backward
        field (np.array): E field on particles E(xi)
        parts (np.array): particles

    Returns:
        parts (np.array): update velocity
    """
    n_parts = len(parts)
    for k in range(n_parts):
        # updating velocity
        if parts[k].mv:
            parts[k].v += direction * field[k] * parts[k].qm * dT / 2.0
    return parts


def moveParticles (field, parts):
    """
    motion eq. integration for particles

    Args:
        field (np.array): E field on particles E(xi)
        parts (np.array): particles

    Returns:
        psrts (np.array): particles with new position and velocity
    """
    n_parts = len(parts)
    for k in range(n_parts):
        if parts[k].mv:
            # updating velocity
            parts[k].v += field[k] * parts[k].qm * dT

            # updating position
            parts[k].x += parts[k].v * dT

            while parts[k].x < 0:
                parts[k].x += SIZE
            while parts[k].x >= SIZE:
                parts[k].x -= SIZE
    
    if VERBOSE:
        print("\n moving of particles done! \n")
    return parts