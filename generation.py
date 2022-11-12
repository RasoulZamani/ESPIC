""" 
this file (generation) contains class for particles
and function for creating two_stream and e_beam. 
"""
# import libs: _______________________________________________
from constants import *
import numpy as np
# ____________________________________________________________
class Particle:

    def __init__(self,
                 pos, # position of particles
                 vel,       # velocity of particles
                 omega_p,    # plasma freq
                 QoverM ,    # q/m
                 move,      # moving(True) or not(False)
                 parts_num,       # number of particles
                 size):
        self.x = pos
        self.v = vel
        # from wp^2 = n*q^2 / m*eps  and n=N/L we drecive q:
        self.q = omega_p**2 * (1 / QoverM) * EPS0 * (size /parts_num)
        self.qm = QoverM
        self.mv = move

# ____________________________________________________________
def two_stream (omega_p, parts_num, size, # for Particle class
                  mode, ampl, # perturbation
                  v0 # velocity of streams 
                  ):
    """ generating parts_num number of particle object informing two streams
        with omega_p in lenthg of size and perturbate it.
    
    Args:
        omega_p   (int)      : plasma freq: ~ 9 sqrt(n)
        parts_num (int)      : number of particles will be generated
        mode(int),ampl(float): for purturbating sreams
        v0 (int)             : velocity of streams
         
    Returns:
        particles (array obj): contains particles information like x, v , ... 
    """
    
    particles = []
    sep = 1.0 * size / (parts_num / 2) #2dX
    for i in range(int(parts_num / 2)):
        # unperturbed position
        x0 = (i + 0.5) * sep

        # perturbation
        theta = 2 * np.pi * mode * x0 / size
        dx = ampl * np.cos(theta)
        x1 = x0 + dx
        x2 = x0 - dx

        # periodic boundaries
        if x1 < 0:
            x1 += size
        if x2 < 0:
            x2 += size
        if x1 >= size:
            x1 -= size
        if x2 >= size:
            x2 -= size

        # add to particles
        particles.append(Particle (pos=x1, vel=-1.0*v0, omega_p=omega_p,
            QoverM = QM_e, move=True, parts_num=parts_num, size=size))
        particles.append(Particle (pos=x2, vel= 1.0*v0, omega_p=omega_p,
            QoverM = QM_e, move=True, parts_num=parts_num, size=size))

    # adding non-movable ions in background __________________
    sep = size / parts_num
    for i in range (parts_num):
        x0 = (i + 0.5) * sep
        particles.append(Particle (pos=x0,vel= 0.0,omega_p=omega_p*(m_e/m_p)**(0.5),
            QoverM=QM_p, move=False, parts_num=parts_num, size=size))

    return particles
# ____________________________________________________________
def e_beam (omega_p, parts_num, size, # for Particle class
                  mode, ampl, # perturbation
                  v0 # velocity of beam
                  ):
    """ generating parts_num number of particle object and form an
    electron beam with omega_p, with velocity v0, in lenthg of size. 
    
    Args:
        omega_p   (int)      : plasma freq: ~ 9 sqrt(n)
        parts_num (int)      : number of particles will be generated
        mode(int),ampl(float): for purturbating sreams
        v0 (int)             : velocity of beam
         
    Returns:
        particles (array obj): contains particles information like x, v , ... 
    """
    
    particles = []
    sep = 1.0 * size / (parts_num ) #dX
    for i in range(int(parts_num)):
        # unperturbed position
        x = (i + 0.5) * sep

        # perturbation
        if mode != None:
            theta = 2 * np.pi * mode * x / size
            dx = ampl * np.cos(theta)
            x = x + dx

        # periodic boundaries
        if x < 0:
            x += size
        if x >= size:
            x -= size

        # add to particles
        particles.append(Particle (pos = x, vel = v0, omega_p = omega_p,
            QoverM = QM_e, move = True, parts_num = parts_num, size = size))

    # adding non-movable ions in background __________________
    sep = size / parts_num
    for i in range (parts_num):
        x0 = (i + 0.5) * sep
        particles.append(Particle (pos=x0, vel=0.0, omega_p = omega_p*(m_e/m_p)**(0.5),
            QoverM = QM_p, move = False, parts_num = parts_num, size = size))

    return particles