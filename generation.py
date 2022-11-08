from parameters import *
import numpy as np

class Particle:

    def __init__(self, pos, # position of particles
                 vel,       # velocity of particles
                 omega_p,    # plasma freq
                 QoverM ,    # q/m
                 move,      # moving(True) or not(False)
                 num,       # number of particles
                 ):
        self.x = pos
        self.v = vel
        # from wp^2 = n*q^2 / m*eps  and n=N/L we drecive q:
        self.q = omega_p**2 * (1 / QoverM) * eps0 * (SIZE / num)
        self.qm = QoverM
        self.mv = move

def twoStream ():
    PARTS = []
    sep = 1.0 * SIZE / (NP / 2) #2dX
    for i in range(int(NP / 2)):
        # unperturbed position
        x0 = (i + 0.5) * sep

        # perturbation
        theta = 2 * np.pi * MODE * x0 / SIZE
        dx = AMPL * np.cos(theta)
        x1 = x0 + dx
        x2 = x0 - dx

        # periodic boundaries
        if x1 < 0:
            x1 += SIZE
        if x2 < 0:
            x2 += SIZE
        if x1 >= SIZE:
            x1 -= SIZE
        if x2 >= SIZE:
            x2 -= SIZE

        # add to PARTS
        PARTS.append(Particle (x1, -1.0 * V0, omega_p, QM_e, True, NP))
        PARTS.append(Particle (x2, 1.0 * V0, omega_p, QM_e, True, NP))

    sep = SIZE / NP
    for i in range (NP):
        x0 = (i + 0.5) * sep
        PARTS.append(Particle (x0, 0.0, omega_p*np.sqrt(m_e/m_p), QM_p, False, NP))

    return PARTS

def ebeam (v0):
    sep = 1.0 * SIZE / NP
    PARTS = []
    for i in range(NP):
        x0 = (i + 0.5) * sep
        x0 += np.random.uniform(-dX, dX)

        if x0 < 0:
            x0 += SIZE
        if x0 >= SIZE:
            x0 -= SIZE
        PARTS.append(Particle (x0, v0, 1.0, -1.0, True, NP))
    for i in range(NP):
        x0 = (i + 0.5) * sep
        PARTS.append(Particle (x0, 0.0, 1.0, 1.0, False, NP))
    return PARTS
