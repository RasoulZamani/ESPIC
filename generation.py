from parameters import *
import numpy as np

class Particle:

    def __init__(self,
                 pos, # position of particles
                 vel,       # velocity of particles
                 omega_p,    # plasma freq
                 QoverM ,    # q/m
                 move,      # moving(True) or not(False)
                 parts_num,       # number of particles
                 size=SIZE):
        self.x = pos
        self.v = vel
        # from wp^2 = n*q^2 / m*eps  and n=N/L we drecive q:
        self.q = omega_p**2 * (1 / QoverM) * eps0 * (size /parts_num)
        self.qm = QoverM
        self.mv = move


def two_stream (omega_p, parts_num, size, # for Particle class
                  mode= MODE , ampl = AMPL, v0 = V0 ):
    
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

    sep = size / parts_num
    for i in range (parts_num):
        x0 = (i + 0.5) * sep
        particles.append(Particle (pos=x0,vel= 0.0,omega_p=Omega_pp,
                                   QoverM=QM_p, move=False, parts_num=parts_num))

    return particles