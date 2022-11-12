"""
this file (e_beam_virtual_cathod) is for simulating virtual cathod formation!
you can alter input file like L,v0, NP, ...
"""
from constants import *

from generation import e_beam

from pic_simulator import Sim


# inpyt params: ______________________________________________
PARTS_GEN = e_beam
Omega_p = 1.0e10             # non normalized plasma frequency ~ 9 sqrt(n0)

CELLS_NUM = 64             # number of cells
#node_num = cells_num + 1       # number of node points
PARTS_NUM = CELLS_NUM * 20
SIZE = 1 * (Omega_p/C)             # size of the system in meter

STEPS = 1000
PLT_STP = int(STEPS / 10) # plotting periode

dX = SIZE / CELLS_NUM       # distance between nodes (spatial step)
    
V0 = 0.4 * C
# verbose get more info of proccess
VERBOSE = False #True

# address of directory for saving plots
RESULT_DIR = "./results"
# Plasma and Beam parameter
Omega_pp= Omega_p * ( (m_e/m_p)**(0.5) )               # plasma freq for ions

dT =  1/(200*Omega_p)                # dt timestep (should be: omegaP * dt < 0.3)

# PERTURBATION
MODE = None # None for no perturbation
AMPL = None #SIZE/5

# instantiate Sim() class for simulation env:_________________
sim = Sim(size = SIZE, cells_num = CELLS_NUM, time_step = dT,
        parts_num = PARTS_NUM,
        dir = RESULT_DIR, verbose = VERBOSE)
# running PIC main loop:
sim.run(parts_gen = PARTS_GEN, steps = STEPS, omega_p = Omega_p,
        mode = MODE , ampl = AMPL, v0 = V0,
        plot = True, plt_stp = PLT_STP)
     