"""
this file (two_stream_instability) is for simulating two stream instability!
you can alter input file like L,v0, NP , ... but remember criterion for this
instability: L*Wp > 2 * pi / sqrt(2)
"""
from constants import *

from generation import two_stream 

from pic_simulator import Sim


# inpyt params: ______________________________________________
PARTS_GEN = two_stream
CELLS_NUM = 200             # number of cells
#node_num = cells_num + 1       # number of node points
PARTS_NUM = CELLS_NUM * 10
SIZE = .2             # size of the system in meter

STEPS = 1000
PLT_STP = int(STEPS / 40) # plotting periode

dX = SIZE / CELLS_NUM       # distance between nodes (spatial step)
dT = 1e-9                # dt timestep (should be: omegaP * dt < 0.3)
    
V0 = 1e5
# verbose get more info of proccess
VERBOSE = False #True

# address of directory for saving plots
RESULT_DIR = "./results"
# Plasma and Beam parameter
Omega_p = 1e8             # non normalized plasma frequency ~ 9 sqrt(n0)
Omega_pp= Omega_p * ( (m_e/m_p)**(0.5) )               # plasma freq for ions

# PERTURBATION
MODE = 1
AMPL = 0.1#SIZE/5

# instantiate Sim() class for simulation env:_________________
sim = Sim(size = SIZE, cells_num = CELLS_NUM, time_step = dT,
        parts_num = PARTS_NUM,
        dir = RESULT_DIR, verbose = VERBOSE)
# running PIC main loop:
sim.run(parts_gen = PARTS_GEN, steps = STEPS, omega_p = Omega_p,
        mode = MODE , ampl = AMPL, v0 = V0,
        plot = True, plt_stp = PLT_STP)
     