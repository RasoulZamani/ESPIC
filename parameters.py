""" 
This file contains constanses and parameters needed for code
also name of used model shoud be choice here
"""
# Constansts__________________________________________________
EPS0 = 8.85e-12
m_e = 9.1091 * 1e-31 # Mass of an electron Kg
m_p = 1.6726 * 1e-27
c = 299790000 # Speed of light meters/second
q = 1.602 * 1e-19  # absulote Charge of an electron coulombs
QM_e = -1 * q/m_e
QM_p = q/ m_p 
# GRIDS ______________________________________________________
CELLS = 100             # number of cells
NODES = CELLS + 1       # number of node points

SIZE = .4             # size of the system in meter

STEPS = 1000
PLT_STP = int(STEPS / 10) # plotting periode

NAME = 'two-stream'
dX = SIZE / CELLS       # distance between nodes (spatial step)
dT = 1e-5                # timestep (should be: omegaP * dt < 0.3)

# Plasma and Beam parameter
omega_p = 1e8             # non normalized plasma frequency ~ 9 sqrt(n0)
n0      = 1e12
eps0 = EPS0              # non normalized vacuum permittivity
keV  = 2                 # beam energy in keV
# getting v in m/c from beam energy
gamma = 1 + keV/511
v_c = (1 - 1/(gamma**2) ) ** (0.5)
V0   = v_c * c


# SOR (pisson solver) parameters
SOR_ERR = 1e-5
SOR_MAX_ITR = int(1e4)
# PARTICLES single
NPpC = 10               # number of particles per species per cell
NP = NPpC * CELLS       # number of particles (of particular specie) (used to determine charge!)

# PERTURBATION
MODE = 2
AMPL = SIZE/5

# verbose get more info of proccess
VERBOSE = False #True

# address of directory for saving plots
RESULT_DIR = "./results"