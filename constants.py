""" 
This file contains constanses and parameters needed for code
also name of used model shoud be choice here
"""

# Constansts__________________________________________________

EPS0 = 8.85e-12
m_e = 9.1091 * 1e-31 # Mass of an electron Kg
m_p = 1.6726 * 1e-27
C = 299790000 # Speed of light meters/second
q = 1.602 * 1e-19  # absulote Charge of an electron coulombs
QM_e = -1 * q/m_e
QM_p = q/ m_p 
# SOR (pisson solver) parameters
SOR_ERR = 1e-5
SOR_MAX_ITR = int(1e4)


"""
keV  = 2                 # beam energy in keV
# getting v in m/c from beam energy
gamma = 1 + keV/511
v_c = (1 - 1/(gamma**2) ) ** (0.5)
V0   = v_c * C

"""



