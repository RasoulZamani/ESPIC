""" 
This file contains constanses needed for code
"""
# physical constansts_________________________________________

EPS0 = 8.85e-12
m_e  = 9.1091 * 1e-31        # Mass of electron Kg
m_p  = 1.6726 * 1e-27        # Mass of proton Kg
C    = 299790000             # Speed of light meters/second
q    = 1.602 * 1e-19         # absulote Charge of an electron coulombs
QM_e = -1 * q/m_e
QM_p = q/ m_p 

# SOR (Poisson solver) parameters ____________________________
SOR_ERR     = 1e-4
SOR_MAX_ITR = int(1e4)


"""
# converting keV to m/s:
keV  = 2                 # beam energy in keV
# getting v in m/c from beam energy
gamma = 1 + keV/511
v_c = (1 - 1/(gamma**2) ) ** (0.5)
V0   = v_c * C
"""
