import os
import numpy as np
from parameters import *
import matplotlib.pyplot as plt

import generation as gen
from cycle import density, sor_solver, fieldOnNodes, \
fieldOnParticles, rewind, moveParticles
import utils

# init _______________________________________________________
# creating directory for savuing results
dir = os.path.join(RESULT_DIR)
os.makedirs(dir, exist_ok=True)
# clear possible files from previouse run:
files = os.listdir(dir) # list of existing file in dir
if len(files)>0:
    for f in files:
       os.remove(os.path.join(dir,f))
       
if VERBOSE:
    print(f"\n directory:  {dir}  created for saving result\n")
# particle generation ________________________________________
#particles = gen.ebeam(v0)
particles = gen.twoStream()
if VERBOSE:
    # ..... TODO ... : print other parameters too
    print( " particles in t=0 generated with ...")
    print(f"\nnumber of cells is: {CELLS}, length of env: {SIZE}m",
          f"\nlength of each grid(dX): {dX}, time step: {dT}, number of steps:{STEPS}",
          f"\nplasma freq for e:{omega_p}, beam energy: {keV} keV equal to {V0} m/s")
# PROGRESS BAR _______________________________________________
stp = 0
utils.printProgress(stp, STEPS, prefix = 'Progress:',
                  suffix = 'Complete', barLength = 50)

# main loop of PIC ___________________________________________
for step in range(STEPS):
    if VERBOSE:
        print(f"\n we are in step {step} \n")
    rho = density(particles)
    PHI = sor_solver(rho)
    EFIELDn = fieldOnNodes(PHI)
    EFIELD = fieldOnParticles(EFIELDn, particles)

    if step == 0:
        particles = rewind(-1, EFIELD, particles)

    particles = moveParticles(EFIELD, particles)

    # write to file
    #  ..... TODO ....
    

    # PROGRESS BAR
    stp += 1
    utils.printProgress(stp, STEPS, prefix = 'Progress:',
                      suffix = 'Complete', barLength = 50)
    
    # plotting _______________________________________________

    if step % PLT_STP == 0:
        n_parts = len(particles)
        xi = []
        vi = []
        for i in range(n_parts):
            if particles[i].mv:
                xi.append(particles[i].x)
                vi.append(particles[i].v)
        
        fig = plt.figure()      
        plt.scatter(x=xi, y=vi)
        plt.xlabel(" x [position] ")
        plt.ylabel(" v (velocity) ")
        plt.title(f" phase space (x-v) of particles in step {step} ")
        fig.savefig(r"results/fig_"+f"step_{step}")
        #plt.show()
        
        # ...... TODO ....... 
        # make t more pretty later!
