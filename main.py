import os
import numpy as np
from parameters import *
import matplotlib.pyplot as plt

from generation import two_stream
from cycle import density, sor_solver, fieldOnNodes, \
fieldOnParticles, rewind, moveParticles
import utils


class Sim():
    """class for simulation of pic
    """
    def __init__(self,size, parts_num, dir, verbose = False):
        self.size = size
        self.parts_num = parts_num
        self.dir = dir
        self.verbose = verbose
        
        # creating directory for savuing results
        dir_path = os.path.join(self.dir)
        os.makedirs(dir_path, exist_ok=True)
        # clear possible files from previouse run:
        files = os.listdir(dir_path) # list of existing file in dir
        if len(files)>0:
            for f in files:
                os.remove(os.path.join(dir_path,f))
       
        if self.verbose:
            print(f"\n directory:  {dir}  created for saving result\n")
    

    def run(self, steps=STEPS, omega_p=Omega_p,
            mode= MODE , ampl = AMPL, v0 = V0,
            plot=True, plt_stp = PLT_STP):
        
        # generating particle:
        self.particles = two_stream(omega_p=omega_p, parts_num = self.parts_num, size = self.size, # for Particle class
                  mode= mode , ampl = ampl, v0 = v0 )
        if self.verbose:
            # ..... TODO ... : print other parameters too
            print( " particles in t=0 generated with ...")
            print(f"\nnumber of cells is: {CELLS}, length of env: {self.sizesize}m",
                f"\nlength of each grid(dX): {dX}, time step: {dT}, number of steps:{STEPS}",
                f"\nplasma freq for e:{omega_p}, beam energy: {keV} keV equal to {V0} m/s")
  
        # PROGRESS BAR _______________________________________________
        stp = 0
        utils.printProgress(stp, steps, prefix = 'Progress:',
                  suffix = 'Complete', barLength = 50)

        # main loop of PIC ___________________________________________
        for step in range(steps):
            if self.verbose:
                print(f"\n we are in step {step} \n")
            rho = density(self.particles)
            PHI = sor_solver(rho)
            EFIELDn = fieldOnNodes(PHI)
            EFIELD = fieldOnParticles(EFIELDn, self.particles)

            if step == 0:
                self.particles = rewind(-1, EFIELD, self.particles)

            self.particles = moveParticles(EFIELD, self.particles)

            # write to file
            #  ..... TODO ....
    

            # PROGRESS BAR
            stp += 1
            utils.printProgress(stp, steps, prefix = 'Progress:',
                      suffix = 'Complete', barLength = 50)
    
            # plotting _______________________________________________

            if (plot) & (step % plt_stp == 0):
                n_parts = len(self.particles)
                xi = []
                vi = []
                for i in range(n_parts):
                    if self.particles[i].mv:
                        xi.append(self.particles[i].x)
                        vi.append(self.particles[i].v)
        
                fig = plt.figure()      
                plt.scatter(x=xi, y=vi)
                plt.xlabel(" x [position] ")
                plt.ylabel(" v (velocity) ")
                plt.title(f" phase space (x-v) of particles in step {step} ")
                fig.savefig(r"results/fig_"+f"step_{step}")
        
                # ...... TODO ....... 
                # make t more pretty later!

# when this file run: 
if __name__ == "__main__": 

    # instantiate Sim() class for simulation env:
    sim = Sim(size = SIZE, parts_num = NP ,
              dir = RESULT_DIR, verbose = VERBOSE)
    
    # running PIC main loop:
    sim.run(steps=STEPS, omega_p=Omega_p,
            mode= MODE , ampl = AMPL, v0 = V0,
            plot=True, plt_stp = PLT_STP)
     