import os
from constants import *
import matplotlib.pyplot as plt

from generation import two_stream as parts_gen
from cycle import density, sor_solver, fieldOnNodes, \
fieldOnParticles, rewind, moveParticles
import utils


class Sim():
    """class for simulation of pic
    """
    def __init__(self, size: float, cells_num, time_step, parts_num, dir, verbose = False):
        self.size = size
        self.cells_num = cells_num
        self.parts_num = parts_num  # NPpc * cells_num
        self.time_step = time_step
        self.dx = size / cells_num
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
    

    def run(self,
            parts_gen, # particle generator such two stream or electron beam
            steps,
            omega_p,
            mode , ampl ,
            v0 ,
            plot, plt_stp ):
        
        # generating particle:
        self.particles = parts_gen(omega_p=omega_p, parts_num = self.parts_num, size = self.size, # for Particle class
                  mode= mode , ampl = ampl, v0 = v0 )
        if self.verbose:
            # ..... TODO ... : print other parameters too
            print( " particles in t=0 generated with ...")
            print(f"\nnumber of cells is: {self.cells_num}, length of env: {self.sizesize}m",
                f"\nlength of each grid(dX): {self.dx}, time step: {self.dt}, number of steps:{STEPS}",
                f"\nplasma freq for e:{omega_p}, beam velocity: {v0} m/s")
  
        # PROGRESS BAR _______________________________________________
        stp = 0
        utils.printProgress(stp, steps, prefix = 'Progress:',
                  suffix = 'Complete', barLength = 50)

        # main loop of PIC ___________________________________________
        for step in range(steps):
            if self.verbose:
                print(f"\n we are in step {step} \n")
            
            rho = density(self.particles, self.cells_num, self.dx)
            if self.verbose: 
                print("\n wheiting dencity on grids from particle position (rho) by dencity funcion done! \n")
    
            PHI = sor_solver(rho, self.cells_num, self.dx)
            if self.verbose:
                print(f"Poisson solved by SOR method")

                
            EFIELDn = fieldOnNodes(PHI, self.cells_num, self.dx)
            if self.verbose:
                print("\n Calculating field from potential on grids (efild) done!\n")

                
            EFIELD = fieldOnParticles(EFIELDn, self.particles, self.cells_num, self.dx)
            if self.verbose:       
                print("\n calculating E field on particles done! \n")

            if step == 0:
                self.particles = rewind(-1, EFIELD, self.particles, dt=self.time_step)

            self.particles = moveParticles(EFIELD, self.particles, self.size, dt= self.time_step)
            if self.verbose:
                print("\n moving of particles done! \n")

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


     