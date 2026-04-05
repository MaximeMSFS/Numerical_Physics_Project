"""
Default parameters used for the simulation if no file is used to resume computation
"""
import numpy as np

eps0 = 1.0
np.random.seed(0)    #for debug, allows for reproductibility

# particles #
Part_num = int(1e4) # Number of particles
q = -1              # Charge of particles (C)
m = 1               # Mass of particles

Pop_A_num = Part_num//2          # Number of particles in population A
Pop_B_num = Part_num - Pop_A_num # Number of particles in population B

Vel_disp = 0.1 # Particles velocity dispertion

# Domain #
L = 1              # Size of the domain (m)
Cell_num = 300     # Number of cells in the domain


# time #
Tmin, Tmax = 0.0, 10.0          # Strating and ending time (s)
dt0 = 1e-3                      # Time step (s)
t = np.arange(Tmin, Tmax, dt0)  # Times array

Relax_step_max = 50
epsilon = 1e-4

eps = dx/2  # dumb schmilblick to be removed later on
save_step = 100
