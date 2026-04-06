""" 
This code computes the 1D dynamics of charged particles in a plasma with no 
external field using relaxation of the Poisson equation or Fourrier techniques
with various integration schemes.
"""

import numpy as np
import matplotlib.pyplot as plt
from colorama import init, Fore, Style
from tqdm import tqdm            #progress bar nothing usefull for the code, just nicer, remove if tqdm not installed or install with 'conda install tqdm' in the spyder terminal (bottom right)
init(autoreset=True)

import modules.plotting as plotting
from modules.user_input import User_input
import modules.init as initialize
import modules.schemes as schemes


################################# Parameters ##################################


parameters, main_file, pos_file, vel_file, pot_file, resume_file, save = User_input()

p, d, t, n = parameters["particles"], parameters["domain"], parameters["time"], parameters["numerical"]

m = p["m"]

L, Cell_num, dx = d["L"], d["Cell_num"],  d["L"] / d["Cell_num"]
Cell_pos = d["Cell_pos"]
eps0 = d["eps0"]

Tmin, Tmax, dt0 = t["Tmin"], t["Tmax"], t["dt0"]
t_arr = t["t"]
initial_step = t["step"]

save_step = n["save_step"]
scheme = n["scheme"]


################################# Initialisation ##############################


Part_pos, Part_vel, Part_charge, Cell_pot, rho0, Part_per_cell = initialize.start(parameters)

plotting.conditions_plot(Cell_pos, Part_per_cell, Part_pos, Part_vel, 'Initial', Tmin, parameters)


################################ Main Program #################################
    

dt = dt0
E_tot = np.zeros_like(t_arr)
        
#for i in range(len(t)):                            # Code to use without tqdm
for i in tqdm(range(initial_step, len(t_arr)), desc="sim progress"):

    Part_force, Part_per_cell, Cell_pot, Cell_field = schemes.Compute_force(Part_pos, Cell_pot, parameters)

    # Integration Scheme 
    if (scheme == 0) :
        Part_pos, Part_vel = schemes.Euler_exp(Part_pos, Part_vel, dt, Part_force, parameters)
    elif  (scheme == 1) :
        Part_pos, Part_vel = schemes.Euler_imp(Part_pos, Part_vel, dt, Part_force, parameters)
    elif  (scheme == 2) :
        Part_pos, Part_vel, Cell_pot, Cell_field = schemes.Leapfrog(Part_pos, Part_vel, dt, Part_force, Cell_pot, parameters)
    elif  (scheme == 3) :
        Part_pos, Part_vel, Cell_pot, Cell_field = schemes.RK4(Part_pos, Part_vel, dt, Part_force, Cell_pot, parameters)

    Part_pos = Part_pos % L                                            # Periodic BC : % L -> 'modulo L'
    
    # Post processing, validation and other
    
    """dt = min(min(dx/Part_vel), dt0)
    print(f'dt = {dt}')"""
    
    E_kin = 0.5 * m * np.sum(Part_vel**2)
    E_pot = 0.5 * eps0 * np.sum(Cell_field**2) * dx
    E_tot[i] = E_kin + E_pot
    
    if (save == 'yes') and (i % save_step == 0):
        main_file.write(f"{i:6d} {t_arr[i]:10.3f} {E_tot[i]:10.3e}\n")
        np.save(pos_file, Part_pos)
        np.save(vel_file, Part_vel)
        np.save(pot_file, Cell_pot)
        
        np.savez(resume_file, step=i, time=t_arr[i], pos=Part_pos, vel=Part_vel, pot=Cell_pot)

        pos_file.flush()
        vel_file.flush()
        pot_file.flush()
    
    #print (f'{t_arr[i]/Tmax*100} %' )                                     # code progress check without tqdm
    
print (Fore.GREEN + Style.BRIGHT + '\nSimulation completed successfully !\n ')
    
if (save == 'yes') :
    main_file.close()
    pos_file.close()
    vel_file.close()
    pot_file.close()
    
    
################################ Plots ########################################


plotting.conditions_plot(Cell_pos, Part_per_cell, Part_pos, Part_vel, 'Final', Tmax, parameters)

plotting.field_plot(Cell_pos, Cell_pot, Cell_field, 'Final', Tmax, parameters)

plotting.energy_plot(E_tot, parameters)

plt.show()