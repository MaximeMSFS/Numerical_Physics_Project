""" 
This code computes the 1D dynamics of charged particles in a plasma with no 
external field using relaxation of the Poisson equation or Fourrier techniques
with various integration schemes.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as cst
from pathlib import Path
from colorama import init, Fore, Style
from tqdm import tqdm            #progress bar nothing usefull for the code, just nicer, remove if tqdm not installed or install with 'conda install tqdm' in the spyder terminal (bottom right)
init(autoreset=True)

import modules.file_handling as fh
import modules.plotting as plotting

################################# Parameters ##################################


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
dx = L / Cell_num  # Size of a cell (m)

Cell_pos = np.linspace(0, L-dx, Cell_num) + dx/2 # position of each cell (center of the cell)
Cell_edges = np.linspace(0, L, Cell_num+1)       # position of each cell borders (used for bins)

# FFT wave numbers #
k = 2*np.pi*np.fft.fftfreq(Cell_num, d=dx)
k2 = k**2
k2[0] = 1.0   # avoids division by zero

# time #
Tmin, Tmax = 0.0, 10.0          # Strating and ending time (s)
dt0 = 1e-2                      # Time step (s)
t = np.arange(Tmin, Tmax, dt0)  # Times array

Relax_step_max = 50
epsilon = 1e-4

eps = dx/2  # dumb schmilblick to be removed later on
save_step = 100


################################# Initialisation ##############################


# here in case we keep track of the history (need for t), otherwise can be moved upward
Part_pos = np.zeros(Part_num)
Part_vel = np.zeros(Part_num)

Part_pos[0:Pop_A_num] = L*np.random.rand(Pop_A_num)
Part_pos[Pop_A_num:] = L*np.random.rand(Pop_B_num)

Part_vel[0:Pop_A_num] = 1.0 + Vel_disp * np.random.rand(Pop_A_num)
Part_vel[Pop_A_num:] = -1.0 + Vel_disp * np.random.rand(Pop_B_num)

Part_charge = np.ones(Part_num) * q

V = np.zeros(Cell_num)

#one particle type
rho0 = Part_num / Cell_num # Intrinsic charge density

#test of initial conditions
Part_per_cell, b = np.histogram(Part_pos, bins=Cell_edges)


################################### Functions #################################


def Compute_force(pos):

    global V, Cell_field

    Part_per_cell, b = np.histogram(pos, bins=Cell_edges)

    #rho = rho0 - Part_per_cell                           # charge density
    rho = q*(Part_per_cell - rho0)/dx
    
    
    # Relaxation scheme #
    if solver == 'p':
        
        for j in range(Relax_step_max):                                        
            Vnew = 0.5 * (np.roll(V,-1) + np.roll(V,1) + dx**2 * (rho / eps0)) # Poisson equation 
        
            delta=np.max(np.abs((Vnew-V)/(V+1e-10)))                           # Convergency check
        
            V = Vnew.copy()
     
            if delta<epsilon :                                                 # Optimisation
                break

        Cell_field = -(np.roll(V,-1) - np.roll(V,1)) / (2*dx)                  # Field in each cell

    
    # Fourier solver #
    else:

        rho_k = np.fft.fft(rho)

        V_k = rho_k / (k2 * eps0)
        V_k[0] = 0.0

        V = np.real(np.fft.ifft(V_k))

        E_k = -1j*k*V_k
    
        Cell_field = np.real(np.fft.ifft(E_k))  # Field in each cell


    Part_field = np.interp(pos, Cell_pos, Cell_field)             # Field felt by each particle (simple interpolation)

    Part_force = q * Part_field                                   # Force on each particle

    return Part_force, Part_per_cell


def Euler_exp(pos, vel, dt, force):

    dpos_dt = vel
    dvel_dt = force / m

    pos_new = pos + dpos_dt * dt
    vel_new = dpos_dt + dvel_dt * dt

    return pos_new, vel_new


def Euler_imp(pos, vel, dt, force):
    
    dpos_dt = vel
    dvel_dt = force / m
    
    vel_new = dpos_dt + dvel_dt * dt
    pos_new = pos + vel_new * dt

    return pos_new, vel_new


def Leapfrog(pos, vel, dt, force):
    
    dpos_dt = vel
    dvel_dt = force / m
    
    pos_new = pos + dpos_dt * dt + 0.5 * dvel_dt * dt**2

    force_new, Part_per_cell = Compute_force(pos_new)

    dvel_dt_new = force_new / m

    vel_new = vel + 0.5 * (dvel_dt + dvel_dt_new) * dt

    return pos_new, vel_new


def RK4(pos, vel, dt, force):

    k1_pos = vel
    k1_vel = force / m
    
    pos_temp = (pos + 0.5 * dt * k1_pos) % L
    vel_temp = vel + 0.5 * dt * k1_vel
    force, Part_per_cell = Compute_force(pos_temp)
    k2_pos = vel_temp
    k2_vel = force / m

    pos_temp = (pos + 0.5 * dt * k2_pos) % L
    vel_temp = vel + 0.5 * dt * k2_vel
    force, Part_per_cell  = Compute_force(pos_temp)
    k3_pos = vel_temp
    k3_vel = force / m

    pos_temp = (pos + dt * k3_pos) % L
    vel_temp = vel + dt * k3_vel
    force, Part_per_cell  = Compute_force(pos_temp)
    k4_pos = vel_temp
    k4_vel = force / m

    pos_new = pos + dt * (k1_pos + 2*k2_pos + 2*k3_pos + k4_pos)/6
    vel_new = vel + dt * (k1_vel + 2*k2_vel + 2*k3_vel + k4_vel)/6

    pos_new = pos_new % L

    return pos_new, vel_new


################################ User Inputs ##################################

Invalid_entry = True
while Invalid_entry :
    init = save = input(Fore.CYAN + "Do you wish to resume calculation from a previously computed file ?\n (Yes, No) :\n").strip().lower()
    if (init == 'yes') or (init == 'no') :
        Invalid_entry = False
        if (init == "yes") :
            input_dir = input(Fore.YELLOW + "/!\ Not Working yet /!\ \nPlease enter the directory of the saved files :\n").strip()
    else :
        print (Fore.RED + 'Invalid entry')

Invalid_entry = True
while Invalid_entry :
    scheme = int(input(Fore.CYAN + "Which scheme do you wish to use ?\n (0 = Explicite Euler, 1 = Implicit Euler, 2 = Leapfrog, 3 = RK4) :\n").strip())
    if (scheme == 0) or (scheme == 1) or (scheme == 2) or (scheme == 3) :
        Invalid_entry = False
    else :
        print (Fore.RED + 'Invalid entry')

Invalid_entry = True
while Invalid_entry :
    solver = input(Fore.CYAN + "Which solver do you wish to use ?\n (F = Fourrier, P = Poisson) :\n").strip().lower()
    if (solver == 'f') or (solver == "p") :
        Invalid_entry = False
    else :
        print (Fore.RED + 'Invalid entry')
        
parameters = dict(physics=dict(q=q, m=m, eps0=eps0),
                  particles=dict(Part_num=Part_num, Pop_A_num=Pop_A_num, Pop_B_num=Pop_B_num, Vel_disp=Vel_disp), 
                  domain=dict(L=L, Cell_num=Cell_num), 
                  time=dict(Tmin=Tmin, Tmax=Tmax, dt=dt0, t=t), 
                  numerical=dict(scheme=scheme, solver=solver))

Invalid_entry = True
while Invalid_entry :
    save = input(Fore.CYAN + "Do you wish to save the results ?\n (Yes, No) :\n").strip().lower()
    if (save == 'yes') or (save == "no") :
        Invalid_entry = False
        if (save == "yes") :
            output_dir = input(Fore.CYAN + "Please enter the output directory name :\n").strip()
            data_dir = Path("data")/output_dir
            data_dir.mkdir(exist_ok=True)
            main_file = open(data_dir/f"{output_dir}_main.dat", 'w')
            pos_file   = open(data_dir /f"{output_dir}_position.dat", "ab")
            vel_file   = open(data_dir /f"{output_dir}_velocity.dat", "ab")
            pot_file = open(data_dir /f"{output_dir}_potential.dat", "ab")
            fh.header(main_file, parameters)
    else :
        print (Fore.RED + 'Invalid entry')

if (scheme == 0) :
    print (Fore.WHITE + Style.BRIGHT + '\nUsing Explicite Euler Scheme ', end="")
elif  (scheme == 1) :
    print (Fore.WHITE + Style.BRIGHT + '\nUsing Implicit Euler Scheme ', end="")
elif  (scheme == 2) :
    print (Fore.WHITE + Style.BRIGHT + '\nUsing Explicite Leapfrog Scheme ', end="")
else :
    print (Fore.WHITE + Style.BRIGHT + '\nUsing Explicite Range-Kutta 4 Scheme ', end="")

    
if (solver == 'f') :
    print (Fore.WHITE + Style.BRIGHT + 'with the Fourrier Solver')
else :
    print (Fore.WHITE + Style.BRIGHT + 'with the Poisson Solver')

if (save == "yes") :
    print (Fore.WHITE + Style.BRIGHT + f'Saving results to {data_dir}')
print ('\n')

    
################################ Main Program #################################

plotting.conditions_plot(Cell_pos, Part_per_cell, Part_pos, Part_vel, 'Initial', Tmin, parameters)
    
dt = dt0
E_tot = np.zeros_like(t)
        
#for i in range(len(t)):                            # Code to use without tqdm
for i in tqdm(range(len(t)), desc="sim progress"):

    Part_force, Part_per_cell = Compute_force(Part_pos)

    # Integration Scheme 
    if (scheme == 0) :
        Part_pos, Part_vel = Euler_exp(Part_pos, Part_vel, dt, Part_force)
    elif  (scheme == 1) :
        Part_pos, Part_vel = Euler_imp(Part_pos, Part_vel, dt, Part_force)
    elif  (scheme == 2) :
        Part_pos, Part_vel = Leapfrog(Part_pos, Part_vel, dt, Part_force)
    elif  (scheme == 3) :
        Part_pos, Part_vel = RK4(Part_pos, Part_vel, dt, Part_force)

    Part_pos = Part_pos % L                                            # Periodic BC : % L -> 'modulo L'
    
    # Post processing, validation and other
    
    """dt = min(min(dx/Part_vel), dt0)
    print(f'dt = {dt}')"""
    
    E_kin = 0.5 * m * np.sum(Part_vel**2)
    E_pot = 0.5 * eps0 * np.sum(Cell_field**2) * dx
    E_tot[i] = E_kin + E_pot
    
    if (save == 'yes') and (i % save_step == 0):
        main_file.write(f"{i:6d} {t[i]:10.3f} {E_tot[i]:10.3e}\n")
        np.save(pos_file, Part_pos)
        np.save(vel_file, Part_vel)
        np.save(pot_file, Cell_field)

        pos_file.flush()
        vel_file.flush()
        pot_file.flush()
    
    #print (f'{t[i]/Tmax*100} %' )                                     # code progress check without tqdm
    
print (Fore.GREEN + Style.BRIGHT + '\nSimulation completed successfully !\n ')
    
if (save == 'yes') :
    main_file.close()
    pos_file.close()
    vel_file.close()
    pot_file.close()
    
    
################################ Plots #########################################


plotting.conditions_plot(Cell_pos, Part_per_cell, Part_pos, Part_vel, 'Final', Tmax, parameters)

plotting.field_plot(Cell_pos, V, Cell_field, 'Final', Tmax, parameters)

plotting.energy_plot(E_tot, parameters)

plt.show()