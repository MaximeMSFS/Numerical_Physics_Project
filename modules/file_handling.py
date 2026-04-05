"""
file handling functions for saving and starting back from an old computation
"""

def header(file, parameters):
    q = parameters["physics"]["q"]
    m = parameters["physics"]["m"]
    eps0 = parameters["physics"]["eps0"]
    
    Part_num = parameters["particles"]["Part_num"]
    Pop_A_num = parameters["particles"]["Pop_A_num"]
    Pop_B_num = parameters["particles"]["Pop_B_num"]
    Vel_disp = parameters["particles"]["Vel_disp"]
    
    L = parameters["domain"]["L"]
    Cell_num = parameters["domain"]["Cell_num"]
    
    Tmin = parameters["time"]["Tmin"]
    Tmax = parameters["time"]["Tmax"]
    dt0 = parameters["time"]["dt"]
    
    scheme = parameters["numerical"]["scheme"]
    solver = parameters["numerical"]["solver"]
    
    
    file.write("==============================================================\n")
    file.write("                  1D PLASMA SIMULATION OUTPUT                 \n")
    file.write("==============================================================\n")
    
    file.write("computed with Project_Main_3.py\n\n")
    
    file.write("--- Simulation Parameters ---\n")
    file.write(f"scheme              : {scheme}\n")
    file.write(f"solver              : {solver}\n\n")
    
    file.write("---- Physical Parameters ----\n")
    file.write(f"particles number    : {Part_num}\n")
    file.write(f"Epsilon0            : {eps0}\n")
    file.write(f"charge              : {q}\n")
    file.write(f"mass                : {m}\n\n")
    
    file.write("----- Domain Parameters -----\n")
    file.write(f"domain length       : {L}\n")
    file.write(f"cell number         : {Cell_num}\n\n")
    
    file.write("------ Time Parameters ------\n")
    file.write(f"Tmin                : {Tmin}\n")
    file.write(f"Tmax                : {Tmax}\n")
    file.write(f"dt                  : {dt0}\n\n")
 
    file.write("--- Particles Populations ---\n")
    file.write(f"population A number : {Pop_A_num}\n")
    file.write(f"population B number : {Pop_B_num}\n")
    file.write(f"velocity dispersion : {Vel_disp}\n\n")

    file.write("-------- Random Seed --------\n")
    file.write("numpy random seed    : 0\n\n")
    
    file.write("==============================================================\n\n")
    
    file.write(f"{'step':>6} {'time':>10} {'energy':>10}\n")
    return
