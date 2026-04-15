"""
Set of functions for validation and extras
"""
import numpy as np

def energy(vel, rho, pot, parameters):
    
    
    # ========================================================================
    
    p, d = parameters["particles"], parameters["domain"]

    m = p["m"]
    dx = d["dx"]

    # ========================================================================
    
    
    E_kin = 0.5 * m * np.sum(vel**2)
    E_pot = 0.5 * np.sum(rho*pot) * dx**2
    E_tot = E_kin + E_pot
    

    return E_tot

def angular_momentum(pos, vel, parameters):
    
    
    # ========================================================================
    
    p, d = parameters["particles"], parameters["domain"]

    m = p["m"]
    L = d["L"]

    # ========================================================================
    
    x = pos[:,0]
    y = pos[:,1]
        
    x_c = np.angle(np.mean(np.exp(1j*2*np.pi*x/L))) * L/(2*np.pi)
    y_c = np.angle(np.mean(np.exp(1j*2*np.pi*y/L))) * L/(2*np.pi)
        
    center = [x_c, y_c]
    
    r = pos - center
    r = r - L * np.round(r/L)
    
    Lz = m * (r[:,0]*vel[:,1] - r[:,1]*vel[:,0])
    
    Lz_tot = np.sum(Lz)
    
    
    return Lz_tot


def Toomre():
    
    
    return 

    
def ostriker_peebles(pos, vel, rho, pot, parameters):
    
    
    # ========================================================================
    
    p, d = parameters["particles"], parameters["domain"]

    m = p["m"]
    L, dx = d["L"], d["dx"]

    # ========================================================================
    
    x = pos[:,0]
    y = pos[:,1]
        
    x_c = np.angle(np.mean(np.exp(1j*2*np.pi*x/L))) * L/(2*np.pi)
    y_c = np.angle(np.mean(np.exp(1j*2*np.pi*y/L))) * L/(2*np.pi)
        
    center = [x_c, y_c]
    
    
    W = 0.5 * np.sum(rho*pot) * dx**2
    
    r = pos - center
    r = r - L * np.round(r/L)
    
    v_theta = (r[:,0] * vel[:,1] - r[:,1] * vel[:,0]) / (np.sqrt(r[:,0]**2 + r[:,1]**2) + dx/2)
    
    T = np.sum(0.5 * m * v_theta**2)
    
    criterion = np.abs(T / W)
   
    
    return criterion

def Jeans_mass():
    
    
    return 

def rotation_curve(pos, vel, parameters, nb=64):
    
    
    # ========================================================================
    
    p, d = parameters["particles"], parameters["domain"]

    m = p["m"]
    L, dx = d["L"], d["dx"]

    # ========================================================================
    
    
    x = pos[:,0]
    y = pos[:,1]
        
    x_c = np.angle(np.mean(np.exp(1j*2*np.pi*x/L))) * L/(2*np.pi)
    y_c = np.angle(np.mean(np.exp(1j*2*np.pi*y/L))) * L/(2*np.pi)
        
    center = [x_c, y_c]
    
    r = pos - center
    r = r - L * np.round(r/L)

    v_theta = (r[:,0] * vel[:,1] - r[:,1] * vel[:,0]) / (np.sqrt(r[:,0]**2 + r[:,1]**2) + dx/2)
    
    b = np.linspace(0, L/2, nb+1)
    bcen = 0.5 * (b[1:] + b[:-1])
    
    r = np.sqrt(r[:,0]**2 + r[:,1]**2)

    Part_per_bin, b = np.histogram(r, bins=b)
    v_rot, b = np.histogram(r, bins=b, weights=v_theta)
    v_rot2, b = np.histogram(r, bins=b, weights=v_theta**2)
    
    div = np.where(Part_per_bin > 0, Part_per_bin, 1)
    v_rot_mean = np.where( Part_per_bin > 0, v_rot / div, 0)
    v_rot2_mean = np.where( Part_per_bin > 0, v_rot2 / div, 0)

    sigma = np.sqrt(np.abs(v_rot2_mean - v_rot_mean**2)) #np.abs not rigorous, here to solve invalid value, to be patched later on


    return bcen, v_rot_mean, sigma