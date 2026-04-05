````markdown
# 1D Plasma Simulation

This project simulates the **1D dynamics of charged particles in a plasma** using a particle-in-cell (PIC) program.  
The electric potential is computed either through **relaxation of the Poisson equation** or using **Fourier techniques**. Particles motion is integrated with several numerical schemes.

---

## Features

- 1D particle-in-cell plasma simulation
- Two particle populations with configurable velocity dispersion
- Periodic boundary conditions
- Multiple time integration schemes:
  - Explicit Euler
  - Implicit Euler
  - Leapfrog
  - Runge–Kutta 4 (RK4)
- Two solvers for the electrostatic potential:
  - Poisson solver (relaxation techniques)
  - Fourier solver (FFT)
- OPtion to save the simulation data during runtime
- Visualization of:
  - Particle distribution
  - Phase-space distribution
  - Electric potential and field
  - Total energy evolution

---

## Usage

Run the simulation script:

```bash
python Project_Main_3.py
```

During execution the program will prompt for several options:

1. Whether to load initial conditions from a previous simulation
2. Which numerical integration scheme is to be uses
3. Which electrostatic solver is to be used
4. Whether to save the results

All simulation data will be stored inside the `data/` directory, to a folder given by the user.

Example structure:

```
data/
└── user_output/
    ├── user_output_main.dat
    ├── user_output_position.dat
    ├── user_output_velocity.dat
    └── user_output_potential.dat
```

---

## Output Files

| File              | Description                                                 |
| ----------------- | ------------------------------------------------------------|
| `*_main.dat`      | Main simulation log containing the set of used parameters,  |
|                   | step number, time, and total energy as ASCII format         |
| `*_position.dat`  | Particle positions as Binary format                         |
| `*_velocity.dat`  | Particle velocities as Binary format                        |
| `*_potential.dat` | Electric potential in each cell as Binary format            |

Binary files are written using **NumPy's binary format**, which allows efficient storage and easy restart of simulations.

---

## Parameters

Simulation parameters are defined in the parameter module and include:

### Physics

* Particle charge `q`
* Particle mass `m`
* Medium permittivity `eps0`

### Particles

* Total number of particles
* Population sizes
* Velocity dispersion

### Domain

* Domain length `L`
* Number of spatial cells `Cell_num`

### Time

* Starting time `Tmin`
* Ending time `Tmax`
* Fixed time step `dt0`

### Numerical Settings

* Integration scheme
* Solver type
* Relaxation parameters
* Save interval

---

