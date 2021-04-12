# ShockTubeModel
CFD model for compressible 1D Shock Tube (SOD) implementing various schemes.
https://en.wikipedia.org/wiki/Sod_shock_tube

## Cases
### Case 1
Variable | Left | Right
---------|------|------
Pressure |10^{5}|10^{4}

### Case 2

## User Inputs
### Initial values for left and right side of the diaphragm.
- Pressure, p_L and p_R
- Density, rho_L and rho_R
- Velocity, u_L and u_R
### Grid Settings
- Number of grid points, nx
- Initial diaphragm position, x0
- Front of plate, xmin
- End of plate, xmanx
### Solver Settings
- Stability condition, cfl
- Solver start time, starttime
- Solver end time, runtime
### Additional settings
- [ ] User designated postproccessing variables, e.g. Mach number, mass flow rate, entropy, etc.
- [ ] Read from input file
- [ ] Output data type

## Implemented Schemes
Explicit schemes:
- [x] Central Difference
- [x] Lax-Wendroff
- [x] Van Leer Flux-Vector Splitting
- [x] MacCormack

## Outputs
### Post-Processing
Outputs pressure, entropy (approx.), velocity, Mach number, density, and mass flow rate.
### Data
Data written to .csv for now.
