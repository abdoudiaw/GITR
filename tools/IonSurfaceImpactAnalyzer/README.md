# IonSurfaceImpactAnalyzer

Analyze the distribution of ion impacts on a flat surface.

The primary goal of PMI (Plasma Material Interaction) is to study ionic transport in fusion devices through Monte Carlo transport.

### Features: 

##### Physics : 
 
- ###### Numerical Integrator:
  ```ruby
  Boris.
  ```
- ###### Collisions:
  ```ruby
  None.
  ```
- ###### Normalization:
```ruby 
     Distances are normalized to the ion Larmor radius (rho_L).
     Time is normalized with the cyclotron frequency (w_c).
```

- ###### Input:
```ruby 
- input.in:  file containing input parameters. This file is organized as followed:
    npart           # Number of particles
    dens_back       # Background density
    zbar_back       # Mean ion charge (background)
    pmass_back      # Mass of background plasma
    tempi_back      # Ion temperature (background)
    tempe_back      # Electron temperature (background)
    pmass_imp       # Mass of impurity
    zbar_imp        # Mean ion charge (impurity)
    tempi_imp       # Ion temperature (impurity)
    dt              # Time step size
    steps           # Number of time steps
    snaps           # Snapshots or outputs
    sheath_factor   # Factor for sheath calculation
    Bmagnitude      # Magnetic field strength
    alpha           # Angle of incidence
    xmin, xmax      # x boundaries
    ymin, ymax      # y boundaries
    zmin, zmax      # z boundaries
```
