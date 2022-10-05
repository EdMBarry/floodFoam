# floodFoam
OpenFOAM-based solver for 2D shallow water equations. 

It is important to note that this solver has been renamed to floodFoam (from shallowFoam).
This branch has been updated to work on the latest version of OpenFOAM from the OpenFOAM Foundation (OpenFOAM-10).
This code should work for the equivalent ESI group OpenFOAM version but requires testing.

Authors:
  - KM-Turbulenz GmbH (www.km-turbulenz.de), 2009
  - Florian Mintgen, 2012
  - Dr. Ed Barry, FloodMapp, 2022

Description:
  - Solves the depth-averaged 2D shallow water equations:

    dH/dt + d(HU_i)/d(x_i) = 0
    
    d(HU_i)/dt + d/dx_j (U_j * HU_i) = - g/2 * dH^2/dx_i - g * H * dz_b/dx_i - tau_bx_i / rho + d^2/dx_j^2 ( nu_t * HU )
    
    with:
      - H: flow depth
      - HU: specific discharge
      - U: depth averaged velocity
      - z_b: bottom elevation
      - tau_b: bottom stresses

  - Bottom stresses are modeled via Strickler-equation (and Manning's equation to be included).
  - Turbulence is captured by an eddy viscosity model
  - Works in parallel
  - Captures wet-dry fronts
  - Mesh should have a height of 1 m in z-direction (see tutorials)

  - Main advantages over shallowWaterFoam (the shallow water solver in the official OpenFOAM repository):
    - Explicit formulation of flow depth and bottom elevation
    - Bottom stresses / surface roughness taken into account
    - Custom  boundary conditions well suited for river hydraulics
    - Avoids divide-by-zero errors for dry areas which shallowWaterFoam cannot handle
