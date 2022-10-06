# floodFoam
OpenFOAM-based solver for 2D shallow water equations. 

It is important to note that this solver has been renamed to floodFoam (from shallowFoam).
This branch has been updated to work on the latest version of OpenFOAM from the OpenFOAM Foundation (OpenFOAM-10).
This code should work for the equivalent ESI group OpenFOAM version but requires testing.

## Authors:
  - KM-Turbulenz GmbH (www.km-turbulenz.de), 2009
  - Florian Mintgen, 2012
  - Dr. Ed Barry, FloodMapp, 2022

## Description:
  - Solves the depth-averaged 2D shallow water equations:

  $$\frac{\partial H}{\partial t} + \frac{\partial (H\mathbf{U}_i)}{\partial x_i} = 0   $$


  $$ 
        \frac{\partial (HU_i)}{\partial t}  + \frac{\partial (U_j HU_i)}{\partial x_j} = -\frac{g}{2} \cdot \frac{\partial H^2}{\partial x_i}  - gH \frac{\partial z_b}{\partial x_i}  - \frac{\tau_b}{\rho} + \frac{\partial ^2(\nu_t HU)}{\partial j^2} $$
    
*with*:
  - $H$: flow depth
  - $H\mathbf{U}$: specific discharge
  - $\mathbf{U}$: depth averaged velocity
  - $z_b$: bottom elevation
  - $\tau_b$: bottom stresses


**Notes**:
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
