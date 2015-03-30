# shallowFoam
OpenFOAM-based solver for 2D shallow water equations

  - Solves the depth-averaged 2D shallow water equations:

    dH/dt + d(HU_i)/d(x_i) = 0
    
    d(HU_i)/dt + d/dx_j (U_j * HU_i) = - g/2 * dH^2/dx_i - g * H * dz_b/dx_i - tau_bx_i / rho + d^2/dx_j^2 ( nu_t * HU )
    
    with:
      - H: flow depth
      - HU: specific discharge
      - U: depth averaged velocity
      - z_b: bottom elevation
      - tau_b: bottom stresses

  - Bottom stresses are modeled via Strickler-equation
  - Turbulence is captured by an eddy viscosity model
  - Works in parallel
  - Captures wet-dry fronts
  - Mesh should have a height of 1 m in z-direction (see tutorials)

