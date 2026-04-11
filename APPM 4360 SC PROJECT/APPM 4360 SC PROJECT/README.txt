APPM 4360 | Schwarz-Christoffel Mapping Project
Authors: Jacob Johnson, Maxwell Meador, William Lampe
================================================

FILES
-----
  main.m              Entry point. Run this script to start the analysis.
  naca_airfoil.m      Generates NACA 4-digit geometry + SC polygon.
  airfoil_analysis.m  Full SC analysis pipeline + five figures.

QUICK START
-----------
  1. Ensure the SC Toolbox (Driscoll & Trefethen) is on your MATLAB path:
       >> addpath(genpath('PATH_TO_SC_TOOLBOX'))

  2. Run:
       >> main

  3. Enter a NACA number (e.g. 2412) and angle of attack (e.g. 5).

  4. Five figures appear:
       Fig A  -- Airfoil geometry
       Fig B  -- Streamlines
       Fig C  -- Velocity magnitude |V|/U_inf
       Fig D  -- Surface pressure coefficient Cp vs x/c
       Fig E  -- Lift curve CL vs alpha

BUGS FIXED
----------
  BUG 1 (naca_airfoil.m) -- Indentation collapse.
    All if/for blocks are now correctly indented. No logic changed,
    but readability and auditability are greatly improved.

  BUG 2 (airfoil_analysis.m) -- Kutta condition prevertex index.
    ORIGINAL: idx_TE = 1 (hardcoded). The SC Toolbox does NOT guarantee
    prevertex ordering, so this produced the wrong trailing-edge angle
    beta_TE and therefore the wrong circulation Gamma for almost every
    airfoil. All downstream quantities (Cp, CL, Fig D, Fig E) were wrong.
    FIX: Evaluate f at every prevertex, compare to known TE location
    w_TE = 1+0i, and use the closest match as idx_TE.

  BUG 3 (airfoil_analysis.m) -- Cp surface sampling.
    ORIGINAL: Sampled theta_surf uniformly on the unit circle in the
    zeta-plane and mapped forward. The resulting w-plane points were
    unordered and not aligned with x/c, so Fig D was scattered or empty.
    FIX: Use the actual airfoil surface vertices from naca_airfoil.m,
    invert the SC map via finverse(), and evaluate velocities there.
    This gives correctly ordered, physically meaningful Cp curves.

  BUG 4 (naca_airfoil.m) -- Polygon closure.
    ORIGINAL: Dropped z_upper(end) AND z_lower(end), which left the
    polygon open (no return to the starting TE vertex).
    FIX: Drop only the shared LE junction duplicate. The SC Toolbox
    polygon() closes the shape implicitly from poly_z(end) to poly_z(1),
    so the TE vertex must appear exactly once.

NOTES
-----
  - Requires MATLAB R2019b or later.
  - No toolboxes beyond SC Toolbox are needed.
  - The SC parameter solve (extermap) takes 10-30 s on a typical laptop.
  - Results are valid for attached flow only (approx. AoA < 15 deg).
  - All figures use non-dimensional coordinates (chord = 1, U_inf = 1).
