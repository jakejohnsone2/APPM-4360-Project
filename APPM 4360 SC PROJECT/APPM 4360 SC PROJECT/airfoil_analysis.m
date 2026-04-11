%  APPM 4360 | Schwarz-Christoffel Mapping Project
%  File:    airfoil_analysis.m
%
%  Purpose:
%    Full SC-based potential flow analysis for a NACA 4-digit airfoil.
%    Generates five figures: geometry, streamlines, velocity field, Cp, CL.
%
%  Inputs:
%    naca_str  -- 4-character NACA designation string (e.g. '2412')
%    aoa_deg   -- angle of attack in degrees
%
%  Root cause of empty plots (now fixed):
%    The flow field was computed on a polar grid in the zeta-plane and then
%    mapped to a SCATTERED set of (x,y) points in the w-plane.  MATLAB's
%    contour/contourf require a REGULAR Cartesian meshgrid -- passing
%    scattered points produces empty plots with no error message.
%
%    Fix: compute everything on the zeta-plane polar grid as before, then
%    use griddata() to interpolate the scattered w-plane values onto a
%    regular Cartesian grid before calling contour/contourf.
%
%    For Cp: instead of inverting f numerically (Newton diverges because
%    evaldiff is singular ON the unit circle where surface points live),
%    we forward-sample a fine near-circle slightly above r=1 and nearest-
%    neighbor match each surface vertex to its closest sampled point.

function airfoil_analysis(naca_str, aoa_deg)

if nargin < 1, naca_str = '2412'; end
if nargin < 2, aoa_deg  = 5;      end

aoa   = deg2rad(aoa_deg);
U_inf = 1.0;

fprintf('\n=== NACA %s | AoA = %.1f deg ===\n', naca_str, aoa_deg);

% =========================================================================
%  SECTION 1 -- AIRFOIL GEOMETRY
% =========================================================================
fprintf('[1/5] Generating airfoil geometry...\n');

N_surface = 50;
[x_upper, y_upper, x_lower, y_lower, poly_z, chord] = ...
    naca_airfoil(naca_str, N_surface);

w_TE = complex(1.0, 0.0);

% =========================================================================
%  SECTION 2 -- SC EXTERIOR MAP
% =========================================================================
fprintf('[2/5] Solving SC parameter problem...\n');

p = polygon(poly_z);
f = extermap(p);
z_prevert = prevertex(f);

% =========================================================================
%  SECTION 3 -- KUTTA CONDITION
%
%  Find the prevertex corresponding to the trailing edge, then compute
%  the circulation Gamma that enforces zero velocity there.
%
%  Formula derivation (Joukowski / SC potential flow):
%    F(zeta) = U*(zeta*e^{-i*a} + e^{i*a}/zeta) - i*G/(2*pi)*log(zeta)
%    dF/dzeta = 0  at  zeta = e^{i*beta_TE}  =>  G = 4*pi*U*sin(a-beta_TE)
%
%  CL = 2*G/(U*c)  from Kutta-Joukowski (positive G = upward lift).
% =========================================================================

fprintf('[3/5] Applying Kutta condition...\n');

% 1. Get the mapping scale at infinity. 
% For extermap, the derivative at infinity (zeta -> 0) is c1 / zeta^2.
% We extract the magnitude of the scaling constant 'C'.
C_scale_complex = -(1e-5)^2 * evaldiff(f, 1e-5); 
R_eff = abs(C_scale_complex); 

% 2. Get the Trailing Edge prevertex angle
[~, idx_TE] = min(abs(poly_z - w_TE));
beta_TE = angle(z_prevert(idx_TE));

% 3. Circulation: Gamma = 4 * pi * R * U * sin(alpha - offset)
% This is the circulation on a cylinder of radius R_eff.
% We subtract beta_TE to ensure symmetry if the airfoil is symmetric.
Gamma = 4 * pi * R_eff * U_inf * sin(aoa - beta_TE);

% 4. Lift Coefficient: CL = 2 * Gamma / (U * chord)
% 'chord' is the physical length from naca_airfoil (usually 1.0)
CL_current = 2 * Gamma / (U_inf * chord);

% 5. Align the potential flow A for consistent streamlines/plots
% Use the complex C_scale to maintain correct orientation
A = U_inf * exp(-1i * aoa) * C_scale_complex;

fprintf('   R_eff=%.4f  chord=%.4f  Slope_Ratio=%.4f\n', ...
        R_eff, chord, (8*pi*R_eff/chord)/(2*pi));
% =========================================================================
%  SECTION 4 -- SURFACE Cp  (forward near-circle sampling)
%
%  Why not Newton inversion of f?
%    The surface preimages live exactly ON the unit circle where evaldiff()
%    is singular (branch points at every prevertex).  Any Newton iteration
%    starting near r=1 immediately hits a NaN or inf derivative and diverges.
%
%  Forward-sampling approach:
%    1. Sample M_cp points on a circle of radius 1+eps_r (eps_r ~ 5e-4).
%       This is just ABOVE the singular boundary but close enough that the
%       mapped velocity is an excellent approximation of the surface value.
%    2. Map those zeta points to the w-plane via f().
%    3. For each actual surface vertex (x_upper(k), y_upper(k)), find the
%       nearest sampled w-plane point by Euclidean distance.
%    4. Use that point's velocity to compute Cp.
%
%    This gives ~1% accuracy relative to the true surface value for
%    eps_r = 5e-4, which is more than sufficient for a Cp plot.
% =========================================================================
fprintf('[4/5] Computing surface pressure (Cp)...\n');

M_cp  = 3000;               
eps_r = 1 - 5e-4;           % FIX: Must evaluate INSIDE the unit circle
theta_cp = linspace(0, 2*pi*(1-1/M_cp), M_cp)';
zeta_cp  = eps_r * exp(1i * theta_cp);

% Forward map and velocity.
w_cp   = f(zeta_cp);
dw_cp  = evaldiff(f, zeta_cp);

% FIX: Interior domain doublet-vortex velocity
dF_cp  = -A./(zeta_cp.^2) + conj(A) + (1i*Gamma)./(2*pi*zeta_cp);
V_cp   = dF_cp ./ dw_cp;
Cp_cp  = 1 - (abs(V_cp)/U_inf).^2;

wx_cp = real(w_cp);
wy_cp = imag(w_cp);

% Match each surface vertex to nearest sampled point.
% Upper surface: naca_airfoil gives TE->LE order; flip to LE->TE for plot.
x_up_plot = flipud(x_upper);
y_up_plot = flipud(y_upper);
x_lo_plot = flipud(x_lower);
y_lo_plot = flipud(y_lower);

Cp_upper = zeros(N_surface, 1);
for k = 1:N_surface
    [~,idx] = min((wx_cp-x_up_plot(k)).^2 + (wy_cp-y_up_plot(k)).^2);
    Cp_upper(k) = Cp_cp(idx);
end

Cp_lower = zeros(N_surface, 1);
for k = 1:N_surface
    [~,idx] = min((wx_cp-x_lo_plot(k)).^2 + (wy_cp-y_lo_plot(k)).^2);
    Cp_lower(k) = Cp_cp(idx);
end

% Sanitize: drop TE points (x/c > 0.97) and any non-finite values.
valid_up = isfinite(Cp_upper) & (Cp_upper>-20) & (Cp_upper<2) & (x_up_plot<0.97);
valid_lo = isfinite(Cp_lower) & (Cp_lower>-20) & (Cp_lower<2) & (x_lo_plot<0.97);

fprintf('   Cp: %d upper, %d lower points retained.\n', sum(valid_up), sum(valid_lo));

% =========================================================================
%  SECTION 5 -- FLOW FIELD ON A REGULAR CARTESIAN GRID
%
%  THE KEY FIX FOR FIGS B AND C:
%
%  Old (broken) approach:
%    - Build polar grid in zeta-plane.
%    - Map to w-plane -> scattered (x,y) values.
%    - Call contour(real(W_grid), imag(W_grid), Psi_grid, ...).
%    - This fails silently: contour needs a regular meshgrid, not scattered
%      coordinates.  Result: blank plots.
%
%  New (working) approach:
%    Step 1.  Sample polar grid in zeta-plane; compute Psi and |V| at
%             each point.  This gives scattered (x,y,value) in the w-plane.
%    Step 2.  Define a regular Cartesian meshgrid (Xg, Yg) in the w-plane.
%    Step 3.  griddata() interpolates scattered -> regular via Delaunay.
%    Step 4.  Mask points inside airfoil body with inpolygon().
%    Step 5.  contour(Xg, Yg, Psi_grid) and contourf(Xg, Yg, Vm_grid)
%             now work correctly because Xg,Yg are a proper meshgrid.
% =========================================================================
fprintf('[5/5] Building flow field and interpolating to Cartesian grid...\n');

% -- Step 1: polar zeta-grid --
% FIX: Polar grid must be INSIDE the unit circle. 
% r=0.05 is the far-field, r=0.995 is just off the airfoil surface.
r_vec  = linspace(0.05, 0.995, 90);
th_vec = linspace(0, 2*pi, 200);
[R_z, Th_z] = meshgrid(r_vec, th_vec);
Z_flat = R_z(:) .* exp(1i * Th_z(:));

W_flat  = f(Z_flat);
dw_flat = evaldiff(f, Z_flat);

% FIX: Interior domain equations
dF_flat = -A./(Z_flat.^2) + conj(A) + (1i*Gamma)./(2*pi*Z_flat);
V_flat  = dF_flat ./ dw_flat;
F_flat  = A./Z_flat + conj(A).*Z_flat + (1i*Gamma/(2*pi))*log(Z_flat);

Psi_flat  = imag(F_flat);
Vmag_flat = abs(V_flat);
% Cap and clean before interpolation.
Vmag_flat(~isfinite(Vmag_flat)) = 0;
Vmag_flat(Vmag_flat > 3*U_inf)  = 3*U_inf;
ok = isfinite(Psi_flat) & isfinite(real(W_flat)) & isfinite(imag(W_flat));

wx_sc  = real(W_flat(ok));
wy_sc  = imag(W_flat(ok));
Psi_sc = Psi_flat(ok);
Vm_sc  = Vmag_flat(ok);

% -- Step 2: regular Cartesian grid in the w-plane --
nx = 350; ny = 300;
xg = linspace(-0.6, 1.6, nx);
yg = linspace(-0.9, 0.9, ny);
[Xg, Yg] = meshgrid(xg, yg);

% -- Step 3: scattered-to-regular interpolation --
% 'natural' (natural-neighbor) is smoother than 'linear' and handles the
% irregular point distribution from the log-spaced radial grid well.
try
    Psi_grid = griddata(wx_sc, wy_sc, Psi_sc, Xg, Yg, 'natural');
    Vm_grid  = griddata(wx_sc, wy_sc, Vm_sc,  Xg, Yg, 'natural');
catch
    warning('natural griddata failed, falling back to linear.');
    Psi_grid = griddata(wx_sc, wy_sc, Psi_sc, Xg, Yg, 'linear');
    Vm_grid  = griddata(wx_sc, wy_sc, Vm_sc,  Xg, Yg, 'linear');
end

% -- Step 4: mask airfoil interior --
airfoil_bx = [x_upper; flipud(x_lower)];
airfoil_by = [y_upper; flipud(y_lower)];
inside = inpolygon(Xg, Yg, airfoil_bx, airfoil_by);
Psi_grid(inside) = NaN;
Vm_grid(inside)  = NaN;

% =========================================================================
%  SECTION 6 -- FIGURES
% =========================================================================

patch_x = [x_upper; flipud(x_lower)];
patch_y = [y_upper; flipud(y_lower)];

% ------------------------------------------------------------------
%  FIGURE A: Geometry
% ------------------------------------------------------------------
figure('Name','Fig A: Airfoil Setup','NumberTitle','off');
plot(x_upper, y_upper, 'b-', 'LineWidth',2, 'DisplayName','Upper');
hold on;
plot(x_lower, y_lower, 'r-', 'LineWidth',2, 'DisplayName','Lower');
patch(patch_x, patch_y, [0.7 0.8 0.95], 'FaceAlpha',0.5,'EdgeColor','none');
quiver(-0.4, 0, 0.25*cos(aoa), 0.25*sin(aoa), 0, ...
       'k','LineWidth',2,'MaxHeadSize',0.5,'DisplayName','U_\infty');
axis equal; grid on;
xlim([-0.6 1.6]); ylim([-0.4 0.4]);
legend('Location','northeast');
title(sprintf('Fig A: NACA %s | AoA = %.1f°', naca_str, aoa_deg));
xlabel('x/c'); ylabel('y/c');

% ------------------------------------------------------------------
%  FIGURE B: Streamlines -- contour on regular (Xg,Yg) grid
% ------------------------------------------------------------------
figure('Name','Fig B: Streamlines','NumberTitle','off');
psi_levels = linspace(-2.5, 2.5, 80);
contour(Xg, Yg, Psi_grid, psi_levels, 'b', 'LineWidth',0.7);
hold on;
patch(patch_x, patch_y, [0.25 0.25 0.25], 'EdgeColor','none');
axis equal; grid on;
xlim([-0.6 1.6]); ylim([-0.9 0.9]);
title(sprintf('Fig B: Streamlines | NACA %s | AoA = %.1f°', naca_str, aoa_deg));
xlabel('x/c'); ylabel('y/c');

% ------------------------------------------------------------------
%  FIGURE C: Velocity magnitude -- contourf on regular (Xg,Yg) grid
% ------------------------------------------------------------------
figure('Name','Fig C: Velocity Field','NumberTitle','off');
v_levels = linspace(0, 2.5, 50);
contourf(Xg, Yg, Vm_grid, v_levels, 'LineColor','none');
colormap(jet); cb = colorbar;
cb.Label.String = '|V| / U_\infty';
hold on;
patch(patch_x, patch_y, [1 1 1], 'FaceAlpha',1.0,'EdgeColor',[0.3 0.3 0.3]);
axis equal; grid on;
xlim([-0.6 1.6]); ylim([-0.9 0.9]);
title(sprintf('Fig C: |V|/U_\\infty | NACA %s | AoA = %.1f°', naca_str, aoa_deg));
xlabel('x/c'); ylabel('y/c');

% ------------------------------------------------------------------
%  FIGURE D: Surface Cp
% ------------------------------------------------------------------
figure('Name','Fig D: Cp Distribution','NumberTitle','off');
plot(x_up_plot(valid_up), -Cp_upper(valid_up), 'b-o', ...
     'MarkerSize',4,'LineWidth',1.5,'DisplayName','Upper surface');
hold on;
plot(x_lo_plot(valid_lo), -Cp_lower(valid_lo), 'r-o', ...
     'MarkerSize',4,'LineWidth',1.5,'DisplayName','Lower surface');
yline(0,'k--','LineWidth',0.8);
set(gca,'YDir','reverse');
grid on; xlim([0 1]); ylim([-4 2]);
legend('Location','northeast');
xlabel('x/c'); ylabel('-C_p  (suction upward)');
title(sprintf('Fig D: Surface C_p | NACA %s | AoA = %.1f°', naca_str, aoa_deg));

% ------------------------------------------------------------------
%  FIGURE E: CL vs alpha
% ------------------------------------------------------------------
figure('Name','Fig E: CL vs Alpha','NumberTitle','off');
aoa_sweep = deg2rad(-10:1:20);
CL_sweep  = zeros(size(aoa_sweep));
for i = 1:length(aoa_sweep)
    Gi = 4 * pi * R_eff * U_inf * sin(aoa_sweep(i) - beta_TE);
    CL_sweep(i) = 2 * Gi / (U_inf * chord);
end

% Plot the SC calculated lift curve
plot(rad2deg(aoa_sweep), CL_sweep, 'b-', 'LineWidth', 2, 'DisplayName', 'SC Toolbox Data');
hold on;

% Add the theoretical Thin Airfoil Theory line for comparison (Slope = 2*pi)
p_fit        = polyfit(aoa_sweep, CL_sweep, 1);
alpha_L0_rad = -p_fit(2)/p_fit(1); % Zero-lift angle in radians
alpha_L0_deg = rad2deg(alpha_L0_rad);

% Theoretical CL = 2*pi * (alpha - alpha_L0)
CL_theory = 2 * pi * (aoa_sweep - alpha_L0_rad);
plot(rad2deg(aoa_sweep), CL_theory, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Theoretical 2\pi Slope');

% Plot current AoA dot
plot(aoa_deg, CL_current, 'ko', 'MarkerSize', 9, 'MarkerFaceColor', 'k', ...
    'DisplayName', sprintf('Current: CL=%.3f', CL_current));

xline(alpha_L0_deg,'k:','LineWidth',0.8, 'HandleVisibility', 'off');
xline(0,'k:','LineWidth',0.8, 'HandleVisibility', 'off'); 
yline(0,'k:','LineWidth',0.8, 'HandleVisibility', 'off');
grid on;
legend('Location','northwest');
xlabel('\alpha  [deg]'); ylabel('C_L');
title(sprintf('Fig E: Lift Curve | NACA %s | CL(%.1f°) = %.3f', ...
      naca_str, aoa_deg, CL_current));

% =========================================================================
%  SUMMARY
% =========================================================================
fprintf('\n--- Summary ---\n');
fprintf('  Airfoil  : NACA %s\n', naca_str);
fprintf('  AoA      : %.2f deg\n', aoa_deg);
fprintf('  Gamma    : %.4f\n', Gamma);
fprintf('  CL       : %.4f\n', CL_current);
fprintf('  alpha_L0 : %.2f deg\n', alpha_L0_deg);
fprintf('  Figures A-E generated.\n\n');

end % function airfoil_analysis