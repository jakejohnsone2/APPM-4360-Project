%  APPM 4360 | Schwarz-Christoffel Mapping Project
%  Function Name: airfoil_analysis.m 
%  Authors: Maxwell Meador
%  Purpose: Calculates a ton of aerodynamic parameters using the SC map
%  (Inviscid, Irrotational, Incompressible flow)
%  Inputs: NACA digits, AOA (main.m),  Sc mapped airfoil polygon
%  (sc_map_polygon)
%  Outputs:
%    Figure A: Airfoil + angle of attack schematic
%    Figure B: Streamlines around the airfoil
%    Figure C: Velocity magnitude field
%    Figure D: Pressure coefficient distribution on airfoil surface
%    Figure E: Lift coefficient vs angle of attack comparison


function airfoil_analysis(naca_str, aoa_deg)
    if nargin < 1, naca_str = '2412'; end
    if nargin < 2, aoa_deg  = 5;      end
    
    aoa = deg2rad(aoa_deg);
    U_inf = 1.0; % Freestream velocity magnitude
    
    fprintf('\n NACA %s at %.1f° AoA\n', naca_str, aoa_deg);
    
    % Generate airfoil 
    N_surface = 150; 
    [x_upper, y_upper, x_lower, y_lower, x_camber, y_camber, poly_z, alphas_interior, params] = naca_airfoil(naca_str, N_surface, true);
    
    % CONVERT TO EXTERIOR ANGLES (phis)
    phis = 1 - alphas_interior;
    
    % Run the Schwarz-Christoffel Solver
    [z_of_t, t_prevert, C_const] = sc_map_polygon(poly_z, phis);
    
    % Airfoil map
    t_surf = 0.5 * (t_prevert(1:end-1) + t_prevert(2:end));
    
    % Circulation for Kutta Condition
    alpha_ZL = -2 * params.M * (0.5 - params.P);        
    Gamma = 4 * pi * U_inf * sin(aoa - alpha_ZL);
    
    % SC Derivative
    dz_dt_surf = zeros(size(t_surf));
    for j = 1:length(t_surf)
        dz_dt_surf(j) = C_const * prod((t_surf(j) - t_prevert).^(-phis(:)));
    end
    
    % Flat velocity
    [~, TE_idx] = min(abs(poly_z - 1.0));
    t_TE = t_prevert(TE_idx);
    dw_dt_surf = U_inf * exp(-1i*aoa) + (1i*Gamma)./(2*pi*(t_surf - t_TE + 1e-6));
    
    % V = (dw/dt) / (dz/dt)
    V_surf = dw_dt_surf ./ (C_const * dz_dt_surf);
    Cp = 1 - (abs(V_surf) / U_inf).^2;
    
    % Streamlines (Help)
    t_x = linspace(min(t_prevert)-0.2, max(t_prevert)+0.5, 35);
    t_y = linspace(0.01, 0.8, 25);
    [T_real, T_imag] = meshgrid(t_x, t_y);
    T_grid = T_real + 1i*T_imag;
    
    Z_grid = zeros(size(T_grid));
    V_grid = zeros(size(T_grid));
    Psi_grid = zeros(size(T_grid));
    
    for r = 1:size(T_grid, 1)
        for c = 1:size(T_grid, 2)
            t_curr = T_grid(r,c);
            
            % t -> z
            Z_grid(r,c) = z_of_t(t_curr);
            
            W = U_inf * exp(-1i*aoa) * t_curr + (1i * Gamma / (2*pi)) * log(t_curr - t_TE);
            Psi_grid(r,c) = imag(W);
            
            % SC Derivative
            dz_dt = C_const * prod((t_curr - t_prevert).^(-phis));
            dw_dt      = U_inf * exp(-1i*aoa) + (1i*Gamma) ./(2*pi*(t_curr - t_TE));
            V_grid(r,c) = dw_dt / dz_dt;
        end
    end
    
    %% Plotting (Thank u gemini)
    % FIGURE A: Airfoil + Setup 
    fig_A = figure('Name','Airfoil Setup','NumberTitle','off');
    fig_A.Position = [50 550 750 350];
    plot(x_upper, y_upper, 'b-', 'LineWidth', 2); hold on;
    plot(x_lower, y_lower, 'b-', 'LineWidth', 2);
    patch([x_upper; flipud(x_lower)], [y_upper; flipud(y_lower)], [0.7 0.8 0.95], 'FaceAlpha', 0.5);
    quiver(-0.4, 0, 0.2*cos(aoa), 0.2*sin(aoa), 0, 'r', 'LineWidth', 2);
    axis equal; grid on; title(['NACA ', naca_str, ' Mapping Setup']);
    
    % FIGURE B: Streamlines
    fig_B = figure('Name','Streamlines','NumberTitle','off');
    fig_B.Position = [100 50 900 500];
    contour(real(Z_grid), imag(Z_grid), Psi_grid, 40, 'b'); hold on;
    patch([x_upper; flipud(x_lower)], [y_upper; flipud(y_lower)], [0.3 0.3 0.3]);
    axis equal; xlim([-0.2 1.2]); ylim([-0.5 0.5]); grid on;
    title('Figure B: Streamlines (via SC Mapping)');
    
    % FIGURE C: Velocity Magnitude
    fig_C = figure('Name','Velocity Field','NumberTitle','off');
    fig_C.Position = [200 50 900 500];
    contourf(real(Z_grid), imag(Z_grid), abs(V_grid), 40, 'LineColor','none');
    colormap(jet); colorbar; hold on;
    patch([x_upper; flipud(x_lower)], [y_upper; flipud(y_lower)], [1 1 1], 'FaceAlpha', 0.5);
    axis equal; xlim([-0.2 1.2]); ylim([-0.5 0.5]);
    title('Figure C: Velocity Magnitude |V|/U_\infty');
    
    % FIGURE D: Pressure Coefficient
    fig_D = figure('Name','Cp Distribution','NumberTitle','off');
    fig_D.Position = [300 50 800 500];
    half = length(x_upper);
    plot(x_upper, -Cp(1:half), 'b-', 'LineWidth', 2); hold on;
    plot(x_lower, -Cp(end-length(x_lower)+1:end), 'r-', 'LineWidth', 2);
    set(gca, 'YDir','reverse'); grid on; 
    xlabel('x/c'); ylabel('-Cp (Suction Up)');
    title('Figure D: Surface Pressure Distribution');
    
    % FIGURE E: CL vs AoA
    fig_E = figure('Name','CL vs Alpha','NumberTitle','off');
    fig_E.Position = [400 50 700 500];
    aoa_sweep = deg2rad(-10:2:20);
    CL_sweep = 8 * pi * sin(aoa_sweep - alpha_ZL) / params.chord;
    plot(rad2deg(aoa_sweep), CL_sweep, 'b-', 'LineWidth', 2); hold on;
    plot(aoa_deg, 2*Gamma/(U_inf*params.chord), 'ko', 'MarkerFaceColor','k');
    grid on; xlabel('Alpha [deg]'); ylabel('C_L');
    title('Figure E: Lift Curve (Predicted by SC Morphing)');
end