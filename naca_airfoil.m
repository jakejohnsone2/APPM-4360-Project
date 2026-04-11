%  APPM 4360 | Schwarz-Christoffel Mapping Project
%  Function Name: naca_airfoil.m
%  Authors: Maxwell Meador
%  Purpose: Generates (x,y) coords of a NACA airfoil shape given 4 digits.
%           Also discretizes it into a polygon for the SC map, and computes
%           the INTERIOR angle
%  Inputs:  naca_str - 4-digit string like '2412'
%           N_pts    - points per surface (default 100)
%           close_te - whether to force TE to close exactly (default true)
%  Outputs:
%    x_upper, y_upper  - upper surface coords (LE to TE)
%    x_lower, y_lower  - lower surface coords (LE to TE)
%    x_camber, y_camber - camber line coords
%    poly_z   - full airfoil boundary as complex polygon (CCW)
%    alphas   - interior angle / pi at every polygon vertex
%    params   - struct with M, P, T, chord, alpha_TE, etc.

function [x_upper, y_upper, x_lower, y_lower, x_camber, y_camber, poly_z, alphas, params] = ...
         naca_airfoil(naca_str, N_pts, close_te)

    % Check
    if length(naca_str) ~= 4 || ~all(isdigit(naca_str))
        error('make sure that you used 4 numbers');
    end

    M = str2double(naca_str(1)) / 100;    % max camber
    P = str2double(naca_str(2)) / 10;     % camber position
    T = str2double(naca_str(3:4)) / 100;  % max thickness

    if nargin < 2 || isempty(N_pts),    N_pts    = 100;  end
    if nargin < 3 || isempty(close_te), close_te = true; end

    params.M        = M;
    params.P        = P;
    params.T        = T;
    params.naca_str = naca_str;
    params.N_pts    = N_pts;
    params.chord    = 1.0;

    fprintf('\n[NACA] Generating NACA %s\n', naca_str);
    if M < 1e-10
        fprintf('[NACA] Symmetric airfoil (no camber)\n');
    end

    % Cosine spacing
    beta  = linspace(pi, 0, N_pts)';
    x_vec = (1 - cos(beta)) / 2;

    % thickness
    y_t = 5 * T * (  0.2969 * sqrt(x_vec) - 0.1260 * x_vec  - 0.3516 * x_vec.^2  + 0.2843 * x_vec.^3 - 0.1015 * x_vec.^4 );

    if close_te
        y_t = y_t - y_t(end) * x_vec;
    end

    % Camber line
    y_c   = zeros(N_pts, 1);
    dy_dx = zeros(N_pts, 1);

    if M < 1e-10 || P < 1e-10
        % Symmetric
        y_c(:)   = 0;
        dy_dx(:) = 0;
    else
        for k = 1:N_pts
            xk = x_vec(k);
            if xk <= P
                y_c(k)   = (M / P^2) * (2*P*xk - xk^2);
                dy_dx(k) = (2*M / P^2) * (P - xk);
            else
                y_c(k)   = (M / (1-P)^2) * (1 - 2*P + 2*P*xk - xk^2);
                dy_dx(k) = (2*M / (1-P)^2) * (P - xk);
            end
        end
    end

    theta = atan(dy_dx);

    x_upper = x_vec - y_t .* sin(theta);
    y_upper = y_c  + y_t .* cos(theta);

    x_lower = x_vec + y_t .* sin(theta);
    y_lower = y_c  - y_t .* cos(theta);

    x_camber = x_vec;
    y_camber  = y_c;

    % Build polygon (US is LE to TE LS is TE to LE)
    z_upper = complex(x_upper, y_upper);
    z_lower = complex(flipud(x_lower), flipud(y_lower)); %flipud just flips the order
    poly_z  = [z_upper(1:end-1); z_lower(1:end-1)];

    n_poly = length(poly_z);
    fprintf('[NACA] Polygon has %d vertices total.\n', n_poly);

   
% angles for SC map
    alphas = zeros(n_poly, 1);

    for k = 1:n_poly
        % Indices of previous and next vertex
        k_prev = mod(k - 2, n_poly) + 1;
        k_next = mod(k,     n_poly) + 1;

        v_in  = poly_z(k)      - poly_z(k_prev);  
        v_out = poly_z(k_next) - poly_z(k);

        turn = angle(v_out / v_in);

        interior_angle = pi - turn;

        alphas(k) = interior_angle / pi;
    end

    angle_sum = sum(alphas);

    params.alpha_TE  = alphas(N_pts - 1); 
    params.alpha_sum = angle_sum;


    % Plot airfoil (Thanks gemini)
    fig = figure('Name', sprintf('NACA %s Geometry', naca_str), 'NumberTitle','off');
    fig.Position = [150 400 900 400];

    plot(x_upper, y_upper, 'b-', 'LineWidth', 2); hold on;
    plot(x_lower, y_lower, 'r-', 'LineWidth', 2);
    plot(x_camber, y_camber, 'k--', 'LineWidth', 1.2);
    plot([0 1], [0 0], 'k:', 'LineWidth', 0.8);

    plot(0, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor','k');
    text(0.01, 0.02, 'LE', 'FontSize', 10);
    plot(1, 0, 'k^', 'MarkerSize', 8, 'MarkerFaceColor','k');
    text(0.98, -0.03, 'TE', 'FontSize', 10, 'HorizontalAlignment','right');

    [~, ip] = max(y_upper);
    plot(x_upper(ip), y_upper(ip), 'b^', 'MarkerSize', 7, 'MarkerFaceColor','b');
    text(x_upper(ip), y_upper(ip)+0.015, sprintf('Max thickness %.1f%%', T*100), ...
        'FontSize', 8, 'Color','b', 'HorizontalAlignment','center');

    patch([x_upper; flipud(x_lower)], [y_upper; flipud(y_lower)], ...
        [0.7 0.8 0.95], 'FaceAlpha', 0.6, 'EdgeColor','none');

    axis equal; grid on;
    xlim([-0.05 1.1]); ylim([-0.25 0.25]);
    xlabel('x/c (chord fraction)');
    ylabel('y/c');
    title(sprintf('NACA %s  |  M=%.2f  P=%.1f  T=%.2f  |  %d polygon vertices', ...
        naca_str, M, P, T, n_poly));
    legend({'Upper surface','Lower surface','Camber line','Chord line'}, ...
        'Location','northeast','FontSize',8);
    saveas(fig, sprintf('naca%s_geometry.png', naca_str));
    fprintf('[NACA] Geometry figure saved.\n');
end


% for checking digits w user inpuy
function result = isdigit(str)
    result = all(str >= '0' & str <= '9');
end
