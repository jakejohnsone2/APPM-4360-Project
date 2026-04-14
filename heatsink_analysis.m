function heatsink_analysis()
    %  APPM 4360 | Schwarz-Christoffel Mapping Project
    %  Purpose: Standalone Heat Sink Thermodynamic Analysis
    
    fprintf('\n=== Heat Sink Thermal Analysis ===\n');
    
    %% Geometry and Mapping
    fprintf('[1/2] Generating Heat Sink Mapping...\n');
    
    hs_verts = [0, 3, 3+1i, 2.5+1i, 2.5+2.5i, 1.5+2.5i, 1.5+1i, ...
                0.5+1i, 0.5+2.5i, -0.5+2.5i, -0.5+1i, 0+1i].';
    p_hs = polygon(hs_verts);
    f_hs = hplmap(p_hs); % Half-Plane Map for the interior
    
    %% Heat Sink Thermo
    fprintf('[2/2] Solving Temperature Distribution...\n');
    figure('Name', 'Heat Sink Thermodynamics', 'Position', [100, 400, 800, 600]);
    
    w_pre = prevertex(f_hs);
    w1 = real(w_pre(1)); 
    w2 = real(w_pre(2)); 
    a = min(w1, w2);
    b = max(w1, w2);
    
    w_pre_finite = w_pre(~isinf(w_pre));
    W_pts = [];
    % Generate sampling points in the canonical plane
    for i = 1:length(w_pre_finite)
        [R_rad, TH_rad] = meshgrid(logspace(-4, 1, 40), linspace(0.01, pi-0.01, 30));
        W_pts = [W_pts; w_pre_finite(i) + R_rad(:).*exp(1i*TH_rad(:))];
    end
    [U_bg, V_bg] = meshgrid(linspace(min(w_pre_finite)-3, max(w_pre_finite)+3, 100), logspace(-3, 2, 80));
    W_pts = [W_pts; U_bg(:) + 1i*V_bg(:)];
    
    % Analytical Temperature in canonical plane
    T_pts = (30 / pi) * abs(atan2(imag(W_pts), real(W_pts) - a) - atan2(imag(W_pts), real(W_pts) - b));
    Z_pts = f_hs(W_pts);
    
    % Interpolation to Physical Plane
    valid = isfinite(Z_pts) & isfinite(T_pts);
    Z_clean = Z_pts(valid); 
    T_clean = T_pts(valid);
    [~, unq_idx] = unique(round(real(Z_clean)*1e4) + 1i*round(imag(Z_clean)*1e4));
    F_temp = scatteredInterpolant(real(Z_clean(unq_idx)), imag(Z_clean(unq_idx)), T_clean(unq_idx), 'linear', 'nearest');
    
    [Xq_hs, Yq_hs] = meshgrid(linspace(-1, 4, 300), linspace(0, 3, 300));
    Tq_hs = F_temp(Xq_hs, Yq_hs);
    in_hs = inpolygon(Xq_hs, Yq_hs, real(hs_verts), imag(hs_verts));
    Tq_hs(~in_hs) = NaN;
    
    % Plotting
    contourf(Xq_hs, Yq_hs, Tq_hs, 40, 'LineColor', 'none'); hold on;
    colormap(jet); colorbar; caxis([0 30]);
    plot(p_hs, 'k', 'LineWidth', 2);
    title('Heat Sink: Temperature Distribution (Base=30^\circC)', 'FontSize', 14);
    xlabel('x'); ylabel('y');
    axis equal; grid on; xlim([-1, 4]); ylim([-0.5, 3]);
    
    fprintf('Heat Sink analysis complete.\n');
end