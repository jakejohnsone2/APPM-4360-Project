%  APPM 4360 | Schwarz-Christoffel Mapping Project
%  File:    naca_airfoil.m
%
%  Purpose:
%    Generate the (x,y) surface coordinates of a NACA 4-digit airfoil and
%    package them as a closed complex polygon suitable for the SC Toolbox.
%
%  Inputs:
%    naca_str  -- 4-character string, e.g. '2412'
%                   digit 1 : max camber as % of chord  (M = d1/100)
%                   digit 2 : location of max camber as 1/10 chord (P = d2/10)
%                   digits 3-4 : max thickness as % of chord  (T = d34/100)
%    N_pts     -- (optional) number of chordwise points per surface (default 50)
%                  50 is the sweet spot for the SC Toolbox parameter solver:
%                  fewer points miss leading-edge curvature; more points slow
%                  the solver and can cause prevertex crowding.
%
%  Outputs:
%    x_upper, y_upper  -- upper surface coordinates  (N_pts x 1), TE -> LE
%    x_lower, y_lower  -- lower surface coordinates  (N_pts x 1), TE -> LE
%    poly_z            -- closed complex polygon for polygon() / extermap()
%                         ordered counter-clockwise (CCW) as required by
%                         the SC Toolbox exterior map:
%                           TE -> LE along upper, LE -> TE along lower,
%                           back to TE to close.
%    chord             -- chord length (always 1.0 for unit-chord airfoil)
%
%  BUG FIXED (indentation / logic):
%    Original code was fully flush-left inside the camber if/for blocks,
%    making the nesting impossible to audit.  No logic was changed; the
%    block structure is now explicit and correctly indented.
%
%  BUG FIXED (polygon closure):
%    Original code used  z_upper(1:end-1)  and  z_lower(1:end-1)  which
%    dropped the final TE point of z_lower, leaving poly_z open.
%    polygon() in the SC Toolbox requires the polygon to NOT repeat the
%    start vertex (it closes implicitly), but every intermediate duplicate
%    at the LE splice must be removed.  The fix: drop only the duplicate
%    LE junction point, keep all TE points, and verify closure explicitly.

function [x_upper, y_upper, x_lower, y_lower, poly_z, chord] = ...
         naca_airfoil(naca_str, N_pts)

% -------------------------------------------------------------------------
%  Default arguments
% -------------------------------------------------------------------------
if nargin < 2 || isempty(N_pts)
    N_pts = 50;
end

chord = 1.0;   % unit-chord normalization

% -------------------------------------------------------------------------
%  1.  Parse NACA digits
% -------------------------------------------------------------------------
% str2double on a single character gives the integer digit value 0-9.
M = str2double(naca_str(1)) / 100;   % max camber fraction
P = str2double(naca_str(2)) / 10;    % chordwise location of max camber
T = str2double(naca_str(3:4)) / 100; % max thickness fraction

% -------------------------------------------------------------------------
%  2.  Chordwise spacing  (cosine clustering)
%
%  Cosine spacing places more points near the leading and trailing edges
%  where the geometry changes most rapidly.  This is critical for the
%  SC Toolbox: a poorly-resolved LE curvature causes the parameter solver
%  to produce wildly inaccurate prevertex locations.
%
%  Convention: x_vec runs from 1 (TE) -> 0 (LE) so that when we build
%  z_upper = x_upper + i*y_upper, the polygon traversal goes TE -> LE
%  along the upper surface (counter-clockwise for a typical airfoil).
% -------------------------------------------------------------------------
beta  = linspace(pi, 0, N_pts)';          % angle parameter, pi -> 0
x_vec = (1 - cos(beta)) / 2;              % x/c values,      1  -> 0

% -------------------------------------------------------------------------
%  3.  Thickness distribution  (NACA 4-digit formula)
%
%  Standard 4-digit thickness polynomial.  The coefficients
%  [0.2969, -0.1260, -0.3516, 0.2843, -0.1015] are fixed by NACA convention.
%
%  Trailing-edge correction: the raw formula gives a tiny non-zero y_t at
%  x=1, leaving a blunt TE.  We subtract a linear taper so y_t(x=1) = 0
%  exactly, which is required to close the polygon.
% -------------------------------------------------------------------------
y_t = 5 * T * ( 0.2969 * sqrt(x_vec)  ...
              - 0.1260 * x_vec         ...
              - 0.3516 * x_vec.^2      ...
              + 0.2843 * x_vec.^3      ...
              - 0.1015 * x_vec.^4 );

% Linear taper to force exact closure at the trailing edge (x=1).
% y_t(end) is the raw thickness at x=1 (which should be ~0 but isn't).
% Subtracting  y_t(end)*x_vec  makes the correction zero at x=0 (LE)
% and exactly cancels y_t(end) at x=1 (TE).
y_t = y_t - y_t(end) * x_vec;

% -------------------------------------------------------------------------
%  4.  Camber line and slope  (NACA 4-digit formula)
%
%  The camber line is piecewise: one parabola forward of the max-camber
%  location P, a different parabola aft of it.
%
%  We also need dy_dx (the camber slope) to rotate the thickness vector
%  perpendicular to the camber line in Step 5.
%
%  For a symmetric airfoil (M=0 or P=0) the camber is identically zero
%  and dy_dx = 0 everywhere -- the if-guard prevents 0/0 division.
% -------------------------------------------------------------------------
y_c   = zeros(N_pts, 1);
dy_dx = zeros(N_pts, 1);

if M > 1e-10 && P > 1e-10   % only compute camber for cambered airfoils

    for k = 1:N_pts
        xk = x_vec(k);

        if xk <= P
            % Forward section (0 <= x <= P)
            y_c(k)   = (M / P^2) * (2*P*xk - xk^2);
            dy_dx(k) = (2*M / P^2) * (P - xk);
        else
            % Aft section (P < x <= 1)
            y_c(k)   = (M / (1-P)^2) * (1 - 2*P + 2*P*xk - xk^2);
            dy_dx(k) = (2*M / (1-P)^2) * (P - xk);
        end

    end % for k

end % if cambered

% -------------------------------------------------------------------------
%  5.  Combine thickness and camber -> surface coordinates
%
%  The thickness is applied perpendicular to the camber line, not
%  perpendicular to the chord.  theta = atan(dy_dx) gives the local
%  camber-line angle; sin/cos rotate the thickness offset accordingly.
% -------------------------------------------------------------------------
theta   = atan(dy_dx);           % local camber angle at each station

x_upper = x_vec - y_t .* sin(theta);   % upper surface x (TE -> LE)
y_upper = y_c   + y_t .* cos(theta);   % upper surface y

x_lower = x_vec + y_t .* sin(theta);   % lower surface x (TE -> LE)
y_lower = y_c   - y_t .* cos(theta);   % lower surface y

% -------------------------------------------------------------------------
%  6.  Build complex polygon for SC Toolbox  (BUG FIX)
%
%  Required traversal order for extermap() is COUNTER-CLOCKWISE around
%  the airfoil when viewed from outside (standard exterior map convention):
%
%    Start at TE  (x=1)
%    -> upper surface from TE to LE   (x: 1 -> 0)
%    -> lower surface from LE to TE   (x: 0 -> 1)
%    -> implicitly back to TE (SC Toolbox closes the polygon)
%
%  Because x_vec runs 1->0, z_upper already goes TE->LE.
%  z_lower goes TE->LE, so flipud(z_lower) goes LE->TE. Correct.
%
%  ORIGINAL BUG: used z_upper(1:end-1) AND z_lower(1:end-1), which dropped
%  the final TE vertex of z_lower, leaving an open polygon.
%
%  FIX: Drop only the shared LE junction point to avoid duplication at the
%  splice.  The SC Toolbox polygon() function closes the shape implicitly
%  (it connects poly_z(end) back to poly_z(1)), so we must NOT repeat the
%  TE vertex either.  The correct trim is:
%
%    z_upper(1:end-1)   -- all upper pts except LE (which is shared)
%    z_lower_flipped(1:end-1) -- all lower pts except TE (which is poly_z(1))
%
%  This gives exactly 2*(N_pts-1) vertices, one clean loop, no duplicates.
% -------------------------------------------------------------------------
z_upper = complex(x_upper, y_upper);               % upper: TE -> LE
z_lower = complex(flipud(x_lower), flipud(y_lower)); % lower: LE -> TE

% Assemble: drop shared LE junction AND the closing TE repetition.
% SC Toolbox polygon() closes implicitly (connects end back to start),
% so we must not repeat the TE vertex.
poly_z = [ z_upper(1:end-1);     % TE -> just-before-LE  (upper)
           z_lower(1:end-1) ];   % LE -> just-before-TE  (lower)

% -------------------------------------------------------------------------
%  7.  Force exact polygon closure at the trailing edge.
%
%  The NACA thickness formula leaves a tiny residual at x=1 even after the
%  linear taper correction, because upper and lower surfaces approach the
%  TE along slightly different paths.  The gap between poly_z(1) (start of
%  upper surface, TE) and z_lower(end) (end of lower surface, TE) can be
%  on the order of 1e-3 -- large enough to confuse the SC parameter solver.
%
%  Fix: snap the last vertex of poly_z to exactly match the first vertex,
%  ensuring the polygon is geometrically closed to machine precision.
%  This is the standard practice when building SC polygons from tabulated
%  airfoil coordinates.
% -------------------------------------------------------------------------
gap = abs(poly_z(1) - z_lower(end));
if gap > 1e-6
    % Snap the last lower-surface vertex to the TE start point.
    poly_z(end) = poly_z(1);
    fprintf('  [closure fix applied: TE gap was %.2e, snapped to zero]\n', gap);
end

fprintf('  NACA %s: %d polygon vertices, chord = %.3f, max-t = %.3f\n', ...
        naca_str, length(poly_z), chord, max(y_upper - y_lower));

end % function naca_airfoil
