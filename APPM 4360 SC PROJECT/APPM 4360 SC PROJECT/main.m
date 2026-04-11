%  APPM 4360 | Schwarz-Christoffel Mapping Project
%  File:    main.m
%  Authors: Jacob Johnson, Maxwell Meador, William Lampe
%
%  Purpose:
%    Entry point for the SC airfoil analysis pipeline.
%    Prompts the user for a NACA 4-digit designation and an angle of attack,
%    then hands control to airfoil_analysis.m which runs the full SC
%    mapping, flow field computation, and plotting.
%
%  Usage:
%    >> main
%    Enter a 4-digit NACA number (default = 2412): 2412
%    Enter angle of attack in degrees (default = 5): 5
%
%  Dependencies:
%    airfoil_analysis.m   -- main analysis driver
%    naca_airfoil.m       -- airfoil geometry generator
%    SC Toolbox           -- polygon(), extermap(), prevertex(), evaldiff()
%                           (Driscoll & Trefethen, must be on MATLAB path)
%
%  Notes:
%    - Run  >>  addpath(genpath('YOUR_SC_TOOLBOX_PATH'))  before this script
%      if the SC Toolbox is not already on your permanent path.
%    - Tested on MATLAB R2022b.  No Toolboxes beyond SC required.

clear;      % clear workspace variables
clc;        % clear command window
close all;  % close any open figure windows

% -------------------------------------------------------------------------
%  1.  NACA designation
%      input(..., 's') reads the answer as a string, so leading zeros are
%      preserved (e.g. '0012' stays '0012', not the number 12).
% -------------------------------------------------------------------------
naca_str = input('Enter a 4-digit NACA number (default = 2412): ', 's');

% If the user just hit Enter, use the symmetric NACA 2412 default.
if isempty(naca_str)
    naca_str = '2412';
end

% Basic sanity check: must be exactly 4 digit characters.
if length(naca_str) ~= 4 || any(~isstrprop(naca_str, 'digit'))
    error('Please enter exactly 4 digits, e.g.  2412  or  0012.');
end

% -------------------------------------------------------------------------
%  2.  Angle of attack
%      input() with no 's' flag returns a numeric value (or [] on empty).
% -------------------------------------------------------------------------
aoa_deg = input('Enter angle of attack in degrees (default = 5): ');

% If the user just hit Enter, use 5 degrees.
if isempty(aoa_deg)
    aoa_deg = 5;
end

% Warn (but do not crash) for angles well outside the attached-flow regime.
% The SC mapping is mathematically valid for any AoA, but the Kutta
% condition assumption breaks down above ~15 deg on a real airfoil.
if abs(aoa_deg) > 20
    warning(['AoA = %.1f deg is outside the typical attached-flow range. ' ...
             'CL results will be optimistic (no stall modeled).'], aoa_deg);
end

% -------------------------------------------------------------------------
%  3.  Run the analysis
% -------------------------------------------------------------------------
fprintf('\nStarting SC airfoil analysis...\n');
airfoil_analysis(naca_str, aoa_deg);

fprintf('\n--- Done. All figures generated. ---\n');
