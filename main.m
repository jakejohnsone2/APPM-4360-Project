%  APPM 4360 | Schwarz-Christoffel Mapping Project
%  Function Name: main.m
%  Authors: Jacob Johnson, Maxwell Meador, William Lampe
%  Purpose: Entry point. Asks the user for a NACA number and angle of
%           attack, then runs the full SC airfoil analysis.
%  Inputs:  User typed input (NACA 4-digit number, angle of attack)
%  Outputs: Calls airfoil_analysis.m which saves all figures

clear; clc; close all;

% Get NACA number from user
naca_str = input('Enter a 4-digit NACA number (default = 2412): ', 's');
if isempty(naca_str)
    naca_str = '2412';
end

% Get angle of attack from user
aoa_deg = input('Enter angle of attack in degrees (default = 5): ');
if isempty(aoa_deg)
    aoa_deg = 5;
end

airfoil_analysis(naca_str, aoa_deg);

heatsink_analysis(); 

fprintf('\nDone\n');
