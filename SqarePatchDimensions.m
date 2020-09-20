%% Author: ield
% This program calculates the patch antenna dimensions for a certain
% frequency and substrate. The program is designed for a square patch
clear all;

% Modify these values
f = 3.5e9;                % Central frequency [Hz]
e = 2.33;                 % Dielectric constant
% h = 0.2e-3;               % Thickness of the substrate [m]: when given comment line 10, 26, 30
bw = 200e6;               % Bandwidth [MHz]: when given, comment line 9, 25
Rin = 50;                 % Input resistance at feeding point [ohm]

c = 3e8;                % Speed of light
n0 = 120*pi;            % Vacuum impedance

% First approximation of L (the patch length)
lambda0 = c / f;                % Wavelength
lambdaG =  lambda0 / sqrt(e);   % Wavelength in the medium
Lpatch = lambdaG / 2            % Patch length [m]
L = lambda0 / 2                 % Length of the substrate [m]

% Feeding point: set for Rin = 50m ohm
%Using the formulas in page 27, first it is needed to know the thickness
%Using the formula in page 11
% bw = 16/(3*sqrt(2))*1/e*h/lambda0 % Bandwidth of the substrate. Since the patch is a square, w/L = 1. Formula taken from diapo 11.
bw = bw / 1e10;                     % The bw is transformed to the appropriate units. Conversion obtained empirically
h = bw*3*sqrt(2)/16*e*lambda0       % Thickness of the subtrate [m]
R0 = n0*lambda0 / (2*pi*bw*(1-(1/6)*(pi*h/lambda0)^2));
L1 = Lpatch/pi*acos(sqrt(Rin / R0)) % Feeding point         

% Othe parameters 
% eEf = (e+1) / 2 + (e - 1) / (2 * sqrt(1 + 12 * h / L));
