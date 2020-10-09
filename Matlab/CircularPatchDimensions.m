%% Author: ield
% This program calculates the patch antenna dimensions for a certain
% frequency and substrate. The program is designed for a circular patch
clear all;

% Modify these values
f = 8.2e9;                % Central frequency [Hz]
e = 2.17;                 % Dielectric constant
h = 0.254e-3;                  % Thickness of the substrate [m] (probably)
Rin = 50;                 % Input resistance at feeding point [ohm]

c = 3e8;                % Speed of light
n0 = 120*pi;            % Vacuum impedance

% First approximation of L (the gnd and substrate length)
lambda0 = c / f;                % Wavelength
L = lambda0 / 2                 % Length of the substrate [m]


% First approximation of r (the patch length) taken from diapo 29
F = 8.794e9/(f*sqrt(e));
r = (F / sqrt(1+2*h/(pi*e*F)*(log(pi*F/(2*h))+1.7726)))/100 % Radius of the patch [m]
% r = 15.9e-3;

% Feeding point: set for Rin = 50m ohm
R0 = n0*lambda0 / (2*pi*2*r*(1-(1/6)*(pi*h/lambda0)^2));
L1 = 2*r/pi*acos(sqrt(Rin / R0)) % Feeding point        