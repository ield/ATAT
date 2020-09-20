%% Author: ield
% This program calculates the patch antenna dimensions for a certain
% frequency and substrate. The program is designed for a circular patch
clear all;

% Modify these values
f = 1e9;                  % Central frequency [Hz]
e = 2.2;                  % Dielectric constant
h = 0.2e-3;               % Thickness of the substrate [m]
Rin = 50;                 % Input resistance at feeding point [ohm]

c = 3e8;                % Speed of light
n0 = 120*pi;            % Vacuum impedance

% First approximation of r (the patch length) taken from diapo 29
F = 8.794e9/(f*sqrt(e));
r = F / sqrt(1+2*h/(pi*e*F)*(log(pi*F/(2*h))+1.7726)) % Radius of the patch [cm]