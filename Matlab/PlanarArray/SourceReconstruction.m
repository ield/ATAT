clear; close all;
path_image = '../../SourceReconstruction/Lab_08_Report/Images/';
path_antennas = 'AntennasSaved/';

% Load the main info
load PlanarArray.mat
% lambda = lambda*1e3;


% Load the file
load unifromDistribution.mat;
lambda = lambda;

% Define x and y axes according to the formula in slide 5
res = 360;
u = linspace(-1, 1, res);
v = linspace(-1, 1, res);
delta_u = u(2) - u(1);
delta_v = v(2) - v(1);
[u, v] = meshgrid(u,v);


delta_x = 1/(delta_u*res)*1e-3;
delta_y = 1/(delta_v*res)*1e-3;

% x = (-0.5*lambda/delta_u):delta_x:((0.5*lambda-delta_x)/delta_u);size(x)
% y = -0.5*lambda/delta_v:delta_y:(0.5*lambda-delta_y)/delta_v; size(y)
x = linspace(-0.5*lambda/delta_u, 0.5*lambda/delta_u, res);
y = linspace(-0.5*lambda/delta_v, 0.5*lambda/delta_v, res);
[x, y] = meshgrid(x, y);

% Define theta and phi fields
theta = real(asin(sqrt(u.^2 + v.^2)));
% theta(imag(theta) ~= 0) = pi/2;

% phi = asin(v./sin(theta));
% [theta, phi] = meshgrid(theta, phi);

phi = asin(v./sin(real(asin(sqrt(u.^2 + v.^2)))));
sin_phi = sin(phi);
E_theta = Etot.*sin_phi;

cos_phi = u./sin(theta);
E_phi = Etot.*cos_phi;

% Define the plane wave spectrum
P_theta = E_theta./cos(theta);
P_phi = E_phi./cos(theta);

P_x = P_theta.*cos(theta).*cos_phi - P_phi.*sin_phi;
P_y = P_theta.*cos(theta).*sin_phi - P_phi.*cos_phi;

% Transform the plane wave into the desired field
E_x = ifftshift(P_x);
E_y = ifftshift(P_y);


plot3Duv(x, y, 20*log10(abs(P_x)), -500, 'Initial radiated field');
figure;
plot3Duv(x, y, 20*log10(abs(P_y)), -500, 'Initial radiated field');