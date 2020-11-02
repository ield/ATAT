%% Constant declaration
clear;
f = 60e9;       % Working frequency (Hz)
c = 3e8;        % Speed of light(m/s)
lambda = c/f;   % Wavelength (m)
D = 30;         % Directivity (dBi)
theta0 = 0;     % Pointting in broadside direction

%% Step 1: Calculate size of the antenna (diapo 26)
BW = 0.77*4*pi/10^(D/10);   % Beamwidth
% Since the antenna is a sqaure BWx = BWy
BWx = sqrt(BW);
BWy = sqrt(BW);

Lx = 0.88*lambda/BWx;       % Length in x
Ly = 0.88*lambda/BWy;       % Length in y: they are sqare so Lx = Ly

fprintf('\nDimensions of the array\n');
fprintf('Lx = %f (mm)     Ly = %f (mm)\n', Lx*1e3, Ly*1e3);

fprintf('\nThe bandwidth in the E plane BWy = %fº\n', BWy*180/pi);
fprintf('The bandwidth in the H plane BWx = %fº\n', BWx*180/pi);

%% Step 2: Calculate the number of elements
% Ideally, the separation between elements is lambda/2
dx = lambda/2;
dy = lambda/2;
M = ceil(Lx/dx);    % Number of elements in x dir. ceil so that the directivity is higher than the specification in stead of lower.
N = ceil(Ly/dy);    % Number of elements in y dir.
totalElem = N*M;    % Total number of elements

fprintf('\nThere are %i elements in the x direction separated %f mm.\n', M, dx*1e3);
fprintf('There are %i elements in the y direction separated %f mm.\n', N, dy*1e3);
fprintf('Therefore there are %i elements.\n', totalElem);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Array Factor Carthesians(Normalized)
% This is done only to plot the array factor
res = 300;         % Resolution
u = linspace(-1, 1, res);
v = linspace(-1, 1, res);
[u, v] = meshgrid(u,v);    

F = zeros(size(u));
phix = 2*pi*dx/lambda*u;    % Phase in x array factor function diapo 19
phiy = 2*pi*dy/lambda*v;    % Phase in y array factor function diapo 19

A = ones(M, N);
F = calcArrayFactor(A, M, N, phix, phiy, u);

plot3Duv(u, v, abs(F), 0, 'Array Factor (lineal)');
plot3Duv(u, v, 20*log10(abs(F)), 0, 'Array Factor (dB)');
plot3Duv(u, v, 20*log10(abs(F/max(max(F)))), -30, 'Array Factor Normalized (dB)');

%% Array Factor Spherical(Normalized)
% It is necessary to do it in theta, phi, so that the sizes match the patch
% later
res = 500;         % Resolution
theta = linspace(0, pi, res);
phi = linspace(0,2*pi, res);
[theta, phi] = meshgrid(theta,phi); 

phix = 2*pi*dx/lambda*sin(theta).*cos(phi);    % Phase in x array factor function diapo 19
phiy = 2*pi*dy/lambda*sin(theta).*sin(phi);    % Phase in y array factor function diapo 19

A = ones(M, N);
F = calcArrayFactor(A, M, N, phix, phiy, theta);

% u = abs(F).*sin(theta).*cos(phi);
% v = abs(F).*sin(theta).*sin(phi);
Fcar = abs(F).*cos(theta);

% plot3Duv(u, v, Fcar, 0, 'Array Factor (lineal)');
% plot3Duv(u, v, 20*log10(abs(Fcar)), 0, 'Array Factor (dB)');
% plot3Duv(u, v, 20*log10(abs(Fcar/max(max(Fcar)))), -30, 'Array Factor Normalized (dB)');

%% Effect of the element Spherical
E = cos(theta).^2;

u = abs(E).*sin(theta).*cos(phi);
v = abs(E).*sin(theta).*sin(phi);
Ecar = abs(E).*cos(theta);

plot3Duv(u, v, Ecar, 0, 'Radiation element (lineal)');

%% Combination element + array factor
% It is important to notice that the radiation diagram depends on the x and
% y axis, but it does not matter how the projection is done.
Etot = F.*abs(E);

%% Parameters BW, D0, SLL
printAndPlotArrayParameters(u, v, Etot, phi, theta, res, 'Uniform amplitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Radiation with failure of some elements (normalized)
failure = .25;      % Percentage of failing elements
A = rand(M);        % Feeding matrix
A = A > failure;    % Matrix with the working elements (1 works, 0 does not).

% Calculate the arry factor and radiation diagram
F = calcArrayFactor(A, M, N, phix, phiy, theta);

% Multiply the array factor and the single element field to obtain the
% radiation pattern.
Etot = F.*abs(E);
printAndPlotArrayParameters(u, v, Etot, phi, theta, res, 'Elements failing');

%% Radiation with different feedings
A = triangDistribution(M, N);   % Triangular distribution

% Calculate the arry factor and radiation diagram
F = calcArrayFactor(A, M, N, phix, phiy, theta);

% Multiply the array factor and the single element field to obtain the
% radiation pattern.
Etot = F.*abs(E);
printAndPlotArrayParameters(u, v, Etot, phi, theta, res, 'Triangular distribution');





