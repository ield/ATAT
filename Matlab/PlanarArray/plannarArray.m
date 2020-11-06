%% Constant declaration
clear;
f = 60e9;       % Working frequency (Hz)
c = 3e8;        % Speed of light(m/s)
lambda = c/f;   % Wavelength (m)
D = 30;         % Directivity (dBi)
theta0 = 0;     % Pointting in broadside direction
res = 360;         % Resolution
path = '../../PlannarArray/Lab06/Images/';   % Path to save the files
%% Step 1: Calculate size of the antenna (diapo 26)
BW = 0.77*4*pi/10^(D/10);   % Beamwidth
% Since the antenna is a sqaure BWx = BWy
BWx = sqrt(BW);
BWy = sqrt(BW);

Lx = 0.88*lambda/BWx;       % Length in x
Ly = 0.88*lambda/BWy;       % Length in y: they are sqare so Lx = Ly

fprintf('\nDimensions of the array\n');
fprintf('Lx = %f (mm)     Ly = %f (mm)\n', Lx*1e3, Ly*1e3);

fprintf('\nThe bandwidth in the E plane BWy = %f�\n', BWy*180/pi);
fprintf('The bandwidth in the H plane BWx = %f�\n', BWx*180/pi);

%% Step 2: Calculate the number of elements
% Ideally, the separation between elements is lambda/2
dx = lambda/2;
dy = lambda/2;
M = ceil(Lx/dx);    % Number of elements in x dir. ceil so that the directivity is higher than the specification in stead of lower.
N = ceil(Ly/dy);    % Number of elements in y dir.
totalElem = N*M;    % Total number of elements

coorX = 0:dx:(M-1)*dx;   % Position of each element
coorY = 0:dx:(M-1)*dx;
[coorX, coorY] = meshgrid(coorX, coorY);

fprintf('\nThere are %i elements in the x direction separated %f mm.\n', M, dx*1e3);
fprintf('There are %i elements in the y direction separated %f mm.\n', N, dy*1e3);
fprintf('Therefore there are %i elements.\n', totalElem);

save PlanarArray.mat BWx BWy Lx Ly dx dy M N; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Array Factor Carthesians(Normalized)
% This is done only to plot the array factor
u = linspace(-1, 1, res);
v = linspace(-1, 1, res);
[u, v] = meshgrid(u,v);    

F = zeros(size(u));
phix = 2*pi*dx/lambda*u;    % Phase in x array factor function diapo 19
phiy = 2*pi*dy/lambda*v;    % Phase in y array factor function diapo 19

A = ones(M, N);
F = calcArrayFactor(A, M, N, phix, phiy, u);

printAndPlotArrayParameters(F, u, v, res, 'Array Factor', path, 'arrayFactor');
%% Effect of the element
thetaUseful = theta<pi/2;

E = thetaUseful.*(cos(theta).^2);

u = sin(theta).*cos(phi);
v = sin(theta).*sin(phi);
Ecar = abs(E).*cos(theta);

figure('Color',[1 1 1]);
plotSingleElement(u, v, Ecar, 0, 'Radiation element (lineal)');

saveas(gca, [path, 'elementField'],'epsc');
hold off;

% Combination element + array factor
% It is important to notice that the radiation diagram depends on the x and
% y axis, but it does not matter how the projection is done.
Etot = F.*E;

% Parameters BW, D0, SLL
Ecar = printAndPlotArrayParameters(Etot, u, v,res, 'Uniform amplitude', path, 'unifDis');
save unifromDistribution.mat A coorX coorY Ecar;
maxNorm = max(max(20*log10(abs(Ecar))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Radiation with failure of some elements (normalized)
failure = .25;      % Percentage of failing elements
A = rand(M);        % Feeding matrix
A = A > failure;    % Matrix with the working elements (1 works, 0 does not).
% plotDistribution(A, M, N)

% Calculate the arry factor and radiation diagram
F = calcArrayFactor(A, M, N, phix, phiy, theta);

% Multiply the array factor and the single element field to obtain the
% radiation pattern.
Etot = F.*abs(E);
Ecar = printAndPlotArrayParameters(Etot, u, v,res, 'Elements failing', path, 'failDis', maxNorm);
save unifromDistribution.mat A coorX coorY Ecar;
%% Radiation with triangular distribution 
A = triangDistribution(M, N);   % Triangular distribution
plotDistribution(A, M, N, path, 'triangDis');


% Calculate the arry factor and radiation diagram
F = calcArrayFactor(A, M, N, phix, phiy, theta);

% Multiply the array factor and the single element field to obtain the
% radiation pattern.
Etot = F.*abs(E);
Ecar = printAndPlotArrayParameters(Etot, u, v,res, 'Triangular distribution', path, 'triangPat', maxNorm);
save triangularDistribution.mat A coorX coorY Ecar;
%% Radiation with binomial distribution 
A = binomialDistribution(M, N);   % Triangular distribution
plotDistribution(A, M, N, path, 'binomDis')

% Calculate the arry factor and radiation diagram
F = calcArrayFactor(A, M, N, phix, phiy, theta);

% Multiply the array factor and the single element field to obtain the
% radiation pattern.
Etot = F.*E;
Ecar = printAndPlotArrayParameters(Etot, u, v,res, 'Binomial distribution', path, 'binomPat', maxNorm);
save binomialDistribution.mat A coorX coorY Ecar;
%% Radiation with cosine distribution 
A = cosDistribution(M, N);   % Triangular distribution
plotDistribution(A, M, N, path, 'cosDis')

% Calculate the arry factor and radiation diagram
F = calcArrayFactor(A, M, N, phix, phiy, theta);

% Multiply the array factor and the single element field to obtain the
% radiation pattern.
Etot = F.*abs(E);
Ecar = printAndPlotArrayParameters(Etot, u, v,res, 'Cosine distribution', path, 'cosPat', maxNorm);
save cosineDistribution.mat A coorX coorY Ecar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Beam Steering x, y
alphax = 45;    % Beam steering in �
alphay = 45;    % Beam stering in �


phix_st = phix + alphax*pi/180;
phiy_st = phiy + alphay*pi/180;

A = ones(M, N); % Uniform distribution

% Calculate the arry factor and radiation diagram
F = calcArrayFactor(A, M, N, phix_st, phiy_st, theta);

% Multiply the array factor and the single element field to obtain the
% radiation pattern.
Etot = F.*abs(E);
Ecar = printAndPlotSteering(Etot, u, v,'Beam Steering x, y', path, 'strxy', maxNorm);
save strxy.mat A coorX coorY Ecar;
%% Beam Steering theta, phi
theta0 = 20;    % Beam steering in �
phi0 = 20;    % Beam stering in �

[alphax, alphay] = beamSteering(theta0, phi0, lambda, dx, dy);  % Input in �, output in rad

phix_st = phix + alphax;
phiy_st = phiy + alphay;

A = ones(M, N); % Uniform distribution

% Calculate the arry factor and radiation diagram
F = calcArrayFactor(A, M, N, phix_st, phiy_st, theta);

% Multiply the array factor and the single element field to obtain the
% radiation pattern.
Etot = F.*E;
Ecar = printAndPlotSteering(Etot, u, v,'Beam Steering theta, phi', path, 'strthph', maxNorm);
save strthph.mat A coorX coorY Ecar;

