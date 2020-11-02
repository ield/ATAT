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

for ii = 0:M-1        % For all the elements in the x direction
    for jj = 0:N-1    % For all the elements in the y direction
        F = F + exp(1i*(ii*phix + jj*phiy));
    end
end

plot3Duv(u, v, abs(F), 0, 'Array Factor (lineal)');
plot3Duv(u, v, 20*log10(abs(F)), 0, 'Array Factor (dB)');
plot3Duv(u, v, 20*log10(abs(F/max(max(F)))), -30, 'Array Factor Normalized (dB)');

%% Array Factor Spherical(Normalized)
% It is necessary to do it in theta, phi, so that the sizes match the patch
% later
res = 300;         % Resolution
theta = linspace(0, pi, res);
phi = linspace(0,2*pi, res);
[theta, phi] = meshgrid(theta,phi); 

F = zeros(size(theta));
phix = 2*pi*dx/lambda*sin(theta).*cos(phi);    % Phase in x array factor function diapo 19
phiy = 2*pi*dy/lambda*sin(theta).*sin(phi);    % Phase in y array factor function diapo 19

for ii = 0:M-1        % For all the elements in the x direction
    for jj = 0:N-1    % For all the elements in the y direction
        F = F + exp(1i*(ii*phix + jj*phiy));
    end
end

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
[~, ~, Etotcar] = sph2cart(phi, pi/2-theta, abs(Etot));
plot3Duv(u, v, abs(Etotcar), 0, 'Radiation element (lineal)');
plot3Duv(u, v, 20*log10(abs(Etotcar)), 0, 'Array Factor (dB)');
plot3Duv(u, v, 20*log10(abs(Etotcar/max(max(Etotcar)))), -30, 'Array Factor Normalized (dB)');

%% Parameters BW, D0, SLL
% First it is calculated the normalized radiation diagram in planes E and H
% In phi = pi/2 it is found the E plane, which is in res/4
xe = theta(round(res/4),:)*180/pi;
ye=20*log10(abs(Etot(round(res/4),:)).*cos(theta(round(res/4),:)));
ye = ye-max(ye);
plotPlane(xe, ye);

% In phi = 0 it is found the E plane, which is in 1
xh = theta(1,:)*180/pi;
yh=20*log10(abs(Etot(1,:)).*cos(theta(1,:)));
yh = yh-max(yh);
plotPlane(xh, yh);

% Second, find the bw
findBw(xe, ye);
findBw(xh, yh);