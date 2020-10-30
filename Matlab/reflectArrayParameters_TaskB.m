% Reflectarray elements
% Measures of the reflectarray in UPM paper
% All units in the app are in mm
clear;
c = 3e8;
%% Frequencies
fmax = 13e9;
fmin = 11e9;
Nfreq = 3;          % Number of frequencies evaluated
f0 = (fmax + fmin)/2;
DiffFreq = (fmax - f0)/floor(Nfreq/2);  % Distance between frequencies
fprintf('The central frequency is %f GHz and the DF is %f GHz.\n', f0/1e9, DiffFreq/1e9);

%% Feed Definition (m)
xf = -40e-3;
yf = 0;
zf = 195e-3;
QF = 13;

%% Cell Dimensions
% Periodic cell definition
lambda = c/f0;
esub = 1.05;    % Dielectric constant of the separator
lambdag = lambda/sqrt(esub);
Lmax = 0.8*lambdag;      % To avoid grating lobes
Px = 15e-3;
Py = 15e-3;
fprintf('The maximum length of the periodic structure is %f mm Px = %f mm, and Py = %f mm\n', Lmax*1e3, Px*1e3, Py*1e3);

%% Number of elements
diam = 18*3e-2;
Nx = floor(diam / Px);
Ny = floor(diam / Py);
fprintf('The number of elements is Nx = %i.   Ny = %i.\n', Nx, Ny);

%% Sandwich definition
% Unities of the separator: the patches are printed in the dielectic and
% the layers are separated by the separator. In the paper the separator is
% Qz and the substrate is Kapton. We are focusing on the separator.
Ts = 0.25e-3;   % Thickness of the separator
ersep = 3.2;    % Dielectric constant of the separator
tandsep = 0.004;
erImag = -ersep * tandsep;   % See Pozar Pag 10
alpha  = 0.85;  % Assumed, since the 3 layers are equally thick
fprintf('ER Imag = %f\n', erImag);

%% Substrate/Superstrate
% If there is no superstrate Ts = 0
Tsub = 3e-3;   % Thickness of the separator
% esub = 3.7;    % Dielectric constant of the separator
tandsub = 0.003;
erImag = -ersep * tandsub;   % See Pozar Pag 10
fprintf('ER sub Imag = %f\n', erImag);

%% Patch Size Range
% Random values, only taking care that the Anmax is not bigger than Px. In
% the fig 3 they are placed between 1 and 2.8

%% Angle of radiation
Theta0x = 20;
Theta0y = 30;

%% Expected gain
eff = 0.65;
G = 4*pi*pi*(diam/2)^2/lambda^2*eff;
fprintf('The expected gain is %f\n', 10*log10(G));