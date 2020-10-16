clear;
f = 15e9;
lambda = 3e8 / f;
k = 2*pi/lambda;

%% Calculating D
% Knowing that bw = 70*lambda/d, where d is the diameter, and that d =
% 4*pi/angle*angle, it is calculated the diameter

direc = 35;                                 % Directivity (dBi)
% D = 10^(direc/20)*7*lambda*sqrt(pi)/36     % Diameter (m)
D = 0.42;
%% Calculating f/d
% Using theta0 = 64º,
theta_0 = 64*pi/180;
f_d = 1/(4*tan(theta_0/2))

%% Calculating 1%F and 3%F
foc = f_d*D
f_01 = 0.01*foc;
f_03 = 0.03*foc;

z_01 = foc + f_01
z_03 = foc + f_03

y_01 = f_01

%% Calculating the aperture efficiency
% Values obtained in grasp
direc = [35.24 35.095 34.92];
sp = [0.4043 0.2044 0.1297];

% Knowing d = 4pi*A/lambda^2*eff
A = pi*(D/2)^2;
eff = lambda^2/4/pi/A*10.^(direc/10)

% Knowing that the etot = esp*eap
esp = 10.^(-sp/10);
eap = eff./esp


