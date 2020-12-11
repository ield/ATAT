clear; close all;
path_image = '../../SourceReconstruction/Lab_08_Report/Images/';
path_antennas = 'AntennasSaved/';

%% Load the main info
load PlanarArray.mat
% lambda = lambda*1e3;
filename = '*Distribution.mat';
allAnt = dir([path_antennas, filename]);

%% Calculate the source for each distribution and plot
close all
for ii = 1:length(allAnt)
% for ii = 1:1
    % Load the antenna
    file = [path_antennas, allAnt(ii).name];
    load(file);
    distribution_name = allAnt(ii).name(1:end-16);
    
    % Define x and y axes according to the formula in slide 5
    res = 360;
    u = linspace(-1, 1, res);
    v = linspace(-1, 1, res);
    delta_u = u(2) - u(1);
    delta_v = v(2) - v(1);
    
    % In order to obtain better results it is convenient to delete the last
    % sample of the u and v components
    [u, v] = meshgrid(u(1:end-1),v(1:end-1));   
    res = res - 1;


    delta_x = 1/(delta_u*res);
    delta_y = 1/(delta_v*res);

%     x = ((-0.5/delta_u):delta_x:(0.5/delta_u-delta_x))*lambda; size(x)
%     y = ((-0.5/delta_v):delta_y:(0.5/delta_v-delta_y))*lambda; size(y)

    x = ((-0.5/delta_u):delta_x:(0.5/delta_u-delta_x)); size(x)
    y = ((-0.5/delta_v):delta_y:(0.5/delta_v-delta_y)); size(y)
    
    [x, y] = meshgrid(x, y);

    % Define theta and phi fields
    theta = asin(sqrt(u.^2 + v.^2));    
    out_unit_circle = find(imag(theta)~=0);    
    theta(out_unit_circle) = pi/2;       % Limit to the unit circle: u^2 + v^2 > 1
    phi = atan2(u, v);

    Ecp_noise = Etot(1:(end-1), 1:(end-1));
    
    E_theta = Ecp_noise.*sin(phi);
    E_phi = Ecp_noise.*cos(phi);

    % Define the plane wave spectrum
    P_theta = E_theta./cos(theta);
    P_theta(out_unit_circle) = 0;
    P_phi = E_phi./cos(theta);
    P_phi(out_unit_circle) = 0;

    P_x = P_theta.*cos(theta).*cos(phi) - P_phi.*sin(phi);
    P_y = P_theta.*cos(theta).*sin(phi) + P_phi.*cos(phi);    % Mistake in the slides: this is a +

    % Transform the plane wave into the desired field
    E_x = fftshift(fft2(ifftshift(P_x)));
    E_y = fftshift(fft2(ifftshift(P_y)));
    
    plotReconstruction(x, y, E_x, E_y, A, path_image, distribution_name);
%     close all;
end

%% Adding Noise
filename = 'uniformDistribution.mat';
file = [path_antennas, filename];
load(file);
distribution_name = 'noise';

Ecp = Etot(1:(end-1), 1:(end-1));
noise = 5*rand(size(Ecp)) + 5*1i*rand(size(Ecp));
Ecp_noise = Ecp + noise;

plotNoisy(Ecp_noise, u, v, res, 'Noisy Elements', path_image, distribution_name);

% Define x and y axes according to the formula in slide 5
res = 360;
u = linspace(-1, 1, res);
v = linspace(-1, 1, res);
delta_u = u(2) - u(1);
delta_v = v(2) - v(1);

% In order to obtain better results it is convenient to delete the last
% sample of the u and v components
[u, v] = meshgrid(u(1:end-1),v(1:end-1));   
res = res - 1;


delta_x = 1/(delta_u*res);
delta_y = 1/(delta_v*res);

%     x = ((-0.5/delta_u):delta_x:(0.5/delta_u-delta_x))*lambda; size(x)
%     y = ((-0.5/delta_v):delta_y:(0.5/delta_v-delta_y))*lambda; size(y)

x = ((-0.5/delta_u):delta_x:(0.5/delta_u-delta_x)); size(x)
y = ((-0.5/delta_v):delta_y:(0.5/delta_v-delta_y)); size(y)

[x, y] = meshgrid(x, y);

% Define theta and phi fields
theta = asin(sqrt(u.^2 + v.^2));    
out_unit_circle = find(imag(theta)~=0);    
theta(out_unit_circle) = pi/2;       % Limit to the unit circle: u^2 + v^2 > 1
phi = atan2(u, v);


E_theta = Ecp_noise.*sin(phi);
E_phi = Ecp_noise.*cos(phi);

% Define the plane wave spectrum
P_theta = E_theta./cos(theta);
P_theta(out_unit_circle) = 0;
P_phi = E_phi./cos(theta);
P_phi(out_unit_circle) = 0;

P_x = P_theta.*cos(theta).*cos(phi) - P_phi.*sin(phi);
P_y = P_theta.*cos(theta).*sin(phi) + P_phi.*cos(phi);    % Mistake in the slides: this is a +

% Transform the plane wave into the desired field
E_x = fftshift(fft2(ifftshift(P_x)));
E_y = fftshift(fft2(ifftshift(P_y)));

plotReconstruction(x, y, E_x, E_y, A, path_image, distribution_name);

% Delete the undesired field (it is only used knowledge of the position of
% the antenna. It is not used anything concerning the radiation.

% 1. Determine the positions in which the antenna is present
axis_min = -0.75;
axis_max = 8.75;

[~, antenna_position_index] = find(x(1, :)>=axis_min & x(1, :)<=axis_max);
antenna_position_index_min = antenna_position_index(1);
antenna_position_index_max = antenna_position_index(end);

% Set the field to 0 outside of that positions.
E_x_corrected = zeros(size(E_x));
E_x_corrected(antenna_position_index_min:antenna_position_index_max, ...
    antenna_position_index_min:antenna_position_index_max) = E_x(...
    antenna_position_index_min:antenna_position_index_max, ...
    antenna_position_index_min:antenna_position_index_max);

E_y_corrected = zeros(size(E_y));
E_y_corrected(antenna_position_index_min:antenna_position_index_max, ...
    antenna_position_index_min:antenna_position_index_max) = E_y(...
    antenna_position_index_min:antenna_position_index_max, ...
    antenna_position_index_min:antenna_position_index_max);

% Plot the reconstructed result
plotReconstruction(x, y, E_x_corrected, E_y_corrected, A, path_image, distribution_name);

% Calculate the P_x and P_y from the corrected field
P_x_corrected = fftshift(fft2(ifftshift(E_x_corrected)));
P_y_corrected = fftshift(fft2(ifftshift(E_y_corrected)));

% Transform to theta and phi
% The equations are worked in paper (atat notes page 9)
parenthesis = cos(theta).*sin(phi) + (cos(theta).*cos(phi)-cos(theta))./sin(phi);
P_theta_corrected = (P_y_corrected + P_x_corrected.*cos(phi)./sin(phi))./parenthesis;
P_phi_corrected = (P_theta_corrected.*cos(theta).*cos(phi) - P_x_corrected)./sin(phi);

% Transform the power wave into the field
E_theta_corrected = P_theta_corrected.*cos(theta);
E_phi_corrected = P_phi_corrected.*cos(theta);

E_cp_1 = E_theta_corrected./sin(phi);
E_cp_2 = E_phi_corrected./cos(phi);

E_cp_corrected = (E_cp_1 + E_cp_2)/2;

plotNoisy(E_cp_2, u, v, res, 'Noisy Elements', path_image, distribution_name);


