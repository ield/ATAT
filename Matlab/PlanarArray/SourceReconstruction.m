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

    x = ((-0.5/delta_u):delta_x:(0.5/delta_u-delta_x))*lambda; size(x)
    y = ((-0.5/delta_v):delta_y:(0.5/delta_v-delta_y))*lambda; size(y)
%     x = linspace(-0.5*lambda/delta_u, 0.5*lambda/delta_u, res);
%     y = linspace(-0.5*lambda/delta_v, 0.5*lambda/delta_v, res);
    [x, y] = meshgrid(x, y);

    % Define theta and phi fields
    theta = asin(sqrt(u.^2 + v.^2));    
    out_unit_circle = find(imag(theta)~=0);    
    theta(out_unit_circle) = pi/2;       % Limit to the unit circle: u^2 + v^2 > 1
    phi = atan2(u, v);

    Ecp = Etot(1:(end-1), 1:(end-1));
    
    E_theta = Ecp.*sin(phi);
    E_phi = Ecp.*cos(phi);

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