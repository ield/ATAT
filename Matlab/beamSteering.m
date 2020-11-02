function [alphax, alphay] = beamSteering(theta0, phi0, lambda, dx, dy)
% Function to solve the beam steering from the notes in diapo 22, solving
% the sytem of equations in ATAT-6 (notes)
theta0 = theta0*pi/180; %Convert to radians
phi0 = phi0*pi/180;
alphay = ((2*pi*dx*dy*sin(theta0))/lambda)/sqrt(dy^2+(dx*dx/(tan(phi0)*dy))^2);
alphax = alphay*dx/(tan(phi0)*dy);
end

