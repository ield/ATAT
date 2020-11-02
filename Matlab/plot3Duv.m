function [] = plot3Duv(u, v, Z, minZ, zLabelText)
% Function to plot in 3D the values received as parameters in u, v
figure('Color',[1 1 1]);
surf(u, v, Z);
xlabel('u');
ylabel('v');    
zlabel(zLabelText);
zlim([minZ max(max(Z))]);
end

