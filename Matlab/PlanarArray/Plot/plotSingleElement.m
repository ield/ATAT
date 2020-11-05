function [] = plotSingleElement(u, v, Z, minZ, zLabelText)
% Function to plot in 3D the values received as parameters in u, v
surf(u, v, Z);
colormap summer
shading interp
caxis([0 1]);
colorbar;
xlabel('u');
ylabel('v');    
zlabel(zLabelText);
zlim([minZ max(max(Z))]);
end
