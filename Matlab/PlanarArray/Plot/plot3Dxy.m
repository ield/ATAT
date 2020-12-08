function [] = plot3Dxy(x, y, Z, minZ, zLabelText)
% Function to plot in 3D the values received as parameters in u, v
surf(x, y, Z);
view(2)
grid off;
colormap winter
shading interp
caxis([minZ 1]);
colorbar;
xlabel('x');
ylabel('y');    
zlabel(zLabelText);
zlim([minZ max(max(Z))]);
end

