function [] = plot3Dxy(x, y, Z, minZ, maxZ, titleText)
% Function to plot in 3D the values received as parameters in u, v
surf(x, y, Z);
view(2)
grid off;
colormap winter
shading interp
caxis([minZ maxZ]);
colorbar;
xlabel('x');
ylabel('y');    
title(titleText);
zlim([minZ max(max(Z))]);
end

