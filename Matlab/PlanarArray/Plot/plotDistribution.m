function [] = plotDistribution(A, M, N, path, fileName)
figure('Color',[1 1 1]);
x = 1:M;
y = 1:N;
[x, y] = meshgrid(x, y);
size(x)
size(y)
pcolor(x, y, A);
colormap winter
colorbar;
saveas(gca, [path, fileName],'epsc');
hold off;
end

