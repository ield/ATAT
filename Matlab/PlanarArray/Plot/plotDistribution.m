function [] = plotDistribution(A, M, N)
figure('Color',[1 1 1]);
x = 1:M;
y = 1:N;
[x, y] = meshgrid(x, y);
size(x)
size(y)
pcolor(x, y, A);
colorbar;
end

