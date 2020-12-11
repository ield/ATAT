function [] = plotNoisy(Etot, u, v, res, titleIn, path, fileName)
    x = linspace(-1, 1, res);
    maxNorm = max(max(20*log10(abs(Etot))));
    
    yu=20*log10(abs(Etot(:, round(res*(1/2)))));
    yu = yu-maxNorm;
    
    yv=20*log10(abs(Etot(round(res*(1/2)), :)));
    yv = yv-maxNorm;
    yv = yv';
    
    plotPlane(x, yu, yv, -50, 0, 0);
    hold off
    
    title(titleIn);
end

