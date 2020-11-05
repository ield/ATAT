function [Etotcar] = printAndPlotSteering(Etot, phi, theta, titleIn, path, fileName, maxNorm)
    %% Print title
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf('\n');
    fprintf(titleIn);
    fprintf('\n');
    figure('Color',[1 1 1]);
    
    hold on;
    %% Plot the radiated field
    %     [~, ~, Etotcar] = sph2cart(phi, pi/2-theta, abs(Etot));
    u=sin(theta).*cos(phi);
    v=sin(theta).*sin(phi);
%     u = linspace(-1, 1, res);
%     v = linspace(-1, 1, res);
    Etotcar=Etot.*cos(theta);
    plot3Duv(u, v, 20*log10(abs(Etotcar))-maxNorm, -50, 'Diagram Normalized (dB)');

    %% Save the file
    saveas(gca, [path, fileName],'epsc');
    hold off;

end

