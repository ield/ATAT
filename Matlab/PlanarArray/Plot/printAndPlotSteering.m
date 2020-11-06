function [Etotcar] = printAndPlotSteering(Etot, u, v, titleIn, path, fileName, maxNorm)
    %% Print title
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf('\n');
    fprintf(titleIn);
    fprintf('\n');
    figure('Color',[1 1 1]);
    
    hold on;
    %% Plot the radiated field
    plot3Duv(u, v, 20*log10(abs(Etot))-maxNorm, -50, 'Diagram Normalized (dB)');

    %% Save the file
    saveas(gca, [path, fileName],'epsc');
    hold off;

end

