function [] = plotReconstruction(x, y, E_x, E_y, A, path, distribution)
    minZ = 0;
    set(0, 'DefaultAxesFontName', 'Times New Roman');
    % Reconstructed field (amplitude)
    maxim_normalize = max(max(sqrt(abs(E_x.^2 + E_y.^2))));
    E_x_treated = abs(E_x)/maxim_normalize;
    E_y_treated = abs(E_y)/maxim_normalize;
    maxZ = 1;
    
    figure('Color',[1 1 1]);
    set(gcf,'position',[100,100,900,600]);
    % First plots notmal and then with zoom
    subplot(2, 2, 1);
    plot3Dxy(x, y, E_x_treated, minZ, maxZ, 'Normalized radiated field E_x');
    
    subplot(2, 2, 2);
    plot3Dxy(x, y, E_y_treated, minZ, maxZ, 'Normalized radiated field E_y');
    
    subplot(2, 2, 3);
    plot3Dxy(x, y, 20*log10(E_x_treated), -100, 0, 'Normalized radiated field E_x (dB) (zoom)');
    xlim([-0.01 0.05]);
    ylim([-0.01 0.05]);
    
    subplot(2, 2, 4);
    plot3Dxy(x, y, E_y_treated, minZ, maxZ, 'Normalized radiated field E_y (zoom)');
    xlim([-0.01 0.05]);
    ylim([-0.01 0.05]);
    
    saveas(gca, [path, distribution, '_reconstructed'],'epsc');
    
    % Reconstructed field (phase)
    phase_x = angle(E_x)*180/pi + 180;
    phase_y = angle(E_y)*180/pi + 180;
    maxZ = 360;
    
    % First plots notmal and then with zoom
    figure('Color',[1 1 1]);
    set(gcf,'position',[100,100,900,600]);
    
    subplot(2, 2, 1);
    plot3Dxy(x, y, phase_x, minZ, maxZ, 'Phase E_x');
    
    subplot(2, 2, 2);
    plot3Dxy(x, y, phase_y, minZ, maxZ, 'Phase E_y');
    
    subplot(2, 2, 3);
    plot3Dxy(x, y, phase_x, minZ, maxZ, 'Phase E_x (zoom)');
    xlim([-0.01 0.05]);
    ylim([-0.01 0.05]);
    
    subplot(2, 2, 4);
    plot3Dxy(x, y, phase_y, minZ, maxZ, 'Phase E_y (zoom)');
    xlim([-0.01 0.05]);
    ylim([-0.01 0.05]);
    
    saveas(gca, [path, distribution, '_phase'],'epsc');
    
    % Original field
    A_mag = abs(A)/max(max(abs(A)));
    A_phase = angle(A)*180/pi+180;
    
    figure('Color',[1 1 1]);
    set(gcf,'position',[100,100,900,300]);
    
    [M, N] = size(A);
    x_ant = 1:M;
    y_ant = 1:N;
    [x_ant, y_ant] = meshgrid(x_ant, y_ant);
    
    subplot(1, 2, 1);    
    pcolor(x_ant, y_ant, A_mag);
    title('Original Field');
    colorbar;
    
    subplot(1, 2, 2);
    pcolor(x_ant, y_ant, A_phase);
    title('Original Phase');
    colorbar;
    
    saveas(gca, [path, distribution, '_original'],'epsc');   

end

