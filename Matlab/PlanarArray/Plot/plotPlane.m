function [] = plotPlane(x, y1, y2, min, alphax, alphay)
%     figure('Color', [1 1 1]);
%     x = [-fliplr(x) x];
%     y1 = [fliplr(y1) y1];
    plot(x, y1, '-');
    hold on;
    
%     y2 = [fliplr(y2) y2];
    plot(x, y2, '-');
    
    
    xlim([-1 1]);
    ylim([min 0]);
    xla = ['u = ' num2str(alphax) 'º, v = ' num2str(alphay) 'º'];
    xlabel(xla);
    ylabel('dB');
    legend('v (H plane)', 'u (E plane)');
end

