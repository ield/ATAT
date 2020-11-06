function [] = plotPlane(x, y1, y2, min)
%     figure('Color', [1 1 1]);
    x = [-fliplr(x) x];
    y1 = [fliplr(y1) y1];
    plot(x, y1, '-');
    hold on;
    
    y2 = [fliplr(y2) y2];
    plot(x, y2, '-');
    
    
    xlim([-1 1]);
    ylim([min 0]);
    xlabel('u, v');
    ylabel('dB');
    legend('E plane', 'H plane');
end

