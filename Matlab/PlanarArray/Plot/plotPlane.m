function [] = plotPlane(x, y1, y2)
%     figure('Color', [1 1 1]);
    x = [-fliplr(x) x];
    y1 = [fliplr(y1) y1];
    plot(x, y1, '-');
    hold on;
    
    y2 = [fliplr(y2) y2];
    plot(x, y2, '-');
    
    
    xlim([-90 90]);
    ylim([-50 0]);
    xlabel('\phi');
    ylabel('dB');
    legend('E plane', 'H plane');
end

