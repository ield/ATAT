function [] = plotPlane(x1, y1, ylab)
    figure('Color', [1 1 1]);
    x1 = [-fliplr(x1) x1];
    y1 = [fliplr(y1) y1];
    plot(x1, y1)
    xlim([-90 90]);
    ylim([max(y1)-50 max(y1)]);
    xlabel('\phi');
    ylabel(ylab);
end

