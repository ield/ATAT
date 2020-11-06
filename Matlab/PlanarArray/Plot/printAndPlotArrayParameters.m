function [] = printAndPlotArrayParameters(Etot, u, v, res, titleIn, path, fileName, maxNorm)
    %% Print title
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf('\n');
    fprintf(titleIn);
    fprintf('\n');
    figure('Color',[1 1 1]);
    set(gcf,'position',[100,100,1500,300]);
    
    min = -60;
    
    hold on;
    %% Plot the radiated field
    
    if(nargin < 8)  %In the case it is uniform distribution
        maxNorm = max(max(20*log10(abs(Etot)))); %maximum used to normalize
    end
    
    subplot(1, 3, 1:2);
    plot3Duv(u, v, 20*log10(abs(Etot))-maxNorm, min, 'Diagram Normalized (dB)');
    


    %% Parameters BW, D0, SLL
    % First it is calculated the normalized radiation diagram in planes E and H
    subplot(1, 3, 3);
    x = linspace(-1, 1, res);
    
    % In phi = pi/2 it is found the E plane, which is in res/4 In the case
    % there is beam sgteering this changes, so in this case the angle is
    % added
    ye=20*log10(abs(Etot(round(res/2),:)));
    ye = ye-maxNorm;
    
    % In phi = 0 it is found the E plane, which is in 1
    yh=20*log10(abs(Etot(:,round(res/2))));
    yh = yh-maxNorm;
    yh = yh';
    
    plotPlane(x, ye, yh, min);
    hold off
    
    title(titleIn);

    % Step 2. Find the bw
    bwe = findBw(x, ye, 'E');
    bwh = findBw(x, yh, 'H');

    % Step 3. Find the sll
    calcSLL(ye, 'E');
    calcSLL(yh, 'H');

    % Step 4. Calculate directivity approximation for plannar arrays balanis
    % pag 51
    d0 = 32400/(bwe*bwh);
%     d0 = 4*pi*max(max(abs(u, v,).^2))/(sum(sum(abs(u, v,).^2.*sin(theta)))/numel(u, v,));
%     d0 = 4*pi*max(max(abs(u, v,Etot).^2))/(sum(sum(abs(u, v,Etot).^2))/numel(u, v,Etot));
   
    D0 = 10*log10(d0);

    fprintf('The directivity is %f dB\n', D0);
    
    %% Save the file
    saveas(gca, [path, fileName],'epsc');
    hold off;

end

