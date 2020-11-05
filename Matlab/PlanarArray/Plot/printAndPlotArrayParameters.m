function [Etotcar] = printAndPlotArrayParameters(Etot, phi, theta, res, titleIn, path, fileName, maxNorm)
    %% Print title
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    fprintf('\n');
    fprintf(titleIn);
    fprintf('\n');
    figure('Color',[1 1 1]);
    set(gcf,'position',[100,100,1500,300]);
    
    hold on;
    %% Plot the radiated field
    %     [~, ~, Etotcar] = sph2cart(phi, pi/2-theta, abs(Etot));
    u=sin(theta).*cos(phi);
    v=sin(theta).*sin(phi);
%     u = linspace(-1, 1, res);
%     v = linspace(-1, 1, res);
    Etotcar=Etot.*cos(theta);
    
    if(nargin < 8)  %In the case it is uniform distribution
        maxNorm = max(max(20*log10(abs(Etotcar)))); %maximum used to normalize
    end
    
    subplot(1, 4, 1:3);
    plot3Duv(u, v, 20*log10(abs(Etotcar))-maxNorm, -50, 'Diagram Normalized (dB)');
    


    %% Parameters BW, D0, SLL
    % First it is calculated the normalized radiation diagram in planes E and H
    subplot(1, 4, 4);
    x = theta(1,:)*180/pi;
    
    % In phi = pi/2 it is found the E plane, which is in res/4 In the case
    % there is beam sgteering this changes, so in this case the angle is
    % added
    ye=20*log10(abs(Etot(round(res/4),:).*cos(theta(round(res/4),:))));
    ye = ye-maxNorm;
    
    % In phi = 0 it is found the E plane, which is in 1
    yh=20*log10(abs(Etot(1,:).*cos(theta(1,:))));
    yh = yh-maxNorm;
    
    plotPlane(x, ye, yh);
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
%     d0 = 4*pi*max(max(abs(Etot).^2))/(sum(sum(abs(Etot).^2.*sin(theta)))/numel(Etot));
%     d0 = 4*pi*max(max(abs(Etotcar).^2))/(sum(sum(abs(Etotcar).^2))/numel(Etotcar));
   
    D0 = 10*log10(d0);

    fprintf('The directivity is %f dB\n', D0);
    
    %% Save the file
    saveas(gca, [path, fileName],'epsc');
    hold off;

end

