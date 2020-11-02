function [] = printAndPlotArrayParameters(u, v, Etot, phi, theta, res, title)
    %% Print title
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n');
    fprintf(title);
    fprintf('\n');

    %% Plot the radiated field
    [~, ~, Etotcar] = sph2cart(phi, pi/2-theta, abs(Etot));
    plot3Duv(u, v, abs(Etotcar), 0, 'Radiation element (lineal)');
    plot3Duv(u, v, 20*log10(abs(Etotcar)), 0, 'Array Factor (dB)');
    plot3Duv(u, v, 20*log10(abs(Etotcar/max(max(Etotcar)))), -30, 'Array Factor Normalized (dB)');

    %% Parameters BW, D0, SLL
    % First it is calculated the normalized radiation diagram in planes E and H
    % In phi = pi/2 it is found the E plane, which is in res/4
    xe = theta(round(res/4),:)*180/pi;
    ye=20*log10(abs(Etot(round(res/4),:).*cos(theta(round(res/4),:))));
    ye = ye-max(ye);
    plotPlane(xe, ye, 'Normalized radiation diagram E plane (dB)');

    % In phi = 0 it is found the E plane, which is in 1
    xh = theta(1,:)*180/pi;
    yh=20*log10(abs(Etot(1,:).*cos(theta(1,:))));
    yh = yh-max(yh);
    plotPlane(xh, yh, 'Normalized radiation diagram H plane (dB)');

    % Step 2. Find the bw
    bwe = findBw(xe, ye, 'E');
    bwh = findBw(xh, yh, 'H');

    % Step 3. Find the sll
    calcSLL(ye, 'E');
    calcSLL(yh, 'H');

    % Step 4. Calculate directivity approximation for plannar arrays balanis
    % pag 51
    d0 = 32400/(bwe*bwh);
    D0 = 10*log10(abs(d0));
    fprintf('The directivity is %f dB\n', D0);

end

