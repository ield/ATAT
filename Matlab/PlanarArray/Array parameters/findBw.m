function [bw] = findBw(x, y, plane, ref, res)
%%%%%%%%%%%%%%%%DOES NOT WORK%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds the beam width of a normalized radiation pattern. 
% ref is u0 for v and v0 for u
y = 10.^(y/20);
maximY = max(y);
y = y/maximY;
indexMax = find(y == 1);

% We separate the array into two symmetric parts
y = circshift(y, length(y)/2-indexMax); % The maximum is moved to the center

y1 = y(1:length(y)/2); % Array separation
y1 = flip(y1);
x1 = x(1:length(x)/2); 
x1= flip(x1);

y2 = y(length(y)/2+1:end); % Array separation
x2 = x(length(x)/2+1:end);

% subplot(1, 2, 1)
% plot(y1)
% subplot(1, 2, 2)
% plot(y2)

% It is calculated the first time the value goes 3db under the maximum
% For y1, it is necessary to turn it around. Since we are in linear
dB3 = 10^(-3/20);

[~, indexBw1] = min(abs(y1-dB3));
% indexBw1 = indexBw1 - 1;   % -1 because the first index is at 0
[~, indexBw2] = min(abs(y2-dB3));

x01 = x1(indexBw1);
x02 = x2(indexBw2);

theta1 = asin(sqrt(x01^2)); %Conversion from u, v to theta
theta2 = asin(sqrt(x02^2));

% theta1*180/pi
% theta2*180/pi

bw = theta1+theta2;
bw = bw*180/pi;    
% fprintf('ref = %f. x01 = %f. y01 = %fThe  %s = %f º\n', plane, bw);
fprintf('The beamwidth in %s = %f º\n', plane, bw);
end

