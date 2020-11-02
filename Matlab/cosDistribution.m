function [A] = cosDistribution(M, N)
% Script to simulate a cos distribution
H = 0.5;
% Feeding in x
Ax = zeros(1, M);

for ii = 0:M-1
    elem = ii-(M-1)/2;
    Ax(ii+1) = 1+H*(cos(pi*elem/(M-1)))^2;
end

% Feeding in y
Ay = zeros(1, N);

for ii = 0:N-1
   elem = ii-(N-1)/2;
    Ay(ii+1) = 1+H*(cos(pi*elem/(N-1)))^2;
end

% Combine feeding
A = zeros(M, N);
for ii = 1:M-1
    for jj = 1:N-1
        A(ii, jj) = Ax(ii)*Ay(jj);
    end
end
end

