function [A] = binomialDistribution(M, N)
% Script to simulate a binomial distribution

% Feeding in x
Ax = zeros(1, M);

for ii = 0:M-1
    Ax(ii+1) = nchoosek(M-1,ii);
end

% Feeding in y
Ay = zeros(1, N);

for ii = 0:N-1
    Ay(ii+1) = nchoosek(N-1,ii);
end

% Combine feeding
A = zeros(M, N);
for ii = 1:M-1
    for jj = 1:N-1
        A(ii, jj) = Ax(ii)*Ay(jj);
    end
end
end

