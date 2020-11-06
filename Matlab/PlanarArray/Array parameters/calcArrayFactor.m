function [F] = calcArrayFactor(A, M, N, phix, phiy, u)
    % Calculates the array facto with a given feeding distribution
    F = zeros(size(u));

    for ii = 0:M-1        % For all the elements in the x direction
        for jj = 0:N-1    % For all the elements in the y direction
            F = F + A(ii+1,jj+1)*exp(1i*(ii*phix + jj*phiy));
        end
    end
end

