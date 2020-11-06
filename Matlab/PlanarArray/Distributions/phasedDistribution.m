function [A] = phasedDistribution(A, M, N, alphax,alphay)

    for ii = 1:M
        for jj = 1:N
            A(ii, jj) = A(ii, jj)*exp(1i*(ii*alphax + jj*alphay));
        end
    end
    
end

