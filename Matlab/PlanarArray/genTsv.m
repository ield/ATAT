load 'PlanarArray.mat';             % Here it is all the information of the planar array (positions, etc)
load 'cosineDistribution.mat';      % Here it is all the information of the distribution

M = length(A(1, :));
N = length(A(:, 1));

phase = zeros(M, N);

for ii = 1:M
    for jj = 1:N
        phase(ii, jj) = (ii*2*pi*dx/lambda*coorX(ii, jj) + jj*2*pi*dy/lambda*coorY(ii, jj))*180/pi;
    end
end
% Now we have all the phase distributions so we write the TSV file

fileID = fopen('arrayPar.tsv','w');
% fprintf(fileID, 'No.\tX\tY\tZ\tAmplitude\tPhase\n');    % First row in Matlab
no = 1;
for ii = 1:M
    for jj = 1:N
%         no_str = num2str(no);
        x = num2str(coorX(ii, jj)*1e3);
        y = num2str(coorY(ii, jj)*1e3);
        z = '0';
        amp = num2str(A(ii, jj));
        ph = num2str(phase(ii, jj));
        
%         newLine = [no_str '\t' x '\t' y '\t' z '\t' amp '\t' ph '\r'];
        newLine = [x '\t' y '\t' z '\t' amp '\t' ph '\r'];
        fprintf(fileID, newLine);
        
        no = no+1;
    end
end

fclose(fileID);

