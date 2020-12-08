load 'PlanarArray.mat';             % Here it is all the information of the planar array (positions, etc)
load 'failingDistribution.mat';      % Here it is all the information of the distribution

M = length(A(1, :));
N = length(A(:, 1));

phase = zeros(M, N);

for ii = 1:M
    for jj = 1:N
        phase(ii, jj) = ((ii)*2*pi*dx/lambda*(coorX(ii, jj)-Lx) + (jj)*2*pi*dy/lambda*(coorY(ii, jj)-Ly))*180/pi;
    end
end
% Now we have all the phase distributions so we write the TSV file

fileID = fopen('arrayPar.tsv','w');
% The file is defined in file:///C:/Program%20Files%20(x86)/CST%20Studio%20Suite%202020/Online%20Help/mergedProjects/DES/userinterface/dialog/array_definition_tsv_file.htm
fprintf(fileID, '# Created by ield\t\t\t\t\t\t\t\t\n');    % First row in Matlab
fprintf(fileID, '# On 11,11,2020\t\t\t\t\t\t\t\t\n');    
fprintf(fileID, '# unit: meters\t\t\t\t\t\t\t\t\n');
fprintf(fileID, '# design frequency: 60000000 Hz\t\t\t\t\t\t\t\t\n');
fprintf(fileID, '# Element\tX\tY\tZ\tMagnitude\tPhase\tPhi\tTheta\tGamma\n');    
for ii = 1:M
    for jj = 1:N
        
        elemStr = [num2str(ii) ',' num2str(jj)];
        x = num2str(coorX(ii, jj)-Lx/2);
        y = num2str(coorY(ii, jj)-Ly/2);
        z = '0';
        amp = num2str(A(ii, jj));
        ph = num2str(phase(ii, jj));
        
        newLine = [elemStr '\t' x '\t' y '\t' z '\t' amp '\t' ph '\t0\t0\t0\n'];
        
        fprintf(fileID, newLine);
        
    end
end

fclose(fileID);

