function [bw] = findBw(x, y)
% Finds the beam width of a normalized radiation pattern
indexBw = find(y<-3, 1);
bw = 2*x(indexBw);
fprintf('The beamwidth in E = %f º\n', bw);
end

