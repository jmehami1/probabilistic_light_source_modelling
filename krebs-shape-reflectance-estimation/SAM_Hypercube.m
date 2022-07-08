function T = SAM_Hypercube(X,Y)
% Sum of Spectral angle mapper (SAM) scale invariant spectral similarity measure.
% Returns the total angle between each hypercube pixel
%
% INPUTS:
%       X - spectral hypercube 
%       Y - spectral hypercube 
% OUTPUT:
%       ang - Angle between radiance vectors. This is a value between 0 and
%           pi/2, where 0 means the vectors are exactly the same.
%
% Author: Jasprabhjit Mehami, 13446277

assert(isequal(size(X), size(Y)), "hypercubes do not have the same size");
assert(isequal(length(size(X)), length(size(Y)), 3), "hypercubes do not have 3 dimensions");

angVec = real(acos(sum(X.*Y, 3)./(vecnorm(X,2,3).*vecnorm(Y,2,3))));
angTotal = sum(angVec, 'all', 'omitnan');
angMean = mean(angVec, 'all', 'omitnan');
angSTD = std(angVec, 0, 'all', 'omitnan');

% normX = normalize(X, 3, 'norm');
% normY = normalize(Y, 3, 'norm');normX = X./sum(X,'all');
normX = X./sum(X,[1,2]);

normY = Y./sum(Y,[1,2]);
angSSE = sumsqr(normX - normY);


T = table(angTotal, angMean, angSTD, angSSE, 'VariableNames', ["Total SAM", "Mean of SAM", "STD of SAM", "SSE of Normalised S"]);
end

