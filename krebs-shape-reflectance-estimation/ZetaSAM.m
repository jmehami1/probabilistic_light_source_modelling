function zeta = ZetaSAM(x, y, r)
% Negative exponential of the Spectral angle mapper (SAM) scale invariant 
% spectral similarity measure angle between the two N-dimensional vectors 
% of spectral radiance.
% 
% INPUTS:
%       x - spectral radiance vector [1 x N]
%       y - spectral radiance vector [1 x N]
%       r - bandwidth parameter which governs how close the vectors must be
%           until they are considered the same. Small values results in
%           stricter tolerances
% OUTPUT:
%       zeta - similiarity metric between 0 and 1, where 1 means the 
%           radiance vectors are exactly the same.
% 
% Author: Jasprabhjit Mehami, 13446277

samAngle = SAM_row(x,y);

zeta = exp(-(samAngle^2)/r);
end

