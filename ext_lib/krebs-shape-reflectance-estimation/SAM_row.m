function ang = SAM_row(x,y)
% Spectral angle mapper (SAM) scale invariant spectral similarity measure.
% Returns the angle between the two N-dimensional row vectors of spectral
% radiances.
%
% INPUTS:
%       x - spectral radiance vector [1 x N]
%       y - spectral radiance vector [1 x N]
% OUTPUT:
%       ang - Angle between radiance vectors. This is a value between 0 and
%           pi/2, where 0 means the vectors are exactly the same.
%
% Author: Jasprabhjit Mehami, 13446277

ang = acos(sum(x.*y, 'all')/(norm(x)*norm(y)));
end

