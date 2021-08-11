function [r, azi, elev] = cart2sphZ(x, y, z)
% Converts a 3D cartesian coordinate to a spherical coordinate defined as
% (r, azimuth, elevation). r is the euclidean distance, azimuth is the angle of
% the projected vector on the XZ-plane from the positive Z-axis,
% and elevation is the angle from the XZ-plane.
% INPUTS:
%       x - component of cartesian coordinate (1xn)
%       y - component of cartesian coordinate (1xn)
%       z - component of cartesian coordinate (1xn)
% OUTPUTS:
%       r - radius of point in spherical coordinates (1xn)
%       azi - azimuth angle of point in spherical coordinates (1xn)
%       elev - elevation angle of point in spherical coordinates (1xn)

% Author: Jasprabhjit Mehami, 13446277

azi = pi/2 - atan2(z,x);
elev = atan2(y,sqrt(x.^2 + z.^2));
r = sqrt(x.^2 + y.^2 + z.^2);

end

