function [x, y, z] = sphZ2cart(r, azi, elev)
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
% 
% Author: Jasprabhjit Mehami, 13446277

x = r .* cos(elev) .* sin(azi);
z = r .* cos(elev) .* cos(azi);
y = r .* sin(elev);

% % y = r .* cos(elev) .* sin(azi);
% % z = r .* sin(elev);
% 
% elev = pi/2 - elev;
% 
% % [x,y,z] = sph2cart(azi,elev,r);
% 
% % x = r .* cos(elev) .* cos(azi);
% % y = r .* cos(elev) .* sin(azi);
% % z = r .* sin(elev);
% 
% x = r .* cos(azi) .* cos(elev);
% y = r .* cos(azi) .* sin(elev);
% z = r .* sin(azi);

end

