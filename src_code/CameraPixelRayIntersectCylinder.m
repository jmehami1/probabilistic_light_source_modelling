function [xInt, valid] = CameraPixelRayIntersectCylinder(normImgPt, radius, xStart, xEnd)
% Determine the 3D location of the intersection of a pixel ray with a
% cylinder defined by an axis and radius. The cylinder exists between a set
% of starting and ending points. The pixels is given as a normalised
% image coordinate, which removes the effects of the intrinsic parameters
% and places the points on a virtual plane that is one unit infront of the
% camera (projective space which has no units).
% All points and vectors should be described in the coordinate
% frame of the camera.
% INPUTS:
%       normImgPt - normalised image coordinate (x,y,1) (1x3)
%       radius - cylinder radius
%       xStart - 3D point on cylinder axis that represents the cylinder's
%               starting point. (1x3)
%       xEnd - 3D point on cylinder axis that represents the cylinder's end
%               point (1x3)
% OUTPUTS:
%       xInt - intersected 3D point between ray and cylinder (1x3) w.r.t to
%       camera's coordinate frame
%       valid - valid intersection found
%Author: Jasprabhjit Mehami, 13446277

%direction vector of ray starts at camera origin towards normalised image
%coordinate
rayDir = normImgPt;


%line cylinder intersection. The camera origin is treated as the starting point of the
%ray
[xIntSol, valid] = LineCylinderIntersection([0, 0, 0, rayDir], [xStart, xEnd, radius]);

if ~valid
    xInt = 0;
    return;
end

%first solution is the closest one to the origin
xInt = xIntSol(1,:);


% cylinder axis
axis = [xStart (xEnd - xStart)];

% compute position on axis
ts = linePosition3d(xInt, axis);

% check bounds
if ~(ts>=0 && ts<=1)
    valid = false;
end

end
