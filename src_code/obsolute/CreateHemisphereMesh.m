function [xHem,yHem,zHem] = CreateHemisphereMesh(noPtsHem, distFromSrc, orien, locCentre)
% Create a mesh of points that lay on the surface of a hemisphere (where z
% is positive). Also rotate about the centre and move the center to a specified
% location of the light source.
% Inputs:
%       noPtsHem - number of points on the hemisphere
%       distFromSrc - distance from the centre of the source/ radius of
%       hemisphere
%       orien - Orientation of the hemisphere
%       locCentre - Location of the centre
% Outputs:
%       X,Y,Z coordinates of the mesh points on hemisphere
% Author: Jasprabhjit Mehami, 13446277

%create hemisphere mesh
[xHem,yHem,zHem] = sphere(noPtsHem);
[rows, ~] = size(zHem);
xHem = xHem((rows-1)/2:end,:);
yHem = yHem((rows-1)/2:end,:);
zHem = zHem((rows-1)/2:end,:);

%rotate hemisphere according to light source direction
[rows, cols] = size(zHem);
pts = zeros(3,cols);

for i = 1:rows
    pts(1,:) = xHem(i, :);
    pts(2,:) = yHem(i, :);
    pts(3,:) = zHem(i, :);

    ptsTrans = orien*pts;
    
    xHem(i, :) = ptsTrans(1,:);
    yHem(i, :) = ptsTrans(2,:);
    zHem(i, :) = ptsTrans(3,:);

end

%move hemisphere to light source location
xHem = distFromSrc.*xHem + locCentre(1);
yHem = distFromSrc.*yHem + locCentre(2);
zHem = distFromSrc.*zHem + locCentre(3);

end

