function symPt = SymmetricPoint_LightSrc_Linescan(pt, lPos, lDir, vpNorm)
%Given a point on a view-plane of a line-scan camera, find the symmetrical 
%point about the principal direction vector of a light source which is also
%on the view-plane. Everything must be in the coordinate frame of the
%line-scan camera.
% INPUTS:
%       pt - point of interest
%       lPos - position of the light source
%       lDir - principal direction vector of light source
%       vpNorm - normal of line-scan view-plane
% OUTPUT:
%       symPt - symmetrical point of pt about lDir on the view-plane
% Author: Jasprabhjit Mehami, 13446277

%project light source location vector to view-plane
proj_lPos = PlaneProjection(lPos, vpNorm);

%project light source direction vector to view-plane
proj_lDir = PlaneProjection(lDir, vpNorm);
proj_lDir = proj_lDir./norm(proj_lDir); %normalise to unit vector

%projected incident vector on view-plane
proj_incid = proj_lPos - pt;

%find the projected reflection vector on the view-plane
proj_refl = ReflectionVector(proj_lDir, proj_incid);

%calculate the symmetrical point on view-plane relative to the line-scan
%camera origin
symPt = proj_lPos + proj_refl;

end

