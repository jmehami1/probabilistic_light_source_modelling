function reflect = ReflectionVector(normal, incident)
%Find the reflected vector at some point given a normal at that point and
%the incident vector. The normal and incident vector need to be in the same
%coordinate frame. The normal must be normalised.

% INPUTS:
%       normal - normal direction vector of surface at some point
%       incident - incidental vector pointing into surface
% OUTPUT:
%       reflect - reflected vector going out of surface
% Author: Jasprabhjit Mehami, 13446277


reflect = incident - 2.*dot(normal, incident).*normal;

end

