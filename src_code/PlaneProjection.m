function projVec = PlaneProjection(vec,normPlane)
% Given a vector and a normal vector for any plane in the same linear
% space, project vector onto the plane defined by that normal vector.

projVec = vec - dot(vec, normPlane).* normPlane;

end
