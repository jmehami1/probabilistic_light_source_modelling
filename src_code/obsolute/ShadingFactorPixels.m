function gVec = ShadingFactorPixels(n, dirL, distS, pixProjArea)
% Calculate the shading factor of pixel locations from a surface given in the
% lighting direction, surface normal, and projected pixel area on the
% surface. Inputs are all vectors, and output is a vector.
% INPUTS:
%       n - surface normal direction vectors at the 3D locations [pixels x 3]
%       dirL - direction light vector from the 3D pixel locations to the source [pixels x 3]
%       distS - distance of points from light source in metres [pixels x 1]
%       pixProjArea - projected area of each pixel on the surface in metres squared [pixels x 1]
% OUTPUT:
%       gVec - normalised shading factor [pixels x 1]

%dot product between light vector and surface normal
dirL_dot_n = dot(dirL, n, 2);
dirL_dot_n(dirL_dot_n < 0) = 0;

%solid angle of light hitting surface as captured by pixel
pixSolidAngle = (pixProjArea .* dirL_dot_n) ./ distS.^2;

%shading factor in absolute units
gVec = (dirL_dot_n.*pixSolidAngle)./(pi);

%normalise between max and min
gVec = mat2gray(gVec);

end

