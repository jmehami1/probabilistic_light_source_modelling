function [outputArg1,outputArg2] = SolveEstimatedReflectancePatch(I, s, g, mask)
% Solve the estimated reflectance for a single patch through a constrainted linear
% least squares minimisation of the dichromatic model.
% INPUTS:
%       I - normalised radiance measurements [v x u x bands]
%       s - radiant intensity magnitude from light source to the 3D
%           locations of I [v x u x bands]
%       n - surface normal direction vectors at the 3D locations [v x u]
%       dirL - direction light vector from the 3D locations to the source [v x u]

%       dirC - direction viewing vector of the camera from the 3D point to its optical
%           centre [pixels x 3]
%       patchSize - number of band-pixels grouped such that they are assumed to
%           have the same reflectances
%       distS - distance of points from light source
%       pixProjArea - projected area of each pixel on the surface.
%       pix - pixel locations in hyperspectral image
%       maxPixDist - max distance allowed between pixels to be considered a
%       patch

maskVec = mask(:)

Ivec = I(:);

%create a mask of all non-zero pixels
pixMask = Ivec > 0;



[numV, numU, numBands] =  size(I);


B = permute(I,[2 1 3]);
out = B(:)



end

