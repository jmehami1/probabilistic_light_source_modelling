function rho = SolveEstimatedReflectance(I, s, n, dirL, dirC, distS, pixProjArea, pix, patchSize, optOptions, maxPixDist)
% Solve the for the estimated reflectance through the constrainted linear
% least squares minimisation of the dichromatic model.
% INPUTS:
%       I - normalised radiance measurements from an image line [bands x pixels]
%       s - radiant intensity magnitude from light source to the 3D
%           locations of I [bands x pixels]
%       n - surface normal direction vectors at the 3D locations [pixels x 3]
%       dirL - direction light vector from the 3D locations to the source [pixels x 3]
%       dirC - direction viewing vector of the camera from the 3D point to its optical
%           centre [pixels x 3]
%       patchSize - number of band-pixels grouped such that they are assumed to
%           have the same reflectances
%       distS - distance of points from light source
%       pixProjArea - projected area of each pixel on the surface.
%       pix - pixel locations in hyperspectral image
%       maxPixDist - max distance allowed between pixels to be considered a
%       patch

[numBands, numPix] = size(I);

%check patchSize
if (patchSize < 2)
    patchSize = 2;
elseif (patchSize > numPix)
    patchSize = numPix;
end

% pixel distance between corresponding points
pixDiffAbs = abs(diff(pix));
canBePatch = pixDiffAbs <= maxPixDist;

patchID = zeros(numPix, 1);
id = 1;
patchID(1) = id;
patchNumPix = zeros(numPix, 1);
curPatCount = 1;

%assigning a patch ID to each pixel location
for i = 2:numPix
    if curPatCount < patchSize
        if canBePatch(i - 1)
            patchID(i) = id;
        else
            patchNumPix(id) = curPatCount;
            id = id + 1;
            patchID(i) = id;
            curPatCount = 0;
        end
    else
        patchNumPix(id) = curPatCount;
        curPatCount = 0;
        id = id + 1;
        patchID(i) = id;
    end

    curPatCount = curPatCount + 1;
end

%number of patches is equal to the last patch ID
numPatch = id;
%Number of unknown reflectance values for each patch
numRho = numPatch*numBands; 



%shading factor for pixel
dirL_dot_n = dot(dirL, n, 2);
pixSolidAngle = (pixProjArea .* dirL_dot_n) ./ distS.^2;
g_coeff = (dirL_dot_n.*pixSolidAngle)./(pi);

% g_coeff = (pixProjArea .* dirL_dot_n.^2)./(pi.*distS.^2);

%specular reflected direction vector (if material was perfect mirror)
dirR = -dirL - 2.*dot(-dirL,n,2).*n;
dirR_dot_dirC = dot(dirR, dirC, 2);
%any negative dot products become zero
dirR_dot_dirC(dirR_dot_dirC < 0) = 0;

% c_coeff = zeros(numBands, numPix);
% for bandLoop = 1:numBands
%     %light vectors for current band
%     l = s(bandLoop, :)'.*dirL;
%     
%     %specular reflected vector (if material was perfect mirror)
%     r = -l + 2.*dot(l,n,2).*n;
%     
%     %radiant intensity along viewing direction. Any viewing direction that is
%     %not inline with the reflected vector direction, will be attenuated.
%     c_coeff(bandLoop,:) = dot(dirC, r, 2);
% end

%solve optimisation across patches

% for patchLoop = 1:numPatch
%     
%     numParam = 
%     
%    for bandLoop = 1:numBands
%        
%        
%        
%     
%    end
% end


%Number of unknown reflectances + unknown specular coefficients + Area of
%pixel
numParam = numRho + numPix;

%formulate Cx - d
C = zeros(numBands*numPix,numParam);
d = zeros(numBands*numPix, 1);

mCount = 0;

for bandLoop = 1:numBands
    for pixLoop = 1:numPix
        mCount = mCount + 1;
        
        rhoCount = (bandLoop - 1)*numPatch + patchID(pixLoop);
        
        C(mCount, rhoCount) = g_coeff(pixLoop)*s(bandLoop, pixLoop);
        C(mCount, numRho + pixLoop) = dirR_dot_dirC(pixLoop)*s(bandLoop, pixLoop);
        
        d(mCount) = I(bandLoop, pixLoop);
    end
end

xOpt = lsqlin(C, d, -eye(numParam,numParam), zeros(numParam,1), [], [], zeros(numParam,1), ones(numParam,1), 0.5.*ones(numParam, 1), optOptions);
% xOpt = lsqlin(C, d, [], [], [], [], zeros(numParam,1), ones(numParam,1), 0.5.*ones(numParam, 1), optOptions);


% xOpt = lsqlin(C, d, [], [], [], [], zeros(numParam, 1), ones(numParam, 1), 0.5.*ones(numParam, 1), optOptions);
% xOpt = C\d;
% xOpt = lsqlin(C, d);

rhoOpt = xOpt(1:numRho);

rho = zeros(numBands, numPix);
% rhoCount = 0;

for bandLoop = 1:numBands
    for pixLoop = 1:numPix
        
        rhoCount = (bandLoop - 1)*numPatch + patchID(pixLoop);
        
        
%         if ~mod(pixLoop-1, patchSize)
%             rhoCount = rhoCount + 1;
%         end
        
        rho(bandLoop, pixLoop) = rhoOpt(rhoCount);
    end
end


end

