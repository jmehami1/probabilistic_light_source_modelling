function imgRadiNorm = NormaliseRadianceBands(imgRad, whiteRef, darkRef, bandStart, bandEnd)
%Calculates the calibrated reflectance of a hyperspectral radiance image
%given the white and dark references through normalisation.
%The output reflectance image will be the same size as the input image.
%All reflectance values will be between 0 and 1.
% INPUTS:
%       imgRad - radiance image
%       whiteRef - white reference image
%       darkRef - dark reference image
%       bandStart - first relevant band (row) in images
%       bandEnd - last relevant band (row) in images
% OUTPUT:
%       imgRadiNorm - reflectance image which is the same size as input
%       radiance image. All of its values are between 0 and 1

%Author: Jasprabhjit Mehami, 13446277

%size of image
[numBands, numPixels] = size(imgRad);

%check if the white reference is brighter than the dark reference for all
%relevant pixel bands
if any(whiteRef(bandStart:bandEnd,:) < darkRef(bandStart:bandEnd, :), 'all')
    imgRadiNorm = -1;
    return;
end

%normalise the radiance image by calculating the numerator and denominator
%separately
numer = imgRad - darkRef;
denom = whiteRef - darkRef;

checkImg = false(numBands, numPixels);

%check that all denominator pixel bands are less than numerator pixel
%bands
checkImg(bandStart:bandEnd,:) = denom(bandStart:bandEnd,:) < numer(bandStart:bandEnd, :);

%if any denominator values are less than numerator values, make the
%numerator values equal to the denominator values
if any(checkImg, 'all')
    numer(checkImg) = denom(checkImg);
end

imgRadiNorm = zeros(numBands, numPixels);
imgRadiNorm(bandStart:bandEnd,:) = double(numer(bandStart:bandEnd, :))./...
    double(denom(bandStart:bandEnd,:));

end

