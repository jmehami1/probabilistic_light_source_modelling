function [S, kImg, gImg] = ReflectanceEstimationKrebs(I, l, mask, lDotnImg, neiPix, rSAM)
% Estimation of reflectance using Krebs methods from the papers:
% 
% Krebs, A., Benezeth, Y., & Marzani, F. (2020). Intrinsic RGB and multispectral 
%   images recovery by independent quadratic programming. PeerJ Computer Science, 
%   6, 1–15. https://doi.org/10.7717/PEERJ-CS.256
% 
% Krebs, A., Benezeth, Y., & Marzani, F. (2017). Quadratic objective functions for
%   dichromatic model parameters estimation. 2017 International Conference on Digital
%   Image Computing: Techniques and Applications (DICTA), 1–8.
% 
% INPUTS:
%       I - radiance hypercube [nY x mX x nB]
%       l - radiant intensity vector [1 x nB] or radiant intensity 
%           hypercube [nY x mX x nB]
%       mask - spaitial pixel mask to only work with relevant pixels [nY x nX]
%       lDotnImg - dot product between light direction vector and surface
%           normal direction at each pixel [nY x nX]
%       neiPix - max pixel distance that is considered apart of neighbourhood
%       rSAM - bandwidth parameter for ZetaSAM
% OUTPUTS:
%       S - Reflectance hypercube from optimisation [nY x mX x nB]
%       gImg - shading factor image [nY x mX]
%       kImg - specular coefficient image [nY x mX][nY, nX, nB] = size(I);
% 
% Author: Jasprabhjit Mehami, 13446277

[nY_l, nX_l, nB_l] = size(l);
[nY, nX, nB] = size(I);


%calculate the R matrix which removes the lighting from radiance
if (nY_l == 1) && (nX_l == nB) && (nB_l == 1)
    % radiant intensity vector
    R = I;
    
    for i = 1:nB
       R(:,:,i) = R(:,:,i)./l(i);
    end
%passed in radiant intensity hypercube    
elseif (nY_l == nY) && (nX_l == nX) && (nB_l == nB)
    % radiant intensity hypercube     
    R = I./l;
else
    error('input radiant intensity not of correct size');
end

%Turn NaN elements to zero
R(isnan(R)) = 0;

%list of relevant pixels given by mask. Each of these pixels has an unknown
%shading factor, reflectance, and specular coefficient
[rows, cols] = find(mask);
pixList = [rows,cols];

%store the index of each relevant pixel from pixList into a ID mask. 
pixImgID = zeros(size(mask));
for i = 1:size(pixList,1)
    pixImgID(pixList(i, 1), pixList(i, 2)) = i;
end

%Estimate reflectance, shading factor, and specular coefficient. Note that
%there is an unknown balancing factor between the S and g.
[S, gImg, kImg] = Estimate_Sgk_WithShape(R, pixList, pixImgID, lDotnImg, neiPix, rSAM);

end
