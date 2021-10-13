function [S, gImg, kImg] = Estimate_Sgk_WithShape(R, pixList, maskID, lDotnImg, neiPix, rSAM)
% Estimate the reflectance, shading factor, and specularity coefficient
% of an input radiance hypercube (that has lighting removed) using Krebs
% method with added shape information. Solve the specularity
% coefficient and shading factor independently through quadratic
% programming formulation.
%
% INPUTS:
%       R - radiance hypercube with lighting removed[nY x mX x nB]
%       pixList - list of relevant pixels where the shading factor needs to
%           be estimated [num pixels x 2]
%       maskID - index of relevant pixels in pixList at their exact pixel
%           location in the hypercube [nY x nX]
%       lDotnImg - Img of dot product between light direction vector from
%       source and surface normal [nY x nX]
%       neiPix - max pixel distance that is considered apart of neighbourhood
%       rSAM - bandwidth parameter for ZetaSAM
% OUTPUTS:
%       S - Reflectance hypercube from optimisation [nY x mX x nB]
%       gImg - shading factor image [nY x mX]
%       kImg - specular coefficient image [nY x mX]
%
% Author: Jasprabhjit Mehami, 13446277

[nY, nX, nB] = size(R);

%STD over bands
sigmaR = std(R,0,3);
ln_sigmaR = log(sigmaR);

%mean over bands
muR = mean(R, 3);
ln_muR = log(muR);

%log of dotproduct image
ln_lDotnImg = log(lDotnImg);

%Convolution kernel applied to extract neighbourhood pixels without centre
%pixel
extractKern = ones(neiPix);
extractKern(ceil(numel(extractKern)/2)) = 0;

%mask which stores the location of the pixel of interest.
pixLocMask = zeros(nY, nX);

%number of relevant pixels
numX = size(pixList,1);

%Shading factor vectors for quadratic programming
fg = zeros(numX,1);
Ig = zeros(numX*neiPix^3,1);
Jg = zeros(numX*neiPix^3,1);
Vg = zeros(numX*neiPix^3,1);
ctg = 0;

%specular coefficient vectors for quadratic programming
fk = zeros(numX,1);
Ik = zeros(numX*neiPix^3,1);
Jk = zeros(numX*neiPix^3,1);
Vk = zeros(numX*neiPix^3,1);
ctk = 0;



%zetaSAM of each pixel with a vector of 1 (how "gray" the pixel is?
onesB = ones(1,nB);
z_1 = zeros(numX,1);

for p = 1:numX
    i = pixList(p,1);
    j = pixList(p,2);
    uR = squeeze(R(i,j,:))';
    z_1(p) = ZetaSAM(uR, onesB, rSAM);
end

%loop through all relevant pixels
for p = 1:numX
    i = pixList(p,1);
    j = pixList(p,2);
    
    %set current location to extract its neighbourhood
    pixLocMask(i,j) = 1;
    u_id = p;
    
    %list of neighbourhood pixel IDs
    neighID = maskID(conv2(pixLocMask,extractKern,'same')>0);
    neighID(neighID < 1) = []; %remove zero ID pixels
    
    uR = squeeze(R(i,j,:))';
    s_u = ln_sigmaR(i,j);
    m_u = ln_muR(i,j);
    d_u = ln_lDotnImg(i,j);
    z_u = z_1(p);
    
    %loop through all neighbours
    for k = 1:length(neighID)
        v_id = neighID(k); %current neighbour ID
        pixLoc = pixList(v_id,:); %current neighbour pixel location
        
        vR = squeeze(R(pixLoc(1),pixLoc(2),:))';
        s_v = ln_sigmaR(pixLoc(1),pixLoc(2));
        m_v = ln_muR(pixLoc(1),pixLoc(2));
        d_v = ln_lDotnImg(pixLoc(1),pixLoc(2));
        
        % Calculate zetaSAM between pixels
        z_uv = ZetaSAM(uR, vR, rSAM);
        z_v = z_1(v_id);
        
        %*******Shading Factor**********
        
        %Calculate the coefficients of the quadratic formulation for
        %shading factor for the current u and v pixel
        coef = ShadingFactorCoeffsUV(s_u, s_v, m_u, m_v, d_u, d_v, z_uv, z_u, z_v);
        
        % u^2 term
        ctg = ctg + 1;
        Ig(ctg) = u_id;
        Jg(ctg) = u_id;
        Vg(ctg) = coef(1);
        
        % v^2 term
        ctg = ctg + 1;
        Ig(ctg) = v_id;
        Jg(ctg) = v_id;
        Vg(ctg) = coef(2);
        
        % u*v term
        ctg = ctg + 1;
        Ig(ctg) = u_id;
        Jg(ctg) = v_id;
        Vg(ctg) = coef(3);
        
        % v*u term
        ctg = ctg + 1;
        Ig(ctg) = v_id;
        Jg(ctg) = u_id;
        Vg(ctg) = coef(3);
        
        fg(u_id) = fg(u_id) + coef(4); % u term
        fg(v_id) = fg(v_id) + coef(5); % v term
        
        
        %******Specular Coefficient*********
        
        %Calculate the coefficients of the quadratic formulation for
        %specular coefficients for the current u and v pixel
        coef = SpecularIntermediateCoeffsUV(s_u, s_v, m_u, m_v, d_u, d_v, z_uv, z_u, z_v);
        
        % u^2 term
        ctk = ctk + 1;
        Ik(ctk) = u_id;
        Jk(ctk) = u_id;
        Vk(ctk) = coef(1);
        
        % v^2 term
        ctk = ctk + 1;
        Ik(ctk) = v_id;
        Jk(ctk) = v_id;
        Vk(ctk) = coef(2);
        
        % u*v term
        ctk = ctk + 1;
        Ik(ctk) = u_id;
        Jk(ctk) = v_id;
        Vk(ctk) = coef(3);
        
        % v*u term
        ctk = ctk + 1;
        Ik(ctk) = v_id;
        Jk(ctk) = u_id;
        Vk(ctk) = coef(3);
        
        fk(u_id) = fk(u_id) + coef(4); % u term
        fk(v_id) = fk(v_id) + coef(5); % v term
        
    end
    
    %unset current location
    pixLocMask(i,j) = 0;
end

%clip arrays to size
Ig = Ig(1:ctg);
Jg = Jg(1:ctg);
Vg = Vg(1:ctg);

Ik = Ik(1:ctk);
Jk = Jk(1:ctk);
Vk = Vk(1:ctk);

%create sparse H matrix
Hg = sparse(Ig,Jg,Vg);
Hk = sparse(Ik, Jk, Vk);

%make H matrix symmetric
Hg = (Hg + Hg')./2;
Hk = (Hk + Hk')./2;

options = optimoptions('quadprog','Algorithm','interior-point-convex',...
    'LinearSolver','sparse','StepTolerance',0, 'MaxIterations', 10000, 'OptimalityTolerance', 1e-10, 'Display', 'none');

%solve for shading factor intermediate variable
x = quadprog(Hg, fg, [], [], [], [], [], [], [], options);

%vector of index position of each relevant pixel location
ind = sub2ind([nY, nX], pixList(:,1)', pixList(:,2)');

muR_vec = muR(ind); %vector of R for all relevant pixels
ln_lDotnVec = ln_lDotnImg(ind); %vector of light source and surface normal dot product for all relevant pixels
R_NAN = R;
R_NAN(R_NAN < 1e-7) = nan; %make any close to zero measurements nan

%minimum measurement for each pixels across the bands
minR = min(R_NAN, [], 3, 'omitnan');
minR_vec = minR(ind); %in vector form

%lower and upper bounds for optimisation
lb = log(muR_vec - minR_vec) - ln_lDotnVec;
ub = log(muR_vec) - ln_lDotnVec;

%solve for specular component intermediate variable
y = quadprog(Hk, fk, [], [], [], [], lb, ub, [], options);

%shading factor and specular coefficient images
gImg = zeros(nY, nX);
kImg = zeros(nY, nX);

%calculate shading factor image from solution of x
gVec = exp(x + ln_lDotnVec');
gImg(ind) = gVec;

%calculate specular cofficient from solution of y
kVec = muR_vec' - exp(y + ln_lDotnVec');
kImg(ind) = kVec;

%calculate reflectance
S = bsxfun(@minus, R, kImg); %subtract specular coefficient
S = bsxfun(@rdivide, S, gImg); % divide by shading factor
S(isinf(S)) = 0; %no inf
S(isnan(S)) = 0; %no nan
S(S < 0 ) = 0; %no negative
end
