function alpha = EstimateBalancingFactor(I, l, lDotnImg, refl)

neiPix = 3;

rSAM = pi*pi/3600;

R = I./l;

[nY, nX, nB] = size(R);

%number of relevant pixels
numX = nY;

maskID = zeros(nY, nX);

for i = 1:numX
   maskID(i) = i; 
end

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



%Shading factor vectors for quadratic programming
fg = zeros(numX,1);
Ig = zeros(numX*neiPix^3,1);
Jg = zeros(numX*neiPix^3,1);
Vg = zeros(numX*neiPix^3,1);
ctg = 0;

%zetaSAM of each pixel with a vector of 1 (how "gray" the pixel is?
onesB = ones(1,nB);
z_1 = zeros(numX,1);

for p = 1:numX
    uR = squeeze(R(p,1,:))';
    z_1(p) = ZetaSAM(uR, onesB, rSAM);
end


%loop through all relevant pixels
for p = 1:numX
    i = p;
    j = 1;
    
    %set current location to extract its neighbourhood
    pixLocMask(i,j) = 1;
    %     u_id = maskID(i,j);
    
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
        pixLoc = [v_id,1]; %current neighbour pixel location
        
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
        Ig(ctg) = p;
        Jg(ctg) = p;
        Vg(ctg) = coef(1);
        
        % v^2 term
        ctg = ctg + 1;
        Ig(ctg) = v_id;
        Jg(ctg) = v_id;
        Vg(ctg) = coef(2);
        
        % u*v term
        ctg = ctg + 1;
        Ig(ctg) = p;
        Jg(ctg) = v_id;
        Vg(ctg) = coef(3);
        
        % v*u term
        ctg = ctg + 1;
        Ig(ctg) = v_id;
        Jg(ctg) = p;
        Vg(ctg) = coef(3);
        
        fg(p) = fg(p) + coef(4); % u term
        fg(v_id) = fg(v_id) + coef(5); % v term
        
    end
    
    %unset current location
    pixLocMask(i,j) = 0;
end

%clip arrays to size
Ig = Ig(1:ctg);
Jg = Jg(1:ctg);
Vg = Vg(1:ctg);

%create sparse H matrix
Hg = sparse(Ig,Jg,Vg);

%make H matrix symmetric
Hg = (Hg + Hg')./2;

options = optimoptions('quadprog','Algorithm','interior-point-convex',...
    'LinearSolver','sparse','StepTolerance',0, 'MaxIterations', 10000, 'OptimalityTolerance', 1e-10, 'Display', 'final-detailed');

%solve for shading factor intermediate variable
x = quadprog(Hg, fg, [], [], [], [], [], [], [], options);

%calculate shading factor image from solution of x
gVec = exp(x + ln_lDotnImg);

S = refl.*ones(nB, nY);
gS = gVec.*S';
T = squeeze(R)./gS;
T = T.*pi;
T = T./lDotnImg;

alpha = mean(T, 'all');

% alpha = gVec.*pi;
% alpha = alpha ./lDotnImg;
end

