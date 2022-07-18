%Initial script for simulating a non-isotropic point light source.

clc;
close all;
clear;

%external library directory
addpath('ext_lib');

%robotics toolbox
run(fullfile('ext_lib', 'rvctools', 'startup_rvc.m'));

%yaml reader package
addpath(genpath(fullfile('ext_lib', 'yamlmatlab-master')));

%GPML toolbox
run(fullfile('ext_lib', 'gpml-matlab-master', 'startup.m'));

%Better export of figures to images
addpath(fullfile('ext_lib', 'export_fig'))

%code for this project
addpath('src_code');

%results directory
resultsDir = "sim_results";
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

%% default settings

%Value of 6 replicates results in paper
rng(6);

%plot defaults
defaultFontSize = 18;
hwRatioDefault = 1;
defaultFont = 'times';

%% Real RID data and plot

%load real rid data
ridData = readmatrix(fullfile('parameter_files', 'real_rid_data.csv'));
thetaRID = ridData(:,2);
radiantIntenRID = ridData(:,1);

%plot RID data polar plot
figRIDPolar = figure('Name', 'Real RID polar plot');

%plot data
polarplot(thetaRID,radiantIntenRID, 'rx', 'MarkerSize',10); hold on;

%fit and plot spline to data
RIDSplineFit = csapi(thetaRID,radiantIntenRID);
thetas = linspace(-pi/2, pi/2, 1000);
radiantIntSpline = fnval(RIDSplineFit, thetas);
polarplot(thetas,radiantIntSpline, '-b', 'LineWidth',1.5);

%have to manually change axes settings for polar plots
ax = gca;
ax.RColor = [0,0,0];
ax.ThetaColor = [0,0,0];
ax.ThetaLim = [-90,90];
ax.ThetaZeroLocation = 'top';
ax.RAxis.Label.String = 'Normalised RID';
ax.RAxis.Label.Interpreter = 'latex';
ax.RAxis.Label.FontName = defaultFont;
ax.RAxis.Label.FontSize = defaultFontSize;
ax.RAxis.Label.Color = [0,0,0];
ax.ThetaAxis.Label.String = 'Angle (\(^{\circ}\))';
ax.ThetaAxis.Label.Position = [0,1.25,0];
ax.ThetaAxis.Label.Interpreter = 'latex';
ax.ThetaAxis.Label.FontSize = defaultFontSize;
ax.ThetaAxis.Label.FontName = defaultFont;
ax.ThetaAxis.Label.Color = [0,0,0];
legend('Measurements', 'Fitted Spline', 'Location', 'southwest');

imgName = fullfile(resultsDir, 'real_rid_polarplot.png');
exportpropergraphic(figRIDPolar, imgName, 0.6);

%% Setting up non-isotropic light source

%Pose of source in world coordinates (frame camera)
locLightSrc = [0.1;0.1;-0.1];
rotLightSrc = eul2rotm(deg2rad([0,-10, 10]), 'ZYX');

T_S_2_F = [rotLightSrc, locLightSrc; 0, 0, 0, 1];

%radiant intensity distribution parameters
maxRadiInt = 10;

%directional component exponent (higher makes the light more laser like) 
mu = 2;

%attentuation
r_d = 1; %effective radius of disk source

%% Setting up light simulator figure

%figure for light simulator
figSim = figure('Name', 'Light Simulator with Frame Camera Origin');
plotCamera('Location', zeros(1,3), 'Orientation', eye(3), 'Size', 0.05, 'AxesVisible', true); hold on;
S = surf(eye(2), 'Visible', 'off');
xlabel('x');
ylabel('y');
zlabel('z');
grid on;
axis equal;

%light source object
lightSrc = LightSimulator(locLightSrc, rotLightSrc, maxRadiInt, r_d, mu, figSim, S, RIDSplineFit);
lightSrc.PlotRID3D(1000);

%% Setting up Line-scan camera

figure(figSim);

%assuming no distortion and no uncertainty in parameters
%variance of pixel intensity noise in line-scan camera measurements
intNoiseSigma = 0.01;

%line-scan intrinsic camera parameters
fyLS  = 800;
v0 = 160;
yPix = 320;
lsIntrinsic = cameraIntrinsics([1, fyLS], [realmin,v0],[yPix, 1]);

rotLS = eul2rotm(deg2rad([90, 0, -10]), 'ZYX');
tLS = [-0.1; 0.1; -0.1];
T_LS_2_F = [rotLS, tLS; 0,0,0,1];
T_F_2_LS = T_LS_2_F \ eye(4);

plotCamera('Location', tLS, 'Orientation', rotLS', 'Size', 0.05, 'AxesVisible', true, 'Color', [0,0,1]); hold on;

%calculate line-scan view-plane FOV
fovLS =  atan((yPix/2)/fyLS)*2;
%minimum/maximum working focal range of the line-scan camera in metres
minRange = 0.2;
maxRange = 0.7;

%Create the polygon shape of the view-plane using the FOV and working focal
%range
theta = linspace(pi/2 - fovLS/2, pi/2 + fovLS/2, 500);

%define the verticies of the shape in the line-scan camera coordinate frame
yPoly = maxRange.*cos(theta);
zPoly = maxRange.*sin(theta);
yPoly = [yPoly, minRange.*cos(theta(end:-1:1))];
zPoly = [zPoly, minRange.*sin(theta(end:-1:1))];
yPoly = [yPoly, yPoly(1)];
zPoly = [zPoly, zPoly(1)];
xPoly = zeros(size(zPoly));
viewPlanePoly = polyshape(yPoly, zPoly);

%transform verticies of view-plane from frame to line-scan coordinate frame
homPts = [xPoly;yPoly;zPoly; ones(size(zPoly))];
xyzTrans = T_LS_2_F*homPts;

%plot plane
x = xyzTrans(1,:);
y = xyzTrans(2,:);
z = xyzTrans(3,:);

vpPatch = patch(x,y,z, [0,0,1], 'FaceAlpha', 0.3);

vpDispCheck = uicontrol('Parent',figSim,'Style','checkbox', 'String', 'Display View-plane', 'Position', [20,25,200,20] , 'Value', true);
vpDispCheck.Callback = @(src, eventData) vp_callback(src, eventData, vpPatch);

%% 2D plot of radiant intensity on view-plane

cab(figSim);

%Pose of the light source w.r.t to the c.f of the line-scan camera
T_S_2_LS = T_F_2_LS * T_S_2_F;

%plot the light intensity produced by the light source on the line-scan
%camera's view-plane
yDist = 0.5; %max distance left/right
zDist = 1; %max distance along z-direction
lsFigYlims = [-0.5,0.5];
lsFigZlims = [T_S_2_LS(3,4),0.8];
irradVarVPMax = 2.5;

%Create mesh
y = linspace(-yDist, yDist, 100);
z = linspace(T_S_2_LS(3,4), zDist, 100);
x = 0;
[X,Y,Z] = meshgrid(x,y,z);
%remove extra unnecessary singular dimension
X = squeeze(X);
Y = squeeze(Y);
Z = squeeze(Z);

%These points are in the line-scan c.f, transform them into
%the frame camera c.f (world coordinates)
pts = [X(:),Y(:),Z(:)]';
ptsHom = [pts; ones(1, size(pts, 2))];
ptsHomFrame = T_LS_2_F*ptsHom;

rows = size(X);
%reshape to mesh
XFrame = reshape(ptsHomFrame(1,:),rows);
YFrame = reshape(ptsHomFrame(2,:),rows);
ZFrame = reshape(ptsHomFrame(3,:),rows);

radIntMag = lightSrc.RadiantIntensityMesh(XFrame, YFrame, ZFrame);

%plot 2D view-plane surface with line-scan camera as origin. Only plotting
%YZ plane.
figViewPlane = figure('Name', 'Irradiance on line-scan view-plane');

[~,cb] = plotlinescanviewplane(Y, Z, radIntMag, maxRadiInt, T_S_2_LS, T_F_2_LS, [yPoly;zPoly], lsFigYlims, lsFigZlims);
cb.Label.String = 'Irradiance';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = defaultFontSize;
cb.Label.FontName = defaultFont;

imgName = fullfile(resultsDir, 'irradiance_slice.png');
exportpropergraphic(figViewPlane, imgName,  hwRatioDefault);

%% Reflective target
% assume that the target is an infinitely tall flat surface which is
% 10cm wide. The origin on the target will be in the middle.
%target will be positioned randonly in the line-scan camera view.
% Z-axis of the target is inline with its normal. i.e. z-axis vector = surface normal
target.Width = 0.2; %width of target in metres
target.Reflectance = 0.95; %reflectance of target for any wavelength
target.LeftEdge = [0; -target.Width/2; 0];%plot the light intensity produced by the light source on the line-scan
target.RightEdge = [0; target.Width/2; 0];

%% Sampling using the target
%Target will be located in random locations in the FOV of the line-scan
%camera to capture light readings.
%The parameters that need to be varied are the following
% - distance from the optical centre of line-scan to the coordinate frame of the target within the FOV of the line-scan
% - the pixel location of the centre of the target in the line-scan image
% - orientation of target or the normal of the target surface.

%different random locations of the target in the view-plane of the
%line-scan camera
noSamples = 50;

KMat = lsIntrinsic.IntrinsicMatrix';
KMat(1,3) = 0;

%store the points in the frame camera c.f and the normal for each sample
targetSamplePtsFrame = cell(1,noSamples);
targetSampleSurfNorm = zeros(3, noSamples);


minTheta = deg2rad(-10);
maxTheta = deg2rad(60);

for sampleLoop = 1:noSamples
    %get random distance between line-scan and target in metres
    distLS2Target = (maxRange - minRange)*rand() + minRange;
    %plot the light intensity produced by the light source on the line-scan
    
    %pixel location of the centre of target in line-scan image
    pixTargetCent = randi([0,yPix]);
    
    %normalise centre pixel location
    pixTargetCentHom = [0 ; pixTargetCent; 1];
    normTargetCent = KMat\pixTargetCentHom;
    normTargetCent = normTargetCent./normTargetCent(3);
    
    %Use the normalised coordinate to determine the direction of the pixel
    %ray. The translation of the target coordinate frame is found by
    %multiplying the distance from the line-scan camera with the direction
    %of the ray
    dirNormTargetCent = normTargetCent./norm(normTargetCent);
    targetCentPos =  distLS2Target.*dirNormTargetCent;

    
    %pick angle about x-axis between 90deg and 270 deg to get the pose of
    %the target. Only assume rotation about x-axis (yaw only)
    %                 a = deg2rad(100);
    %                 b = deg2rad(260);
    %                 xAngle = (b-a)*rand()+(a);
    %     %initially assume a constant normal angle
    xAngle = pi;
    rotTargetLS = eul2rotm([0, 0, xAngle], 'ZYX');
    
    %pose of target w.r.t to line-scan camera c.f
    T_target_2_LS = [rotTargetLS, targetCentPos; 0, 0, 0, 1];

    %normal vector of the target surface (z-axis) w.r.t to LS
    normTargetLS = T_target_2_LS(1:3,3);
    
    %left and right edge of the target in its coordinate frame
    targetEdgesPtsHom = [[target.LeftEdge; 1], [target.RightEdge; 1]];
    
    %project target edges to line-scan image
    imgTargetEdges = projectPoints(targetEdgesPtsHom', KMat, T_target_2_LS, [], lsIntrinsic.ImageSize);
    
    %extract v coordinate of pixels
    vImgTargetEdges = imgTargetEdges(:,2);
    
    %Check if the edge pixels are within the size of the line-scan image, else
    %set the edge pixels to be the limits of the line-scan image
    for j = 1:2
        vImgTargetEdges(j) = round(vImgTargetEdges(j));
        
        %outside on the left of the image line
        if vImgTargetEdges(j) < 1
            vImgTargetEdges(j) = 1;
            %outside on the right of the image line
        elseif vImgTargetEdges(j) > yPix
            vImgTargetEdges(j) = yPix;
        end
    end
    
    
    %vector of all pixels that see the target
    if vImgTargetEdges(1) > vImgTargetEdges(2)
        vtargetPix = vImgTargetEdges(2):vImgTargetEdges(1);
    else
        vtargetPix = vImgTargetEdges(1):vImgTargetEdges(2);
    end
    
    %transform to normalized coordinates
    targetPixHom  = [zeros(size(vtargetPix)); vtargetPix; ones(size(vtargetPix))];
    normTargetPixHom = KMat\targetPixHom;
    
    %3D points on target plane
    targetPtsLS = zeros(3, length(vtargetPix));
    ptCount = 0;
    
    %trace pixel ray from the line-scan camera to the target frame to
    %determine its 3D location w.r.t LS
    for pixelLoop = 1:length(vtargetPix)
        pnt = normTargetPixHom(:,pixelLoop);
                
        [pntOnTarget, valid] = CameraPixelRayIntersectPlane(pnt', normTargetLS', targetCentPos');
        
        %ray intersection was not valid
        if ~valid
           continue; 
        end
        
        ptCount = ptCount + 1;
        targetPtsLS(:, ptCount) = pntOnTarget;
    end
    
    targetPtsLS = targetPtsLS(:,1:ptCount);
    
    %transform points from line-scan to frame camera
    targetPtsFrame = T_LS_2_F*[targetPtsLS; ones(1,length(targetPtsLS(1,:)))];
    targetPtsFrame = targetPtsFrame(1:3,:);
    
    %save the target pts in the c.f of the frame camera
    targetSamplePtsFrame(sampleLoop) = {targetPtsFrame};
    
    %pose of target w.r.t frame camera c.f
    T_target_2_F = T_LS_2_F*T_target_2_LS;
    %surface normal of target w.r.t to F
    targetSampleSurfNorm(:,sampleLoop) = T_target_2_F(1:3,3);
    
    
    %plot line and axes of target in the simulator figure
    figure(figSim);
    line(targetPtsFrame(1,:), targetPtsFrame(2,:), targetPtsFrame(3,:), 'Color', [0,0,0], 'LineWidth', 2); hold on;
    trplot(T_target_2_F, 'rviz', 'length', (maxRange - minRange)*0.05);
    
    %plot line in the view-plane intensity plot
    figure(figViewPlane);
    line(targetPtsLS(2,:), targetPtsLS(3,:), maxRadiInt.*ones(size(targetPtsLS(2,:))), 'Color', [0,0,0], 'LineWidth', 2);    
end

drawnow();

%% Measure pixel intensity at the 3D location on the target which is relative to the frame camera

%stores the intensity measured at a specific 3D location. The measured
%intensity at the symmetric target would be the exact same
targetSampleNormRadiance = cell(1, noSamples);

for sampleLoop = 1:noSamples
    targetPtsFrame = targetSamplePtsFrame{sampleLoop};
    targetNormal = targetSampleSurfNorm(:, sampleLoop);
    
    %intensity measured at the 3D location by the line-scan camera in the
    %coordinate frame of the frame camera.
    %***** NOTE: Darkening effects, vignetting, sensor irregularties have
    %not been considered yet
    
    numPts = size(targetPtsFrame, 2);
    
    targetInten = zeros(1,numPts);
    
    for pt = 1:numPts
        targetInten(pt) = lightSrc.RadianceOutMaterialPoint(targetPtsFrame(:,pt), targetNormal, target.Reflectance);
    end
    

    %add gaussian white noise to target pixel measurements
    targetIntenNoise = targetInten + normrnd(0, intNoiseSigma, size(targetInten));
        
    targetSampleNormRadiance(sampleLoop) = {targetIntenNoise};
end

%% Combine the samples such that they can be analysed by regressional approaches

%number of trials is equal to the number of samples. Each consective trial
%will add an extra sample

% {noSamples} x {targetL, targetPntNormals, targetPntFrame, targetPntNormalsSymm, targetPntFrameSymm}
targetTrials = cell(noSamples, 3);

for trials = 1:noSamples
    targetL = zeros(1,2*trials*yPix);
    targetPntNormals = zeros(3,2*trials*yPix);
    targetPntFrame = zeros(3, 2*trials*yPix);
    
    noPts = 0;
    
    %combines all samples into a single array for the current trials
    for sampleLoop = 1:trials
        noCurPts = length(targetSampleNormRadiance{sampleLoop});
        
        targetL(noPts+1:noPts+noCurPts) = targetSampleNormRadiance{sampleLoop};
        targetPntNormals(:, noPts+1:noPts+noCurPts) = repmat(targetSampleSurfNorm(:, sampleLoop), [1, noCurPts]);
        targetPntFrame(:, noPts+1:noPts+noCurPts) = targetSamplePtsFrame{sampleLoop};

        noPts = noPts + noCurPts;
    end
    
    %clip the arrays to correct size
    targetL = targetL(1:noPts);
    targetPntNormals = targetPntNormals(:, 1:noPts);
    targetPntFrame = targetPntFrame(:, 1:noPts);
    
    targetTrials(trials,:) = {targetL, targetPntNormals, targetPntFrame};    
end

%% Estimate parameters through least squares

optOptions = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'SpecifyObjectiveGradient',true, 'CheckGradients', false, ...
    'MaxIterations', 1000000000, 'FunctionTolerance',1e-6, 'MaxFunctionEvaluations',1000000000, 'StepTolerance',1e-6, ...
    'FiniteDifferenceType', 'central', 'ScaleProblem','none');
optOptions.Display = 'none';

optPhiTrials = zeros(noSamples, 3);
resNormTrials = zeros(noSamples, 1);

eigFIM = zeros(noSamples, 3);
rankFIM = zeros(noSamples, 1);

Phi0 = [1, 1, 1];

for trials = 1:noSamples
    
    targetL = targetTrials{trials,1};
    targetPntFrame = targetTrials{trials,3};
    targetPntNormals = targetTrials{trials,2};
    
    [optPhi,resNorm, curRankFIM, curEigFIM] = LightSrcOptmLS(lightSrc, Phi0, targetL, targetPntFrame, targetPntNormals, target.Reflectance, optOptions, intNoiseSigma);
    
    %store optimised parameters and residual
    optPhiTrials(trials, :) = optPhi;
    resNormTrials(trials) = resNorm;
    
    rankFIM(trials) = curRankFIM;
    eigFIM(trials,:) = curEigFIM;
end

%% Results of least squares

%close all figures except for the simulator and the view-plane ground truth
%cab(figSim, figViewPlane, figRIDPolar);

%plot residual
figure('Name', 'Optimisation Residuals');
plot(1:noSamples, resNormTrials);
xlabel('samples');
ylabel('Residual after optimisation');
ylim([0, max(resNormTrials)])
grid on;

%plot eigenvalues of FIM
figFIM = figure('Name', 'Eigenvalues of NLS Formulation');
plot(1:noSamples, eigFIM(:,1), 'Color', [1,0,0], 'LineWidth', 1); hold on;
plot(1:noSamples, eigFIM(:,2), 'Color', [0,1,0], 'LineWidth', 1); hold on;
plot(1:noSamples, eigFIM(:,3), 'Color', [0,0,1], 'LineWidth', 1); hold on;
legend('\phi0', '\mu', 'r_d')
xlabel('samples');
ylabel('Eigenvalue Magnitude');
ylim([min(eigFIM, [], 'all'), max(eigFIM, [], 'all')])
grid on;

%Compare Radiant intensity view-plane image of ground-truth parameters and
%estimated parameters

lightSrcLeastSqrs = LightSimulator(locLightSrc, rotLightSrc, optPhiTrials(end, 1), optPhiTrials(end, 2), optPhiTrials(end, 3));

%plot the light intensity produced by the light source on the line-scan
%camera's view-plane
radIntMagGroundtruth = lightSrc.RadiantIntensityMesh(XFrame, YFrame, ZFrame);
radIntMagLeastSqr = lightSrcLeastSqrs.RadiantIntensityMesh(XFrame, YFrame, ZFrame);
diffRad = abs(radIntMagGroundtruth - radIntMagLeastSqr);

hfig = figure();
plotlinescanviewplane(Y, Z, radIntMagLeastSqr, maxRadiInt, T_S_2_LS, T_F_2_LS, [yPoly;zPoly], lsFigYlims, lsFigZlims);
imgName = fullfile(resultsDir, 'sim_ls_a.png');
exportpropergraphic(hfig, imgName,  hwRatioDefault);

hfig = figure();
plotlinescanviewplane(Y, Z, diffRad, maxRadiInt, T_S_2_LS, T_F_2_LS, [yPoly;zPoly], lsFigYlims, lsFigZlims, bone(1000));
imgName = fullfile(resultsDir, 'sim_ls_c.png');
exportpropergraphic(hfig, imgName,  hwRatioDefault);


%% GP Testing Points

%transform points from frame camera C.F to light source C.F
ptsHomLightSrc = T_S_2_F\ptsHomFrame;
ptsLightSrc = ptsHomLightSrc(1:3, :);

[ptsRadius,ptsTheta] = cart2sphZ(ptsLightSrc(1,:), ptsLightSrc(2,:), ptsLightSrc(3,:));
testingX = [ptsRadius', ptsTheta'];

%% GP Training Points
    
targetL = targetTrials{end,1};
targetPntFrame = targetTrials{end,3};
targetPntNormals = targetTrials{end,2};

%calculate the direction light vector (point to light source)
pntLightVec = locLightSrc - targetPntFrame;
dirPntLightVec = pntLightVec./vecnorm(pntLightVec);

%Calculate radiant intensity magnitude used for building model
radIntMagPnt = (targetL.*pi)./(target.Reflectance.*dot(targetPntNormals, dirPntLightVec,1));

%find the point w.r.t to light source c.f
targetPntFrameHom = [targetPntFrame; ones(1, length(targetPntFrame(1,:)))];
targetPntLightSrcHom = T_S_2_F\targetPntFrameHom;
targetPntLightSrc = targetPntLightSrcHom(1:3, :);

[ptsRadius,ptsTheta] = cart2sphZ(targetPntLightSrc(1,:), targetPntLightSrc(2,:), targetPntLightSrc(3,:));
gpTrainingdata = [ptsRadius', ptsTheta', radIntMagPnt'];
downSamplingGP = 30;

%% Building model with GP Zero-mean non-symmetric

%Build zero-mean GP model
[mu, varMu] = LightSrcOptmGP(testingX, gpTrainingdata, 0, intNoiseSigma, downSamplingGP, 10000, false);

radIntGP = reshape(mu,rows);
radIntGPZero = radIntGP;
radVar = reshape(varMu,rows);
diffRad = abs(radIntGP - radIntMagGroundtruth);

hfig = figure();
plotlinescanviewplane(Y, Z, radIntGP, maxRadiInt, T_S_2_LS, T_F_2_LS, [yPoly;zPoly], lsFigYlims, lsFigZlims);
imgName = fullfile(resultsDir, 'sim_gp_z_a.png');
exportpropergraphic(hfig, imgName,  hwRatioDefault);

hfig = figure();
plotlinescanviewplane(Y, Z, radVar, irradVarVPMax, T_S_2_LS, T_F_2_LS, [yPoly;zPoly], lsFigYlims, lsFigZlims, turbo(1000));
imgName = fullfile(resultsDir, 'sim_gp_z_b.png');
exportpropergraphic(hfig, imgName,  hwRatioDefault);

hfig = figure();
plotlinescanviewplane(Y, Z, diffRad, maxRadiInt, T_S_2_LS, T_F_2_LS, [yPoly;zPoly], lsFigYlims, lsFigZlims, bone(1000));
imgName = fullfile(resultsDir, 'sim_gp_z_c.png');
exportpropergraphic(hfig, imgName,  hwRatioDefault);

%% Building model with GP non-Symmetric Constant Mean

%Build zero-mean GP model
[mu, varMu] = LightSrcOptmGP(testingX, gpTrainingdata, 1, intNoiseSigma, downSamplingGP, 10000, false);

radIntGP = reshape(mu,rows);
radVar = reshape(varMu,rows);
diffRad = abs(radIntGP - radIntMagGroundtruth);

hfig = figure();
plotlinescanviewplane(Y, Z, radIntGP, maxRadiInt, T_S_2_LS, T_F_2_LS, [yPoly;zPoly], lsFigYlims, lsFigZlims);
imgName = fullfile(resultsDir, 'sim_gp_const_a.png');
exportpropergraphic(hfig, imgName,  hwRatioDefault);

hfig = figure();
plotlinescanviewplane(Y, Z, radVar, irradVarVPMax, T_S_2_LS, T_F_2_LS, [yPoly;zPoly], lsFigYlims, lsFigZlims, turbo(1000));
imgName = fullfile(resultsDir, 'sim_gp_const_b.png');
exportpropergraphic(hfig, imgName,  hwRatioDefault);

hfig = figure();
plotlinescanviewplane(Y, Z, diffRad, maxRadiInt, T_S_2_LS, T_F_2_LS, [yPoly;zPoly], lsFigYlims, lsFigZlims, bone(1000));
imgName = fullfile(resultsDir, 'sim_gp_const_c.png');
exportpropergraphic(hfig, imgName,  hwRatioDefault);


%% Building model with GP Non-Symmetric light source mean function

%Build zero-mean GP model
[mu, varMu] = LightSrcOptmGP(testingX, gpTrainingdata, 3, intNoiseSigma, downSamplingGP, 10000, false);

radIntGP = reshape(mu,rows);
radVar = reshape(varMu,rows);
diffRad = abs(radIntGP - radIntMagGroundtruth);

hfig = figure();
plotlinescanviewplane(Y, Z, radIntGP, maxRadiInt, T_S_2_LS, T_F_2_LS, [yPoly;zPoly], lsFigYlims, lsFigZlims);
imgName = fullfile(resultsDir, 'sim_gp_lightsrc_a.png');
exportpropergraphic(hfig, imgName,  hwRatioDefault);

hfig = figure();
plotlinescanviewplane(Y, Z, radVar, irradVarVPMax, T_S_2_LS, T_F_2_LS, [yPoly;zPoly], lsFigYlims, lsFigZlims, turbo(1000));
imgName = fullfile(resultsDir, 'sim_gp_lightsrc_b.png');
exportpropergraphic(hfig, imgName,  hwRatioDefault);

hfig = figure();
plotlinescanviewplane(Y, Z, diffRad, maxRadiInt, T_S_2_LS, T_F_2_LS, [yPoly;zPoly], lsFigYlims, lsFigZlims, bone(1000));
imgName = fullfile(resultsDir, 'sim_gp_lightsrc_c.png');
exportpropergraphic(hfig, imgName,  hwRatioDefault);

%% 

%calculate absolute differences w.r.t to ground-truth
diffRadLS = abs(radIntMagGroundtruth - radIntMagLeastSqr);
diffRadGP = abs(radIntMagGroundtruth - radIntGP);
diffRadZero = abs(radIntMagGroundtruth - radIntGPZero);

%vectorise 2d matrices
Yvec = Y(:);
Zvec = Z(:);
diffLSvec = diffRadLS(:);
diffGPvec = diffRadGP(:);
diffGPZvec = diffRadZero(:);

%view-plane polygon verticies
YpolyVP = viewPlanePoly.Vertices(:,1);
ZpolyVP = viewPlanePoly.Vertices(:,2);

%determine if points are inside/outside view-plane
ptsInVPMask = inpolygon(Yvec, Zvec, YpolyVP, ZpolyVP);

%separate inside/outside points
diffLS_in_VP = diffLSvec(ptsInVPMask);
diffLS_out_VP = diffLSvec(~ptsInVPMask);
diffGP_in_VP = diffGPvec(ptsInVPMask);
diffGP_out_VP = diffGPvec(~ptsInVPMask);
diffGPZ_out_VP = diffGPZvec(~ptsInVPMask);
diffGPZ_in_VP = diffGPZvec(ptsInVPMask);


radDiffCell = {diffLS_in_VP, diffLS_out_VP,  diffGPZ_in_VP, diffGPZ_out_VP, diffGP_in_VP, diffGP_out_VP};
varLabel = {'Inside', 'Outside', 'Inside', 'Outside', 'Inside', 'Outside'};


hfig = figure();

%violin plot colours
violinCols = [
0.5	0.5	1.0;
0.5	0.5	1.0;
0	1	1;
0	1	1;
1	1	0;
1	1	0;
];

violin(radDiffCell, 'xlabel', varLabel, 'ylabel', 'Fitted Absolute Difference', 'facecolor', violinCols, 'facealpha', 1); hold on;
ylim([-inf, maxRadiInt])
ylabel('Fitted Absolute Difference');
set(gca, 'YGrid', 'on', 'XGrid', 'off')
annotation('textbox',...
    [0.23 0.015 0.06 0.04],...
    'String',{'Least Squares'},...  
    'FitBoxToText','on', 'EdgeColor', 'none', 'HorizontalAlignment','center');

annotation('textbox',...
    [0.49 0.015 0.06 0.04],...
    'String',{'GP Zero-Mean'},...  
    'FitBoxToText','on', 'EdgeColor', 'none', 'HorizontalAlignment','center');

annotation('textbox',...
    [0.75 0.015 0.06 0.04],...
    'String',{'GP Light-Src-Mean'},...  
    'FitBoxToText','on', 'EdgeColor', 'none', 'HorizontalAlignment','center');


imgName = fullfile(resultsDir, 'violin.png');
exportpropergraphic(hfig, imgName,  0.7);

%% CALLBACK functions

function vp_callback(src, ~, vpPatch)
value = get(src,'Value');

if value
    vpPatch.Visible = true;
else
    vpPatch.Visible = false;
end

end