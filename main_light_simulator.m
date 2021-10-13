%Initial script for simulating a non-isotropic point light source.

close all;
clear;

%external library directory
addpath('ext_lib');

%robotics toolbox
run(['ext_lib', filesep, 'rvctools', filesep, 'startup_rvc.m']);

%yaml reader package
addpath(genpath(['ext_lib', filesep, 'yamlmatlab-master']));

%GPML toolbox
run(['ext_lib', filesep, 'gpml-matlab-master', filesep, 'startup.m']);

%code for this project
addpath('src_code');

%% Parameters to change




%% Real RID data and plot
%load real rid data
data = readmatrix('parameter_files/real_rid_data.csv');
thetaRID = data(:,2);
radiantIntenRID = data(:,1);

figRIDPolar = figure('Name', 'Real RID');
polarplot(thetaRID,radiantIntenRID, 'rx'); hold on;
ax = gca;
ax.ThetaLim = [-90,90];
ax.ThetaZeroLocation = 'top';
hold on;
ax.RAxis.Label.String = 'Normalised Radiant Intensity';
% ax.ThetaAxis.Label.String = 'Angle (Deg)';
title('Real RID with Fitted Spline');

RIDSplineFit = csapi(thetaRID,radiantIntenRID);
thetas = linspace(-pi/2, pi/2, 1000);

radiantIntSpline = fnval(RIDSplineFit, thetas);
polarplot(thetas,radiantIntSpline, '-b');

legend('Data Points', 'Fitted Spline', 'Location', 'southwest');

%% Setting up non-isotropic light source

%Pose of source in world coordinates (frame camera)
locLightSrc = [0.1;0.1;-0.1];
rotLightSrc = eul2rotm(deg2rad([0,-10, 0]), 'ZYX');

% locLightSrc = [0.18; -0.05; 0.05];
% rotLightSrc = [0.792897128208011,0.00375573235981344,-0.609343941098892;0.00375573235981344,0.999931891212147,0.0110502222303317;0.609343941098892,-0.0110502222303317,0.792829019420159];
T_S_2_F = [rotLightSrc, locLightSrc; 0, 0, 0, 1];

%radiant intensity distribution parameters
maxRadiInt = 10;
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

%% Setting up Line-scan camera

%assuming no distortion and no uncertainty in parameters
maxPixInt = 256;
darkPixInt = 10;
respCurveGrad = 10;

%variance of pixel intensity noise in line-scan camera measurements
intNoiseSigma = 0.01;
% 
% maxRadiLS = 1.5*maxRadiInt;
% darkRadLS = 0.1;

%line-scan intrinsic camera parameters
fyLS  = 800;
v0 = 160;
yPix = 320;

lsIntrinsic = cameraIntrinsics([1, fyLS], [realmin,v0],[yPix, 1]);

rotLS = eul2rotm(deg2rad([90, 0, -10]), 'ZYX');
% rotLS = rotLightSrc;
tLS = [-0.1; 0.1; -0.1];
% tLS = locLightSrc + [0; 0.1; 0];



% tLS = [0.125101976645449;-0.0169881345205247;-0.110138876985863];
% rotLS = [-0.970289293330054,-0.0609119653632981,-0.234154691869808;0.0563354557071021,-0.998068323042120,0.0261904366165156;-0.235297691614978,0.0122210889642028,0.971846511186408];

%pose of line-scan w.r.t frame camera
% rotLS = rotLightSrc;
% tLS = [0.125101976645449;-0.0169881345205247;-0.10138876985863];

T_LS_2_F = [rotLS, tLS; 0,0,0,1];
T_F_2_LS = T_LS_2_F \ eye(4);
 
% tLS = tform2trvec(poseLS);
% rotLS = tform2rotm(extLS);

% poseLS = [rotLS, tLS; 0, 0, 0, 1];
% extLS = poseLS \ eye(4);

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

%normal vector of the line-scan view-plane in its coordinate frame (normal
%vector is inline with the x-axis)
vpNormal = [1;0;0];

%% 2D plot of radiant intensity on view-plane

%Pose of the light source w.r.t to the c.f of the line-scan camera
T_S_2_LS = T_F_2_LS * T_S_2_F;

%plot the light intensity produced by the light source on the line-scan
%camera's view-plane
yDist = 0.5; %max distance left/right
zDist = 1; %max distance along z-direction

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
figViewPlane = figure('Name', 'Radiant intensity on line-scan view-plane');
surf(Y, Z, radIntMag, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('View-plane RADIANT INTENSITY')
view(2);
scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(hot(1000));
colorbar; caxis([0, maxRadiInt]);
axis equal;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);

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
                
        [pntOnTarget, valid] = camerapixelrayintersectplane(pnt', normTargetLS',  targetCentPos');
        
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
    drawnow();
    
    %plot line in the view-plane intensity plot
    figure(figViewPlane);
    line(targetPtsLS(2,:), targetPtsLS(3,:), maxRadiInt.*ones(size(targetPtsLS(2,:))), 'Color', [0,0,0], 'LineWidth', 2);    
    drawnow();
end

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
    
%     targetNormRadiance = (respCurveGrad.*targetInten)./(maxPixInt - darkPixInt) + normrnd(0, intNoiseSigma, size(targetInten));
    
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
cab(figSim, figViewPlane, figRIDPolar);

% %calculate absolute errors in the estimated parameters
% abs_phi0 = abs(optPhiTrials(:,1) - maxRadiInt);
% abs_mu = abs(optPhiTrials(:,2) - mu);
% abs_r_d = abs(optPhiTrials(:,3) - r_d);
% 
% %plot absolute error
% figAbsLS = figure('Name', 'Absolute Erorr in Estimated Parameters');
% subplot(1,3,1);
% plot(1:noSamples, abs_phi0);
% xlabel('samples');
% ylabel('absolute error in \Phi_0');
% ylim([0, max(abs_phi0)])
% grid on;
% 
% subplot(1,3,2);
% plot(1:noSamples, abs_mu);
% xlabel('samples');
% ylabel('absolute error in \mu');
% ylim([0, max(abs_mu)])
% grid on;
% 
% subplot(1,3,3);
% plot(1:noSamples, abs_r_d);
% xlabel('samples');
% ylabel('absolute error in r_d');
% ylim([0, max(abs_r_d)])
% grid on;

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
yDist = 0.5; %max distance left/right
zDist = 1; %max distance along z-direction

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

radIntMagGroundtruth = lightSrc.RadiantIntensityMesh(XFrame, YFrame, ZFrame);
% radNormGroundTruth = (respCurveGrad.*radIntMagGroundtruth)./(maxPixInt - darkPixInt);

radIntMagLeastSqr = lightSrcLeastSqrs.RadiantIntensityMesh(XFrame, YFrame, ZFrame);
% radNormLeastSqr = (respCurveGrad.*radIntMagLeastSqr)./(maxPixInt - darkPixInt);
% diffRad = abs(radNormGroundTruth - radNormLeastSqr);


% radNormMax = (respCurveGrad.*maxRadiInt)./(maxPixInt - darkPixInt);

diffRad = abs(radIntMagGroundtruth - radIntMagLeastSqr);


%plot 2D view-plane surface with line-scan camera as origin. Only plotting
%YZ plane.
figLeastSqrViewPlane = figure('Name', 'Radiant intensity view-plane from least squares');
% s1 = subplot(1,3,1);
% surf(Y, Z, radIntMagGroundtruth, 'EdgeColor', 'none'); hold on;
% xlabel('y');
% ylabel('z');
% title('Ground-truth')
% view(2);
% scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
% colormap(s1, hot);
% caxis([0, maxRadiInt]);
% % colorbar;
% % axis equal;
% 
% %plot light source in the line-scan camera coordinate frame (green)
% scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
% 
% %plot frame camera in line-scan coordinate frame
% scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
% 
% %plot fov of line-scan on view-plane
% plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);



s2 = subplot(1,2,1);
surf(Y, Z, radIntMagLeastSqr, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Least Squares View-Plane')
view(2);
scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s2, hot);
caxis([0, maxRadiInt]);
colorbar;
% axis equal;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);



maxDiff = max(diffRad, [], 'all');

s3 = subplot(1,2,2);
surf(Y, Z, diffRad, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Absolute Difference')
view(2);
scatter3(0, 0, (maxDiff+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s3, bone);
caxis([0, maxDiff]);
colorbar;
% axis equal;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);


%% GP Testing Points

%transform points from frame camera C.F to light source C.F
ptsHomLightSrc = T_S_2_F\ptsHomFrame;
ptsLightSrc = ptsHomLightSrc(1:3, :);

[ptsRadius,ptsTheta] = cart2sphZ(ptsLightSrc(1,:), ptsLightSrc(2,:), ptsLightSrc(3,:));


% [ptsRadius,ptsTheta] = cart2rtlightsrc(ptsLightSrc);

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

% [ptsRadius,ptsTheta] = cart2rtlightsrc(targetPntLightSrc);
[ptsRadius,ptsTheta] = cart2sphZ(targetPntLightSrc(1,:), targetPntLightSrc(2,:), targetPntLightSrc(3,:));

gpTrainingdata = [ptsRadius', ptsTheta', radIntMagPnt'];


downSamplingGP = 50;

%% Building model with GP Zero-mean non-symmetric

cab(figSim, figRIDPolar, figViewPlane, figLeastSqrViewPlane);

plotMaxHeight = 100;

%Build zero-mean GP model
[mu, varMu, hypOpt] = LightSrcOptmGP(0, gpTrainingdata, testingX, downSamplingGP, 1000, intNoiseSigma,false);

figGP_ZeroMean = figure('Name','GP Zero-Mean');
s1 = subplot(1,3,1);
surfGP = surf(Y, Z, zeros(size(Y)), 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Zero-mean GP View-Plane')
view(2);
scatter3(0, 0, (plotMaxHeight+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s1, hot);
caxis([0, maxRadiInt]);
colorbar;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (plotMaxHeight+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (plotMaxHeight+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (plotMaxHeight+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);

s2 = subplot(1,3,2);
surfGP_var = surf(Y, Z, zeros(size(Y)), 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Variance')
view(2);
scatter3(0, 0, (plotMaxHeight+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s2, turbo);
caxis([0, 1]);
colorbar;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (plotMaxHeight+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (plotMaxHeight+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (plotMaxHeight+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);



figure(figGP_ZeroMean);

%Plotting training points on figure
targetPntLSHom = (T_LS_2_F\targetPntFrameHom)';

lsTrainPlotPt = (downsample(targetPntLSHom, downSamplingGP))';
scatter3(lsTrainPlotPt(2,:), lsTrainPlotPt(3,:), (maxRadiInt+1)*ones(size(lsTrainPlotPt(3,:))), 20, [0,1,1], 'filled','MarkerEdgeColor', [0,0,0], 'LineWidth', 0.1);

radIntGP = reshape(mu,rows);
radVar = reshape(varMu,rows);

radIntGPZero = radIntGP; 

%update figure with GP model
surfGP.ZData = radIntGP;
surfGP_var.ZData = radVar;

colormap(s2, turbo);
caxis([0, max(radVar, [], 'all')]);
drawnow();




diffRad = abs(radIntGP - radIntMagGroundtruth);
maxDiff = max(diffRad, [], 'all');

% Absolute difference plot
s3 = subplot(1,3,3);
surf(Y, Z, diffRad, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Absolute Difference')
view(2);
scatter3(0, 0, (maxDiff+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s3, bone);
caxis([0, maxDiff]);
colorbar;
% axis equal;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
drawnow();












%plot 2D view-plane surface with line-scan camera as origin. Only plotting
%YZ plane.
figDiffGP_ZeroMean = figure('Name', 'Radiant intensity view-plane from GP');
% s1 = subplot(1,3,1);
% surf(Y, Z, radIntMagGroundtruth, 'EdgeColor', 'none'); hold on;
% xlabel('y');
% ylabel('z');
% title('Ground-truth')
% view(2);
% scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
% colormap(s1, hot);
% caxis([0, maxRadiInt]);
% % colorbar;
% % axis equal;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);


s2 = subplot(1,2,1);
surf(Y, Z, radIntGP, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Zero-mean GP View-Plane')
view(2);
scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s2, hot);
caxis([0, maxRadiInt]);
colorbar;
% axis equal;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);


diffRad = abs(radIntGP - radIntMagGroundtruth);
maxDiff = max(diffRad, [], 'all');

% Absolute difference plot
s3 = subplot(1,2,2);
surf(Y, Z, diffRad, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Absolute Difference')
view(2);
scatter3(0, 0, (maxDiff+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s3, bone);
caxis([0, maxDiff]);
colorbar;
% axis equal;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
drawnow();


% %% Building model with GP Zero-mean Symmetric
% 
% cab(figSim, figViewPlane, figAbsLS, figLeastSqrViewPlane, figGP_ZeroMean, figDiffGP_ZeroMean);
% 
% 
% %Build zero-mean GP model
% [mu, varMu, hypOpt] = LightSrcOptmGP(0, gpTrainingdataSYM, testingX, downSamplingGP, 1000, false);
% 
% figGP_ZeroMeanSYM = figure('Name','GP Zero-Mean Symmetric');
% s1 = subplot(1,2,1);
% surfGP = surf(Y, Z, zeros(size(Y)), 'EdgeColor', 'none'); hold on;
% xlabel('y');
% ylabel('z');
% title('Zero-mean Symmetric View-Plane GP')
% view(2);
% scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
% colormap(s1, hot);
% caxis([0, maxRadiInt]);
% colorbar;
% 
% %plot light source in the line-scan camera coordinate frame (green)
% scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
% %plot frame camera in line-scan coordinate frame
% scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
% %plot fov of line-scan on view-plane
% plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
% 
% s2 = subplot(1,2,2);
% surfGP_var = surf(Y, Z, zeros(size(Y)), 'EdgeColor', 'none'); hold on;
% xlabel('y');
% ylabel('z');
% title('Variance')
% view(2);
% scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
% colormap(s2, turbo);
% caxis([0, 1]);
% colorbar;
% 
% %plot light source in the line-scan camera coordinate frame (green)
% scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
% %plot frame camera in line-scan coordinate frame
% scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
% %plot fov of line-scan on view-plane
% plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
% 
% 
% 
% figure(figGP_ZeroMeanSYM);
% 
% %Plotting training points on figure
% targetPntLSHom = (poseLS\targetPntFrameHomSYM)';
% 
% lsTrainPlotPt = (downsample(targetPntLSHom, downSamplingGP))';
% scatter3(lsTrainPlotPt(2,:), lsTrainPlotPt(3,:), (maxRadiInt+1)*ones(size(lsTrainPlotPt(3,:))), 20, [0,1,1], 'filled','MarkerEdgeColor', [0,0,0], 'LineWidth', 0.1);
% 
% radIntGP = reshape(mu,rows);
% radVar = reshape(varMu,rows);
% 
% %update figure with GP model
% surfGP.ZData = radIntGP;
% surfGP_var.ZData = radVar;
% drawnow();
% 
% 
% %plot 2D view-plane surface with line-scan camera as origin. Only plotting
% %YZ plane.
% figDiffGP_ZeroMeanSYM = figure('Name', 'Radiant intensity view-plane from GP');
% % s1 = subplot(1,3,1);
% % surf(Y, Z, radIntMagGroundtruth, 'EdgeColor', 'none'); hold on;
% % xlabel('y');
% % ylabel('z');
% % title('Ground-truth')
% % view(2);
% % scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
% % colormap(s1, hot);
% % caxis([0, maxRadiInt]);
% % colorbar;
% % axis equal;
% 
% %plot light source in the line-scan camera coordinate frame (green)
% scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
% %plot frame camera in line-scan coordinate frame
% scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
% %plot fov of line-scan on view-plane
% plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
% 
% 
% s2 = subplot(1,2,1);
% surf(Y, Z, radIntGP, 'EdgeColor', 'none'); hold on;
% xlabel('y');
% ylabel('z');
% title('Zero-mean GP Non-Symmetric')
% view(2);
% scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
% colormap(s2, hot);
% caxis([0, maxRadiInt]);
% colorbar;
% % axis equal;
% 
% %plot light source in the line-scan camera coordinate frame (green)
% scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
% %plot frame camera in line-scan coordinate frame
% scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
% %plot fov of line-scan on view-plane
% plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
% 
% 
% diffRad = abs(radIntGP - radIntMagGroundtruth);
% maxDiff = max(diffRad, [], 'all');
% 
% % Absolute difference plot
% s3 = subplot(1,2,2);
% surf(Y, Z, diffRad, 'EdgeColor', 'none'); hold on;
% xlabel('y');
% ylabel('z');
% title('Absolute Difference')
% view(2);
% scatter3(0, 0, (maxDiff+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
% colormap(s3, bone);
% caxis([0, maxDiff]);
% colorbar;
% % axis equal;
% 
% %plot light source in the line-scan camera coordinate frame (green)
% scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
% 
% %plot frame camera in line-scan coordinate frame
% scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
% 
% %plot fov of line-scan on view-plane
% plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
% drawnow();

%% Building model with GP non-Symmetric Constant Mean

cab(figSim, figViewPlane, figLeastSqrViewPlane, figGP_ZeroMean, ...
    figDiffGP_ZeroMean);

%Build zero-mean GP model
[mu, varMu, hypOpt] = LightSrcOptmGP(1, gpTrainingdata, testingX, downSamplingGP, 1000, intNoiseSigma, false);

figGP_ConstMean = figure('Name','GP Constant Mean');
s1 = subplot(1,3,1);
surfGP = surf(Y, Z, zeros(size(Y)), 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Const-mean GP View-Plane')
view(2);
scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s1, hot);
caxis([0, maxRadiInt]);
colorbar;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);

s2 = subplot(1,3,2);
surfGP_var = surf(Y, Z, zeros(size(Y)), 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Variance')
view(2);
scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s2, turbo);
caxis([0, 1]);
colorbar;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);



figure(figGP_ConstMean);

%Plotting training points on figure
targetPntLSHom = (T_LS_2_F\targetPntFrameHom)';

lsTrainPlotPt = (downsample(targetPntLSHom, downSamplingGP))';
scatter3(lsTrainPlotPt(2,:), lsTrainPlotPt(3,:), (maxRadiInt+1)*ones(size(lsTrainPlotPt(3,:))), 20, [0,1,1], 'filled','MarkerEdgeColor', [0,0,0], 'LineWidth', 0.1);

radIntGP = reshape(mu,rows);
radVar = reshape(varMu,rows);

%update figure with GP model
surfGP.ZData = radIntGP;
surfGP_var.ZData = radVar;

colormap(s2, turbo);
caxis([0, max(radVar, [], 'all')]);
drawnow();

diffRad = abs(radIntGP - radIntMagGroundtruth);
maxDiff = max(diffRad, [], 'all');

% Absolute difference plot
s3 = subplot(1,3,3);
surf(Y, Z, diffRad, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Absolute Difference')
view(2);
scatter3(0, 0, (maxDiff+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s3, bone);
caxis([0, maxDiff]);
colorbar;
% axis equal;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
drawnow();



%plot 2D view-plane surface with line-scan camera as origin. Only plotting
%YZ plane.
figDiffGP_ConstMean = figure('Name', 'Radiant intensity view-plane from GP');
% s1 = subplot(1,3,1);
% surf(Y, Z, radIntMagGroundtruth, 'EdgeColor', 'none'); hold on;
% xlabel('y');
% ylabel('z');
% title('Ground-truth')
% view(2);
% scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
% colormap(s1, hot);
% caxis([0, maxRadiInt]);
% colorbar;
% axis equal;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);


s2 = subplot(1,2,1);
surf(Y, Z, radIntGP, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Const-mean GP View-Plane')
view(2);
scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s2, hot);
caxis([0, maxRadiInt]);
colorbar;
% axis equal;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);


diffRad = abs(radIntGP - radIntMagGroundtruth);
maxDiff = max(diffRad, [], 'all');

% Absolute difference plot
s3 = subplot(1,2,2);
surf(Y, Z, diffRad, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Absolute Difference')
view(2);
scatter3(0, 0, (maxDiff+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s3, bone);
caxis([0, maxDiff]);
colorbar;
% axis equal;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
drawnow();

%% Building model with GP Non-Symmetric light source mean function

cab(figSim, figViewPlane, figLeastSqrViewPlane, figGP_ZeroMean, ...
    figDiffGP_ZeroMean, figGP_ConstMean, ...
    figDiffGP_ConstMean);


%Build zero-mean GP model
[mu, varMu, hypOpt] = LightSrcOptmGP(2, gpTrainingdata, testingX, downSamplingGP, 10000, intNoiseSigma, false);

figGP_LightSrcMean = figure('Name','GP Light Source Mean Non-Symmetric');
s1 = subplot(1,3,1);
surfGP = surf(Y, Z, zeros(size(Y)), 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Light-Source-mean GP View-Plane')
view(2);
scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s1, hot);
caxis([0, maxRadiInt]);
colorbar;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);

s2 = subplot(1,3,2);
surfGP_var = surf(Y, Z, zeros(size(Y)), 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Variance')
view(2);
scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s2, turbo);
caxis([0, 1]);
colorbar;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);



figure(figGP_LightSrcMean);

%Plotting training points on figure
targetPntLSHom = (T_LS_2_F\targetPntFrameHom)';

lsTrainPlotPt = (downsample(targetPntLSHom, downSamplingGP))';
scatter3(lsTrainPlotPt(2,:), lsTrainPlotPt(3,:), (maxRadiInt+1)*ones(size(lsTrainPlotPt(3,:))), 20, [0,1,1], 'filled','MarkerEdgeColor', [0,0,0], 'LineWidth', 0.1);

radIntGP = reshape(mu,rows);
radVar = reshape(varMu,rows);

%update figure with GP model
surfGP.ZData = radIntGP;
surfGP_var.ZData = radVar;

colormap(s2, turbo);
caxis([0, max(radVar, [], 'all')]);
drawnow();

diffRad = abs(radIntGP - radIntMagGroundtruth);
maxDiff = max(diffRad, [], 'all');

% Absolute difference plot
s3 = subplot(1,3,3);
surf(Y, Z, diffRad, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Absolute Difference')
view(2);
scatter3(0, 0, (maxDiff+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s3, bone);
caxis([0, maxDiff]);
colorbar;
% axis equal;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
drawnow();



%plot 2D view-plane surface with line-scan camera as origin. Only plotting
%YZ plane.
figDiffGP_LightSrcMean = figure('Name', 'Radiant intensity view-plane from GP');
% s1 = subplot(1,3,1);
% surf(Y, Z, radIntMagGroundtruth, 'EdgeColor', 'none'); hold on;
% xlabel('y');
% ylabel('z');
% title('Ground-truth')
% view(2);
% scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
% colormap(s1, hot);
% caxis([0, maxRadiInt]);
% colorbar;
% axis equal;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);


s2 = subplot(1,2,1);
surf(Y, Z, radIntGP, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Light-Source-mean View-Plane')
view(2);
scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s2, hot);
caxis([0, maxRadiInt]);
colorbar;
% axis equal;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);


diffRad = abs(radIntGP - radIntMagGroundtruth);
maxDiff = max(diffRad, [], 'all');

% Absolute difference plot
s3 = subplot(1,2,2);
surf(Y, Z, diffRad, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Absolute Difference')
view(2);
scatter3(0, 0, (maxDiff+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s3, bone);
caxis([0, maxDiff]);
colorbar;
% axis equal;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(T_S_2_LS(2,4), T_S_2_LS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(T_F_2_LS(2,4), T_F_2_LS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
drawnow();

%% 

figure('Name', 'ABS of LS vs GP');

diffRadLS = abs(radIntMagGroundtruth - radIntMagLeastSqr);
diffRadGP = abs(radIntMagGroundtruth - radIntGP);
diffRadZero = abs(radIntMagGroundtruth - radIntGPZero);


maxDiff = max(diffRadLS, [], 'all');
maxDiff = max([maxDiff, max(diffRadGP, [], 'all')], [], 'all');


% Absolute difference plot
s1 = subplot(1,2,1);
surf(Y, Z, diffRadLS, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Absolute Difference in Least Squares')
view(2);
scatter3(0, 0, (maxDiff+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s1, bone);
caxis([0, maxDiff]);
colorbar;



% Absolute difference plot
s2 = subplot(1,2,2);
surf(Y, Z, diffRadGP, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Absolute Difference in GP')
view(2);
scatter3(0, 0, (maxDiff+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s2, bone);
caxis([0, maxDiff]);
colorbar;


Yvec = Y(:);
Zvec = Z(:);
diffLSvec = diffRadLS(:);
diffGPvec = diffRadGP(:);
diffGPZvec = diffRadZero(:);

YpolyVP = viewPlanePoly.Vertices(:,1);
ZpolyVP = viewPlanePoly.Vertices(:,2);

ptsInVPMask = inpolygon(Yvec, Zvec, YpolyVP, ZpolyVP);

diffLS_in_VP = diffLSvec(ptsInVPMask);
% diffLS_in_VP(end+1:10000) = 0;

diffLS_out_VP = diffLSvec(~ptsInVPMask);
% diffLS_out_VP(end+1:10000) = 0;


diffGP_in_VP = diffGPvec(ptsInVPMask);
% diffGP_in_VP(end+1:10000) = 0;

diffGP_out_VP = diffGPvec(~ptsInVPMask);
% diffGP_out_VP(end+1:10000) = 0;

diffGPZ_out_VP = diffGPZvec(~ptsInVPMask);
diffGPZ_in_VP = diffGPZvec(ptsInVPMask);



% radDiffCell = {diffLS_in_VP, diffLS_out_VP, diffGP_in_VP, diffGP_out_VP};
% varLabel = {'LS in VP', 'LS out VP', 'GP in VP', 'GP out VP'};


radDiffCell = {diffLS_in_VP, diffLS_out_VP, diffGP_in_VP, diffGP_out_VP, diffGPZ_in_VP, diffGPZ_out_VP};
varLabel = {'LS in VP', 'LS out VP', 'GP in VP', 'GP out VP', 'GP ZERO in VP', 'GP ZERO out VP'};
% 
figure();
violin(radDiffCell, 'xlabel', varLabel);


% radDiffTable = table(diffLS_in_VP, diffLS_out_VP, diffGP_in_VP, diffGP_out_VP, 'VariableNames', {'LS in VP', 'LS out VP', 'GP in VP', 'GP out VP'});
% 
% 
% figure();
% 
% vioPl = violinplot(radDiffTable);
% vioPl.Sh
% viol1 = violinplot(diffLS_in_VP, 'LS in VP'); hold on;
% viol2 = violinplot(diffLS_out_VP, 'LS out VP'); 



grid on;

%%


function vp_callback(src, ~, vpPatch)
value = get(src,'Value');

if value
    vpPatch.Visible = true;
else
    vpPatch.Visible = false;
end

end