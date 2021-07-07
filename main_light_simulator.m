%Initial script for simulating a non-isotropic point light source.

close all;
clear;

%robotics toolbox
run(['rvctools' filesep 'startup_rvc.m'])

%code for this project
addpath('src_code');

%GPML toolbox
run(['gpml-matlab-master', filesep, 'startup.m']);


%% Parameters to change

%variance of pixel intensity noise in line-scan camera measurements
intNoiseSigma = 0.1;


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

RIDSplineFit = csapi(thetaRID,radiantIntenRID);
thetas = linspace(-pi/2, pi/2, 1000);

radiantIntSpline = fnval(RIDSplineFit, thetas);
polarplot(thetas,radiantIntSpline, '-b');

%% Setting up non-isotropic light source

%Pose of source in world coordinates (frame camera)
% locLightSrc = [0.1;0.1;-0.1];
% rotLightSrc = eul2rotm(deg2rad([0,-10, 0]), 'ZYX');

locLightSrc = [0.18; -0.05; 0.05];
rotLightSrc = [0.792897128208011,0.00375573235981344,-0.609343941098892;0.00375573235981344,0.999931891212147,0.0110502222303317;0.609343941098892,-0.0110502222303317,0.792829019420159];
poseLightSrc = [rotLightSrc, locLightSrc; 0, 0, 0, 1];

%radiant intensity distribution parameters
maxRadiInt = 10;
mu = 2;

%attentuation
r_d = 1; %effective radius of disk source

%Starting distance from the source in metres for visualisation
distFromSrc = 0.2;

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
lightSrc = LightSimulator(locLightSrc, rotLightSrc, maxRadiInt, mu, r_d, distFromSrc, figSim, S, RIDSplineFit);

%% Setting up Line-scan camera

%assuming no distortion and no uncertainty in parameters

%line-scan intrinsic camera parameters
fyLS  = 800;
v0 = 160;
yPix = 320;

lsIntrinsic = cameraIntrinsics([1, fyLS], [realmin,v0],[yPix, 1]);

rotLS = eul2rotm(deg2rad([90, 0, -10]), 'ZYX');
rotLS = rotLightSrc;
tLS = [-0.1; 0.1; -0.1];
tLS = locLightSrc + [0; 0.1; 0];



% tLS = [0.125101976645449;-0.0169881345205247;-0.110138876985863];
% rotLS = [-0.970289293330054,-0.0609119653632981,-0.234154691869808;0.0563354557071021,-0.998068323042120,0.0261904366165156;-0.235297691614978,0.0122210889642028,0.971846511186408];

%pose of line-scan w.r.t frame camera
% rotLS = rotLightSrc;
% tLS = [0.125101976645449;-0.0169881345205247;-0.10138876985863];

poseLS = [rotLS, tLS; 0,0,0,1];
extLS = inv(poseLS);

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
xyzTrans = poseLS*homPts;

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
poseLightSrcLS = extLS * poseLightSrc;

%plot the light intensity produced by the light source on the line-scan
%camera's view-plane
yDist = 0.5; %max distance left/right
zDist = 1; %max distance along z-direction

%Create mesh
y = linspace(-yDist, yDist, 100);
z = linspace(poseLightSrcLS(3,4), zDist, 100);
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
ptsHomFrame = poseLS*ptsHom;

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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);

%% Reflective target
% assume that the target is an infinitely tall flat surface which is
% 10cm wide. The origin on the target will be in the middle.
%target will be positioned randonly in the line-scan camera view.
% Z-axis of the target is inline with its normal. i.e. z-axis vector = surface normal
target.Width = 0.1; %width of target in metres
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
targetSamplePosFrame = cell(1,noSamples);
targetSampleNormal = zeros(3, noSamples);

%for storing symmetric points
targetSamplePosSymmFrame = cell(1,noSamples);
targetSampleSymmNormal = zeros(3, noSamples);


locLightSrcLS = tform2trvec(poseLightSrcLS)';
rotLightSrcDirLS = tform2rotm(poseLightSrcLS)*[0;0;1];

minTheta = deg2rad(-10);
maxTheta = deg2rad(60);

for sample = 1:noSamples
    %get random distance between line-scan and target in metres
    distLS2Target = (maxRange - minRange)*rand() + minRange;
    %plot the light intensity produced by the light source on the line-scan
    
    %pixel location of the centre of target in line-scan image
    pixTargetCent = randi([0,yPix]);
    pixTargetCentHom = [0 ; pixTargetCent; 1];
    normTargetCent = KMat\pixTargetCentHom;
    normTargetCent = normTargetCent./normTargetCent(3);
    
    %find the 3D line equation between the optical centre and the point on
    %the normalized plane. equation of line is of the form r = r0 + t*r_dir
    %distance between the optical centre and target centre is taken to be
    %z-coordinate. Use that to calculate the t
    t = distLS2Target - 1;
    %coordinate of the target centre from the line-scan camera, i.e. translation
    targetCentPos = normTargetCent + t.*normTargetCent;
    
    %pick angle about x-axis between 90deg and 270 deg to get the pose of
    %the target. Only assume rotation about x-axis (yaw only)
    %                 a = deg2rad(100);
    %                 b = deg2rad(260);
    %                 xAngle = (b-a)*rand()+(a);
    %     %initially assume a constant normal angle
    xAngle = pi;
    rotTargetLS = eul2rotm([0, 0, xAngle], 'ZYX');
    
    %pose of target w.r.t to line-scan camera c.f
    poseTargetLS = [rotTargetLS, targetCentPos; 0, 0, 0, 1];
    extTargetLS = poseTargetLS / eye(4);
    
    %pose of target w.r.t frame camera c.f
    poseTargetFrame =  poseLS*poseTargetLS;
    
    %normal vector of the target surface (z-axis)
    normTargetLS = poseTargetLS(1:3,3);
    targetSampleNormal(:,sample) = poseTargetFrame(1:3,3);
    
    %left and right edge points of the target
    targetEdgesPts = [[target.LeftEdge; 1], [target.RightEdge; 1]];
    
    %project target edges to line-scan image
    imgTargetEdges = projectPoints(targetEdgesPts', KMat, extTargetLS, [], lsIntrinsic.ImageSize);
    
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
    normTargetPixHom = normTargetPixHom./normTargetPixHom(3);
    
    %reproject these pixels to find their 3D location relative to the
    %line-scan camera
    targetPtsLS = zeros(3, length(vtargetPix));
    
    %find the origin of the symmetric target about the light source direction vector
    targetSymmPtsLS = zeros(3, length(vtargetPix));
    targetCentSymmPos = SymmetricPoint_LightSrc_Linescan(targetCentPos, locLightSrcLS, rotLightSrcDirLS, vpNormal);
    
    
    for pixel = 1:length(vtargetPix)
        pnt = normTargetPixHom(:,pixel);
        
        dirPnt = pnt./norm(pnt);
        
        %calculate 3D point of pixel on target relative to line-scan camera
        %using line-plane intersection
        [pntOnTarget, rc] = line_plane_intersection(dirPnt, pnt, normTargetLS, targetCentPos);
        
        %the point should be on the plane
        if rc ~= 1
            error('line-plane intersection was not on plane');
        end
        
        targetPtsLS(:, pixel) = pntOnTarget;
        
        %find symmetric point about light source direction on view-plane
        pntOnTargetSymm = SymmetricPoint_LightSrc_Linescan(pntOnTarget, locLightSrcLS, rotLightSrcDirLS, vpNormal);
        targetSymmPtsLS(:, pixel) = pntOnTargetSymm;
        
    end
    
    %calculate the symmetric normal, which is actually finding the rotation
    %of the symmetrical target and taking its z-directional axis
    %y direction vector found using edge points
    ydir = targetSymmPtsLS(:,end) - targetSymmPtsLS(:,1);
    ydir = ydir./norm(ydir);
    
    %x direction vector is same as the original
    xdir = poseTargetLS(1:3,1);
    zdir = cross(xdir, ydir);
    
    rot = [xdir, ydir, zdir];
    targetPoseSymLS = [rot, targetCentSymmPos; 0, 0, 0, 1];
    targetSymmPoseFrame = poseLS*targetPoseSymLS;
    
    targetSymmPtsFrame = poseLS*[targetSymmPtsLS; ones(1,length(targetSymmPtsLS(1,:)))];
    targetSymmPtsFrame = targetSymmPtsFrame(1:3,:);
    targetSampleSymmNormal(:,sample) = targetSymmPoseFrame(1:3,3);
    
    %transform points from line-scan to frame camera
    targetPtsFrame = poseLS*[targetPtsLS; ones(1,length(targetPtsLS(1,:)))];
    targetPtsFrame = targetPtsFrame(1:3,:);
    
    
    %plot line and axes of target in the simulator figure
    figure(figSim);
    line(targetPtsFrame(1,:), targetPtsFrame(2,:), targetPtsFrame(3,:), 'Color', [0,0,0], 'LineWidth', 2); hold on;
    trplot(poseTargetFrame, 'rviz', 'length', (maxRange - minRange)*0.1);
    
    %plot line and axes of symmetric target in the simulator figure
    line(targetSymmPtsFrame(1,:), targetSymmPtsFrame(2,:), targetSymmPtsFrame(3,:), 'Color', [0,0,0], 'LineWidth', 2);
    trplot(targetSymmPoseFrame, 'rviz', 'length', (maxRange - minRange)*0.1);
    
    
    drawnow();
    
    %plot line in the view-plane intensity plot
    figure(figViewPlane);
    line(targetPtsLS(2,:), targetPtsLS(3,:), maxRadiInt.*ones(size(targetPtsLS(2,:))), 'Color', [0,0,0], 'LineWidth', 2);
    
    %plot symmetric line in the view-plane intensity plot
%     line(targetSymmPtsLS(2,:), targetSymmPtsLS(3,:), maxRadiInt.*ones(size(targetSymmPtsLS(2,:))), 'Color', [0,0,0], 'LineWidth', 2);
    
    drawnow();
    
    %save the target pose in the c.f of the frame camera
    targetSamplePosFrame(sample) = {targetPtsFrame};
    
    %save symmetric target pose in the c.f of the frame camera
    targetSamplePosSymmFrame(sample) = {targetSymmPtsFrame};
end

%% Measure pixel intensity at the 3D location on the target which is relative to the frame camera

%stores the intensity measured at a specific 3D location. The measured
%intensity at the symmetric target would be the exact same
targetSampleInten = cell(1, noSamples);

for sample = 1:noSamples
    targetPtsFrame = targetSamplePosFrame{sample};
    targetNormal = targetSampleNormal(:, sample);
    
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
    targetSampleInten(sample) = {targetIntenNoise};
end

%% Combine the samples such that they can be analysed by regressional approaches



%number of trials is equal to the number of samples. Each consective trial
%will add an extra sample

% {noSamples} x {targetL, targetPntNormals, targetPntFrame, targetPntNormalsSymm, targetPntFrameSymm}
targetTrialsWithSymm = cell(noSamples, 5);

for trials = 1:noSamples
    targetL = zeros(1,2*trials*yPix);
    targetPntNormals = zeros(3,2*trials*yPix);
    targetPntFrame = zeros(3, 2*trials*yPix);
    targetPntNormalsSymm = zeros(3,2*trials*yPix);
    targetPntFrameSymm = zeros(3, 2*trials*yPix);
    
    noPts = 0;
    
    %combines all samples into a single array for the current trials
    for sample = 1:trials
        noCurPts = length(targetSampleInten{sample});
        
        targetL(noPts+1:noPts+noCurPts) = targetSampleInten{sample};
        targetPntNormals(:, noPts+1:noPts+noCurPts) = repmat(targetSampleNormal(:, sample), [1, noCurPts]);
        targetPntFrame(:, noPts+1:noPts+noCurPts) = targetSamplePosFrame{sample};
        targetPntNormalsSymm(:, noPts+1:noPts+noCurPts) = repmat(targetSampleSymmNormal(:, sample), [1, noCurPts]);
        targetPntFrameSymm(:, noPts+1:noPts+noCurPts) = targetSamplePosSymmFrame{sample};
        
        noPts = noPts + noCurPts;
    end
    
    %clip the arrays to correct size
    targetL = targetL(1:noPts);
    targetPntNormals = targetPntNormals(:, 1:noPts);
    targetPntFrame = targetPntFrame(:, 1:noPts);
    targetPntNormalsSymm = targetPntNormalsSymm(:,1:noPts);
    targetPntFrameSymm = targetPntFrameSymm(:,1:noPts);
    
    targetTrialsWithSymm(trials,:) = {targetL, targetPntNormals, targetPntFrame, targetPntNormalsSymm, targetPntFrameSymm};    
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
    
    targetL = targetTrialsWithSymm{trials,1};
    targetPntFrame = targetTrialsWithSymm{trials,3};
    targetPntNormals = targetTrialsWithSymm{trials,2};
    
    [optPhi,resNorm, curRankFIM, curEigFIM] = LightSrcOptmLS(lightSrc, Phi0, targetL, targetPntFrame, targetPntNormals, target.Reflectance, optOptions, intNoiseSigma);
    
    %store optimised parameters and residual
    optPhiTrials(trials, :) = optPhi;
    resNormTrials(trials) = resNorm;
    

    
    rankFIM(trials) = curRankFIM;
    eigFIM(trials,:) = curEigFIM;
end

%% Results of least squares

%close all figures except for the simulator and the view-plane ground truth
cab(figSim, figViewPlane);

%calculate absolute errors in the estimated parameters
abs_phi0 = abs(optPhiTrials(:,1) - maxRadiInt);
abs_mu = abs(optPhiTrials(:,2) - mu);
abs_r_d = abs(optPhiTrials(:,3) - r_d);

%plot absolute error
figAbsLS = figure('Name', 'Absolute Erorr in Estimated Parameters');
subplot(1,3,1);
plot(1:noSamples, abs_phi0);
xlabel('samples');
ylabel('absolute error in \Phi_0');
ylim([0, max(abs_phi0)])
grid on;

subplot(1,3,2);
plot(1:noSamples, abs_mu);
xlabel('samples');
ylabel('absolute error in \mu');
ylim([0, max(abs_mu)])
grid on;

subplot(1,3,3);
plot(1:noSamples, abs_r_d);
xlabel('samples');
ylabel('absolute error in r_d');
ylim([0, max(abs_r_d)])
grid on;

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
grid on;cab(figSim, figViewPlane);


%Compare Radiant intensity view-plane image of ground-truth parameters and
%estimated parameters

lightSrcLeastSqrs = LightSimulator(locLightSrc, rotLightSrc, optPhiTrials(end, 1), optPhiTrials(end, 2), optPhiTrials(end, 3));

%plot the light intensity produced by the light source on the line-scan
%camera's view-plane
yDist = 0.5; %max distance left/right
zDist = 1; %max distance along z-direction

%Create mesh
y = linspace(-yDist, yDist, 100);
z = linspace(poseLightSrcLS(3,4), zDist, 100);
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
ptsHomFrame = poseLS*ptsHom;

rows = size(X);
%reshape to mesh
XFrame = reshape(ptsHomFrame(1,:),rows);
YFrame = reshape(ptsHomFrame(2,:),rows);
ZFrame = reshape(ptsHomFrame(3,:),rows);

radIntMagGroundtruth = lightSrc.RadiantIntensityMesh(XFrame, YFrame, ZFrame);
radIntMagLeastSqr = lightSrcLeastSqrs.RadiantIntensityMesh(XFrame, YFrame, ZFrame);
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
% scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);


%% GP Testing Points

%transform points from frame camera C.F to light source C.F
ptsHomLightSrc = poseLightSrc\ptsHomFrame;
ptsLightSrc = ptsHomLightSrc(1:3, :);

[ptsRadius,ptsTheta] = cart2rtlightsrc(ptsLightSrc);

testingX = [ptsRadius', ptsTheta'];

%% GP Training Points
    
targetL = targetTrialsWithSymm{end,1};
targetPntFrame = targetTrialsWithSymm{end,3};
targetPntNormals = targetTrialsWithSymm{end,2};

targetL_SYM = [targetL,targetL];
targetPntFrame_SYM = [targetTrialsWithSymm{end,3}, targetTrialsWithSymm{end,5}];
targetPntNormals_SYM = [targetTrialsWithSymm{end,2}, targetTrialsWithSymm{end,4}];

%calculate the direction light vector (point to light source)
pntLightVec = locLightSrc - targetPntFrame;
dirPntLightVec = pntLightVec./vecnorm(pntLightVec);

%Calculate radiant intensity magnitude used for building model
radIntMagPnt = (targetL.*pi)./(target.Reflectance.*dot(targetPntNormals, dirPntLightVec,1));

%find the point w.r.t to light source c.f
targetPntFrameHom = [targetPntFrame; ones(1, length(targetPntFrame(1,:)))];
targetPntLightSrcHom = poseLightSrc\targetPntFrameHom;
targetPntLightSrc = targetPntLightSrcHom(1:3, :);

[ptsRadius,ptsTheta] = cart2rtlightsrc(targetPntLightSrc);
gpTrainingdata = [ptsRadius', ptsTheta', radIntMagPnt'];



%calculate the direction light vector (point to light source)
pntLightVec = locLightSrc - targetPntFrame_SYM;
dirPntLightVec = pntLightVec./vecnorm(pntLightVec);

%Calculate radiant intensity magnitude used for building model
radIntMagPnt = (targetL_SYM.*pi)./(target.Reflectance.*dot(targetPntFrame_SYM, dirPntLightVec,1));

%find the point w.r.t to light source c.f
targetPntFrameHomSYM = [targetPntFrame_SYM; ones(1, length(targetPntFrame_SYM(1,:)))];
targetPntLightSrcHom = poseLightSrc\targetPntFrameHomSYM;
targetPntLightSrc = targetPntLightSrcHom(1:3, :);

[ptsRadius,ptsTheta] = cart2rtlightsrc(targetPntLightSrc);
gpTrainingdataSYM = [ptsRadius', ptsTheta', radIntMagPnt'];

downSamplingGP = 50;

%% Building model with GP Zero-mean non-symmetric

cab(figSim, figViewPlane, figAbsLS, figLeastSqrViewPlane);

plotMaxHeight = 100;

%Build zero-mean GP model
[mu, varMu, hypOpt] = LightSrcOptmGP(0, gpTrainingdata, testingX, downSamplingGP, 1000, false);

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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (plotMaxHeight+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (plotMaxHeight+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (plotMaxHeight+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (plotMaxHeight+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (plotMaxHeight+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);



figure(figGP_ZeroMean);

%Plotting training points on figure
targetPntLSHom = (poseLS\targetPntFrameHom)';

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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

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
% scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
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
% scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
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
% scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
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
% scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
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
% scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
% 
% %plot frame camera in line-scan coordinate frame
% scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
% 
% %plot fov of line-scan on view-plane
% plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
% drawnow();

%% Building model with GP non-Symmetric Constant Mean

cab(figSim, figViewPlane, figAbsLS, figLeastSqrViewPlane, figGP_ZeroMean, ...
    figDiffGP_ZeroMean);

%Build zero-mean GP model
[mu, varMu, hypOpt] = LightSrcOptmGP(1, gpTrainingdata, testingX, downSamplingGP, 1000, false);

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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);



figure(figGP_ConstMean);

%Plotting training points on figure
targetPntLSHom = (poseLS\targetPntFrameHom)';

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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
drawnow();

%% Building model with GP Non-Symmetric light source mean function

cab(figSim, figViewPlane, figAbsLS, figLeastSqrViewPlane, figGP_ZeroMean, ...
    figDiffGP_ZeroMean, figGP_ConstMean, ...
    figDiffGP_ConstMean);


%Build zero-mean GP model
[mu, varMu, hypOpt] = LightSrcOptmGP(2, gpTrainingdata, testingX, downSamplingGP, 10000, false);

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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);



figure(figGP_LightSrcMean);

%Plotting training points on figure
targetPntLSHom = (poseLS\targetPntFrameHom)';

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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
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
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
drawnow();

%% 

figure()

diffRad = abs(radIntMagGroundtruth - radIntMagLeastSqr);


% Absolute difference plot
s1 = subplot(1,2,1);
surf(Y, Z, diffRad, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Absolute Difference in Least Squares')
view(2);
scatter3(0, 0, (maxDiff+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s1, bone);
caxis([0, maxDiff]);
colorbar;

diffRad = abs(radIntGP - radIntMagGroundtruth);


% Absolute difference plot
s2 = subplot(1,2,2);
surf(Y, Z, diffRad, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Absolute Difference in GP')
view(2);
scatter3(0, 0, (maxDiff+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s2, bone);
caxis([0, maxDiff]);
colorbar;
%%


function vp_callback(src, ~, vpPatch)
value = get(src,'Value');

if value
    vpPatch.Visible = true;
else
    vpPatch.Visible = false;
end

end