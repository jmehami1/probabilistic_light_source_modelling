%Initial script for simulating a non-isotropic point light source.

close all;
clear;

%robotics toolbox
run(['rvctools' filesep 'startup_rvc.m'])

%code for 3D line plane intersection
addpath(genpath('line_plane_intersection'));

%code for this project
addpath('code');


%use the CMAES Optimisation instead of non-linear least squares
useCMAES = false;

%variance of pixel intensity noise in line-scan camera measurements
intNoiseVar = 0.01;

%% Setting up non-isotropic light source

%Pose of source in world coordinates (frame camera)
ligSourLoc = [0.1;0.1;-0.1];
ligSourOrien = eul2rotm(deg2rad([0,-10, 0]), 'ZYX');
lightSrcPoseFrame = [ligSourOrien, ligSourLoc; 0, 0, 0, 1];

%radiant intensity distribution parameters
maxRadiInt = 10;
mu = 1.5;

%attentuation
rAtt = 1; %effective radius of disk source

%light source object
lightSrc = LightSimulator(ligSourLoc, ligSourOrien, maxRadiInt, mu, rAtt);

%% Setting up light simulator figure

%figure for plotting
fig = figure('Name', 'Light Simulator with Frame Camera Origin');
plotCamera('Location', zeros(1,3), 'Orientation', eye(3), 'Size', 0.05, 'AxesVisible', true); hold on;
S = surf(eye(2), 'Visible', 'off');
xlabel('x');
ylabel('y');
zlabel('z');
grid on;
axis equal;

distFromSrc = 0.2; %Starting distance from the source in metres

%all for distFromSrc slider
b = uicontrol('Parent',fig,'Style','slider','Position',[81,100,419,23],...
    'value',distFromSrc, 'min',0, 'max',1);
bgcolor = fig.Color;
bl1 = uicontrol('Parent',fig,'Style','text','Position',[50,70,23,23],...
    'String','0','BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',fig,'Style','text','Position',[500,70,23,23],...
    'String','1','BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',fig,'Style','text','Position',[240,70,200,23],...
    'String','Distance from the source','BackgroundColor',bgcolor);
%callback function at the end of the script
b.Callback = @(src, eventData) distFromSrc_callback(src, eventData, lightSrc, fig, S);

lightSrc.PlotLightSource(fig);

%% Setting up Line-scan camera

%assuming no distortion and no uncertainty in parameters

%line-scan intrinsic camera parameters
fyLS  = 800;
v0 = 160;
yPix = 320;

lsIntrinsic = cameraIntrinsics([1, fyLS], [realmin,v0],[yPix, 1]);

%combination of the ceoefficient of clarity and optical coefficient
%betaLS = 0.95;targetNormal

%pose of line-scan w.r.t frame camera
rotLS = eul2rotm(deg2rad([90, 0, -10]), 'ZYX');
tLS = [-0.1, 0.1, -0.1]';
plotCamera('Location', tLS, 'Orientation', rotLS, 'Size', 0.05, 'AxesVisible', true, 'Color', [0,0,1]); hold on;

%calculate line-scan view-plane FOV
fovLS =  atan((yPix/2)/fyLS)*2;

%minimum/maximum working focal range of the line-scan camera in metres
minRange = 0.3;
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

%transformation from line-scan to frame camera coordinate frame
T_LS2Frame = [rotLS', tLS; 0, 0, 0, 1];

%transform verticies of view-plane from frame to line-scan coordinate frame
homPts = [xPoly;yPoly;zPoly; ones(size(zPoly))];
xyzTrans = T_LS2Frame*homPts;

%plot plane
x = xyzTrans(1,:);
y = xyzTrans(2,:);
z = xyzTrans(3,:);
patch(x,y,z, [0,0,1], 'FaceAlpha', 0.3);

%% 2D plot of radiant intensity on view-plane

%plot it as a rectangle starting at where the line-scan is positioned and
%only plot the YZ-plane
%size of plane in each dimension
yDist = 0.5;
zDist = 1;

%Create mesh
y = linspace(-yDist, yDist, 100);
z = linspace(0, zDist, 100);
x = 0;

[X,Y,Z] = meshgrid(x,y,z);

%remove extra necessary singular dimension
X = squeeze(X);
Y = squeeze(Y);
Z = squeeze(Z);

%These points are in the line-scan coordinate frame, transform them into
%the frame camera coordinate frame (world coordinates
%rotation
XYZrot = [X(:),Y(:),Z(:)]*rotLS;

rows = size(X);
%reshape to mesh and translate
XFrame = reshape(XYZrot(:,1),rows) + tLS(1);
YFrame = reshape(XYZrot(:,2),rows) + tLS(2);
ZFrame = reshape(XYZrot(:,3),rows) + tLS(3);

radIntMag = lightSrc.RadiantIntensityMesh(XFrame, YFrame, ZFrame);

%plot 2D view-plane surface with line-scan camera as origin. Only plotting
%YZ plane.
figGroundTruth = figure('Name', 'Radiant intensity on line-scan view-plane');
surf(Y, Z, radIntMag, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('View-plane RADIANT INTENSITY')
view(2);
scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled'); %line-scan origin
colormap(jet(100));
colorbar; caxis([0, maxRadiInt]);
axis equal;

%plot green light source
lightSrcLS = T_LS2Frame \ [ligSourLoc; 1]; %transform to line-scan coordinate system
scatter3(lightSrcLS(2), lightSrcLS(3), (maxRadiInt+1), 200, [0,1,0], 'filled');

%plot frame camera
frameCameraLS = inv(T_LS2Frame);
scatter3(frameCameraLS(2,4), frameCameraLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled');


%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);


%% Reflective target
% assume that the target is an infinitely tall flat surface which is
% 10cm wide. The origin on the target will be in the middle.
%target will be positioned randonly in the line-scan camera view.
% Z-axis of the target is inline with its normal. i.e. z-axis vector = surface normal
target.Width = 0.1; %width of target in metres
target.Reflectance = 0.95; %reflectance of target for any wavelength
target.LeftEdge = [0; -target.Width/2; 0];
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
noSamples = 100;

KMat = lsIntrinsic.IntrinsicMatrix';
KMat(1,3) = 0;

targetSamplePosFrame = cell(1,noSamples);
targetSamplePolarLS = cell(1,noSamples);
targetSampleNormal = zeros(3, noSamples);

lightSrcPoseLS = T_LS2Frame\lightSrcPoseFrame;

for sample = 1:noSamples
    %get random distance between line-scan and target in metres
    distLS2Target = (maxRange - minRange)*rand() + minRange;
    
    %pixel location of the centre of target in line-scan image
    pixTargetCent = randi([1,yPix]);
    pixTargetCentHom = [0 ; pixTargetCent; 1];
    normTargetCent = KMat\pixTargetCentHom;
    normTargetCent = normTargetCent./normTargetCent(3);
    
    %find the 3D line equation between the optical centre and the point on
    %the normalized plane. equation of line is of the form r = r0 + t*r_dir
    %distance between the optical centre and target centre is taken to be
    %z-coordinate. Use that to calculate the t
    t = distLS2Target - 1;
    %coordinate of the target centre from optical centre, i.e. translation
    targetCentPos = normTargetCent + t.*normTargetCent;
    
    %pick angle about x-axis between 90deg and 270 deg to get the pose of
    %the board. Only assume rotation about x-axis (yaw only)
    %         a = deg2rad(120);
    %         b = deg2rad(250);
    %         xAngle = (b-a)*rand()+(a);
    %     %initially assume a constant normal angle
    xAngle = pi;
    rotTarget = eul2rotm([0, 0, xAngle], 'ZYX');
    
    
    %     rotTarget = eul2rotm([0,0,pi], 'ZYX')*tform2rotm(lightSrcPoseLS);
    %     rotTarget = tform2rotm(T_LS2Frame\rotm2tform(rotTarget));
    
    %Get normal of plane from the rotation
    normTarget = rotTarget(1:3,3);
    targetSampleNormal(:,sample) = normTarget;
    
    %pose of target relative to line-scan camera
    targetPosLS = [rotTarget, targetCentPos; 0, 0, 0, 1];
    
    %Pose of target relative to frame camera
    T_target2Frame = T_LS2Frame*targetPosLS;
    
    %transform line edges to the frame camera coordinate system
    targetEdge = [[target.LeftEdge; 1], [target.RightEdge; 1]];
    targetEdgeLS = targetPosLS*targetEdge;
    targetEdgeFrame = T_target2Frame*targetEdge;
    
    %plot line in the simulator figure
    figure(fig);
    line(targetEdgeFrame(1,:), targetEdgeFrame(2,:), targetEdgeFrame(3,:), 'Color', [0,0,0], 'LineWidth', 2);
    %     scatter3(targetEdgeFrame(1,1), targetEdgeFrame(2,1), targetEdgeFrame(3,1), 50, [1,1,0], 'filled');
    %     scatter3(targetEdgeFrame(1,2), targetEdgeFrame(2,2), targetEdgeFrame(3,2), 50, [1,0,1], 'filled');
    text(T_target2Frame(1,4) + 0.05*T_target2Frame(1,4), T_target2Frame(2,4) + 0.05*T_target2Frame(2,4), ...
        T_target2Frame(3,4) + 0.05*T_target2Frame(3,4), num2str(sample));
    %plot coordinate frame of target w.r.t frame camera (world coordinates)
    trplot(T_target2Frame, 'rviz', 'length', (maxRange - minRange)*0.1);
    drawnow();
    
    %plot on view-plane
    %     figure(figGroundTruth);
    %     line(targetEdgeLS(2,:), targetEdgeLS(3,:), (maxRadiInt+1)*ones(size(targetEdgeFrame(2,:))), 'Color', [0,0,0], 'LineWidth', 2);
    %     text(targetPosLS(2,4) + 0.05*targetPosLS(2,4), targetPosLS(3,4) + 0.05*targetPosLS(3,4), ...
    %         (maxRadiInt+1), num2str(sample));
    %
    %project target edge points to line-scan image
    imgCoorTarget = projectPoints(targetEdge', KMat, targetPosLS, [], lsIntrinsic.ImageSize);
    
    %extract v coordinate of pixels
    vImgCoorTarget = imgCoorTarget(:,2);
    
    %Check if the edge pixels are within the size of the line-scan image, else
    %set the edge pixels to be the limits of the line image
    for j = 1:2
        
        vImgCoorTarget(j) = round(vImgCoorTarget(j));
        
        if vImgCoorTarget(j) < 1
            vImgCoorTarget(j) = 1;
        elseif vImgCoorTarget(j) > yPix
            vImgCoorTarget(j) = yPix;
        end
    end
    
    %vector of all pixels that see the target
    if vImgCoorTarget(1) > vImgCoorTarget(2)
        vtargetPix = vImgCoorTarget(2):vImgCoorTarget(1);
    else
        vtargetPix = vImgCoorTarget(1):vImgCoorTarget(2);
    end
    
    %transform to normalized coordinates
    targetPixHom  = [zeros(size(vtargetPix)); vtargetPix; ones(size(vtargetPix))];
    normTargetPixHom = KMat\targetPixHom;
    normTargetPixHom = normTargetPixHom./normTargetPixHom(3);
    
    %reproject these pixels to find their 3D location relative to the
    %line-scan camera
    targetPtsFrame = zeros(3, length(vtargetPix));
    targetPolarLS = zeros(2, length(vtargetPix));
    
    
    for pixel = 1:length(vtargetPix)
        pnt = normTargetPixHom(:,pixel);
        
        %calculate 3D point of pixel on target relative to line-scan camera
        %using line-plane intersection
        [pntOnTarget, rc] = line_plane_intersection(pnt, pnt, normTarget, targetCentPos);
        
        %the point should be on the plane
        if rc ~= 1
            error('line-plane intersection was not on plane');
        end
        
        %transform points from cartesian coordinates to polar
        [theta, rho] = cart2pol(pntOnTarget(2), pntOnTarget(3));
        targetPolarLS(:, pixel) = [theta - pi/2; rho];
        
        %transform target points to the frame camera coordinate system
        pntTargetFrame = T_LS2Frame*[pntOnTarget; 1];
        targetPtsFrame(:,pixel) = pntTargetFrame(1:3);
    end
    
    
    targetSamplePolarLS(sample) = {targetPolarLS};
    targetSamplePosFrame(sample) = {targetPtsFrame};
end

%% 2D plot of radiant intensity on view-plane dot product with normal

%plot it as a rectangle starting at where the line-scan is positioned and
%only plot the YZ-plane
%size of plane in each dimension
% yDist = 0.5;
% zDist = 1;

%Create mesh
y = linspace(-yDist, yDist, 100);
z = linspace(0, zDist, 100);
x = 0;

[X,Y,Z] = meshgrid(x,y,z);

%remove extra necessary singular dimension
X = squeeze(X);
Y = squeeze(Y);
Z = squeeze(Z);

%These points are in the line-scan coordinate frame, transform them into
%the frame camera coordinate frame (world coordinates
%rotation
XYZrot = [X(:),Y(:),Z(:)]*rotLS;

rows = size(X);
%reshape to mesh and translate
XFrame = reshape(XYZrot(:,1),rows) + tLS(1);
YFrame = reshape(XYZrot(:,2),rows) + tLS(2);
ZFrame = reshape(XYZrot(:,3),rows) + tLS(3);

% [~, radIntMagVec] = lightSrc.RadiantIntensityMesh(XFrame, YFrame, ZFrame);

radianceInTarget = zeros(size(XFrame));

for i = 1:size(XFrame,1)
    for j = 1:size(XFrame, 2)
        pnt = [XFrame(i,j), YFrame(i,j), ZFrame(i,j)];
        
        radianceInTarget(i, j) = lightSrc.RadianceInMaterialPoint(pnt, normTarget) + ...
            normrnd(0, intNoiseVar);
    end
end

%plot 2D view-plane surface with line-scan camera as origin. Only plotting
%YZ plane.
figGroundTruthDotTarget = figure('Name', 'Radiant intensity on line-scan view-plane');
surf(Y, Z, radianceInTarget, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title(['View-plane RADIANCE hitting target, noise of VAR ', num2str(intNoiseVar)]);
view(2);
scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled'); %line-scan origin
colormap(jet(100));
colorbar; caxis([0, maxRadiInt]);
axis equal;

%plot green light source
lightSrcLS = T_LS2Frame \ [ligSourLoc; 1]; %transform to line-scan coordinate system
scatter3(lightSrcLS(2), lightSrcLS(3), (maxRadiInt+1), 200, [0,1,0], 'filled');

%plot frame camera
frameCameraLS = inv(T_LS2Frame);
scatter3(frameCameraLS(2,4), frameCameraLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled');

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);

%% Measure pixel intensity at the 3D location on the target which is relative to the frame camera

targetSampleInten = cell(1, noSamples);


for sample = 1:noSamples
    targetPtsFrame = targetSamplePosFrame{sample};
    targetNormal = targetSampleNormal(:, sample);
    targetPolarLS = targetSamplePolarLS{sample};
    
    %intensity measured at the 3D location by the line-scan camera in the
    %coordinate frame of the frame camera.
    %***** NOTE: Darkening effects, vignetting, sensor irregularties have
    %not been considered yet
    
    numPts = size(targetPtsFrame, 2);
    
    targetInten = zeros(1,numPts);
    
    for pt = 1:numPts
        %         angleLightRayLS = targetPolarLS(1,pt);
        %         vignet = betaLS*cos(angleLightRayLS)^4;
        vignet = 1;
        targetInten(pt) = vignet*lightSrc.RadianceOutMaterialPoint(targetPtsFrame(:,pt), targetNormal, target.Reflectance);
    end
    
    %add gaussian white noise to target pixel measurements
    targetIntenNoise = targetInten + normrnd(0, 0.01, size(targetInten));
    
    %     figure();
    %     plot(targetIntenNoise); hold on;
    %     plot(targetInten)
    
    targetSampleInten(sample) = {targetIntenNoise};
end


%% Combine samples into single array for training
% tic()
% targetInten = [];
% targetPolarLS = [];
%
% for sample = 1:noSamples
%     targetInten = [targetInten, targetSampleInten{sample}];
%     targetPolarLS = [ targetPolarLS, targetSamplePolarLS{sample}];
% end
% toc();

% targetInten = zeros(1,noSamples*yPix);
% targetPolarLS = zeros(2,noSamples*yPix);
% targetPntNormals = zeros(3,noSamples*yPix);
% targetPntFrame = zeros(3, noSamples*yPix);
% noPts = 0;

if useCMAES
    cmaesOptions.LBounds = [0, 0, 0]';
    cmaesOptions.UBounds = [100,100,100]';
    
    cmaesOptions.SaveVariables = 'off';
    cmaesOptions.DispFinal = 'off';
    cmaesOptions.LogModulo = 0;
    cmaesOptions.DispModulo = 0;
    
else
    
    optOptions = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'SpecifyObjectiveGradient',true, 'CheckGradients', false, ...
        'MaxIterations', 1000000000, 'FunctionTolerance',1e-6, 'MaxFunctionEvaluations',1000000000, 'StepTolerance',1e-6);
    %     ,'ScaleProblem', ...
    %         'jacobian', 'InitDamping', 0.01, 'FiniteDifferenceType', 'central');
    
    optOptions.Display = 'none';
    
end


% optOptions = optimoptions('lsqnonlin', 'Algorithm', 'trust-region-reflective', 'SpecifyObjectiveGradient',true, 'CheckGradients', false, ...
%     'MaxIterations', 1000000000, 'FunctionTolerance',1e-6, 'MaxFunctionEvaluations',1000000000, 'StepTolerance',1e-15, 'OptimalityTolerance', 1e-6);
% optOptions.Display = 'none';

optPhiTrials = zeros(noSamples, 3);
resNormTrials = zeros(noSamples, 1);

Phi0 = [maxRadiInt, mu, rAtt];


for trials = 1:noSamples
    
    targetL = zeros(1,trials*yPix);
    targetPolarLS = zeros(2,trials*yPix);
    targetPntNormals = zeros(3,trials*yPix);
    targetPntFrame = zeros(3, trials*yPix);
    noPts = 0;
    
    %combines all samples into a single array for the current trials
    for samples = 1:trials
        noCurPts = length(targetSampleInten{sample});
        
        targetL(noPts+1:noPts+noCurPts) = targetSampleInten{sample};
        targetPolarLS(:, noPts+1:noPts+noCurPts) = targetSamplePolarLS{sample};
        targetPntNormals(:, noPts+1:noPts+noCurPts) = repmat(targetSampleNormal(:, sample), [1, noCurPts]);
        targetPntFrame(:, noPts+1:noPts+noCurPts) = targetSamplePosFrame{sample};
        
        noPts = noPts + noCurPts;
    end
    
    targetL = targetL(1:noPts);
    
    %     targetE = targetE + normrnd(0, 0.1, size(targetE));
    
    targetPolarLS = targetPolarLS(:, 1:noPts);
    targetPntNormals = targetPntNormals(:, 1:noPts);
    targetPntFrame = targetPntFrame(:, 1:noPts);
    
    %The optimisation function
    f = @(H)LightSourceLeastSqrObj(H, targetL', targetPntFrame, ...
        targetPntNormals, ligSourLoc, lightSrc.get_SourDirVec(), target.Reflectance);
    
    
    [optPhi,resNorm] = lsqnonlin(f,Phi0, [], [], optOptions);
    
    % [optPhi,resNorm] = cmaes('LightSourceResidual', Phi0, [], cmaesOptions, targetL', targetPntFrame, ...
    %     targetPntNormals, ligSourLoc, lightSrc.get_SourDirVec(), target.Reflectance);
    
    optPhi = optPhi';
    
    optPhiTrials(trials, :) = optPhi;
    resNormTrials(trials) = resNorm;
end

%%
abs_phi0 = abs(optPhiTrials(:,1) - maxRadiInt);
abs_mu = abs(optPhiTrials(:,2) - mu);
abs_r_d = abs(optPhiTrials(:,3) - rAtt);

figure('Name', 'Absolute error in Phi_0');
plot(1:noSamples, abs_phi0);
xlabel('samples');
ylabel('absolute error in \Phi_0');
ylim([0, max(abs_phi0)])
grid on;

figure('Name', 'Absolute error in mu');
plot(1:noSamples, abs_mu);
xlabel('samples');
ylabel('absolute error in \mu');
ylim([0, max(abs_mu)])
grid on;

figure('Name', 'Optimisation Residuals');
plot(1:noSamples, resNormTrials);
xlabel('samples');
ylabel('Residual after optimisation');
ylim([0, max(resNormTrials)])
grid on;




% %% Non-linear optimisation using parameter light source model
%
% %Calculate incoming radiance recieved through the knowledge of the material
% targetL = (pi/target.Reflectance) .* targetInten;
%
% %The optimisation function
% f = @(H)LightSourceLeastSqrObj(H, targetL', targetPntFrame, ...
%     targetPntNormals, ligSourLoc, lightSrc.get_SourDirVec());
%
% Phi0 = [10, 1, 1];
%
% optOptions = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'SpecifyObjectiveGradient',true, 'CheckGradients', true, ...
%     'MaxIterations', 1000000000, 'FunctionTolerance',1e-6, 'MaxFunctionEvaluations',1000000000, 'StepTolerance',1e-6);%,'ScaleProblem', ...
% %'jacobian', 'InitDamping', 0.01, 'FiniteDifferenceType', 'central');
% optOptions.Display = 'none';
%
% % optOptions = optimoptions('lsqnonlin', 'Algorithm', 'trust-region-reflective', 'SpecifyObjectiveGradient',false, 'CheckGradients', false, ...
% %     'MaxIterations', 1000000000, 'FunctionTolerance',1e-6, 'MaxFunctionEvaluations',1000000000, 'StepTolerance',1e-15, 'OptimalityTolerance', 1e-6);
% % optOptions.Display = 'none';
%
% %Perform optimisation
% [optPhi,resNorm] = lsqnonlin(f,Phi0, [] , [], optOptions);
%
% %% Set up for training
% xAngle = pi;
% rotTarget = eul2rotm([0, 0, xAngle], 'ZYX');
%
% predictorNames = {'theta', 'rho'};
%
% % Convert input to table
% inputTable = array2table(targetPolarLS', 'VariableNames', predictorNames);
%
% predictors = inputTable(:, predictorNames);
% response = targetInten(:);
% isCategoricalPredictor = [false, false];
%
% %% Train a linear regression model
% % This code specifies all the model options and trains the model.
% concatenatedPredictorsAndResponse = predictors;
% concatenatedPredictorsAndResponse.targetInten = response;
% linearModel = fitlm(concatenatedPredictorsAndResponse, 'linear', 'RobustOpts', 'off');
% disp(linearModel);
% % coeffLine = linearModel.Coefficients.Estimate;
%
% %% Train a quadratic regression model
% % This code specifies all the model options and trains the model.
% % concatenatedPredictorsAndResponse = predictors;
% % concatenatedPredictorsAndResponse.targetInten = response;
% QuadModel = fitlm(concatenatedPredictorsAndResponse, 'quadratic', 'RobustOpts', 'off');
% disp(QuadModel);
%
% % coeffQuad = QuadModel.Coefficients.Estimate;
%
% %% Train a GP using maternal kernal
%
% % Train a regression model
% % This code specifies all the model options and trains the model.
% GPmatModel = fitrgp(predictors, ...
%     response, ...
%     'BasisFunction', 'constant', ...
%     'KernelFunction', 'matern52', ...
%     'Standardize', true);
% disp(GPmatModel);
%
% % GPmatModel = fitlm(concatenatedPredictorsAndResponse, 'g
%
% %% 2D plot of radiant intensity on view-plane
%
% %plot it as a rectangle starting at where the line-scan is positioned and
% %only plot the YZ-plane
% %size of plane in each dimension
% % yDist = 0.5;
% % zDist = 1;
%
% %Create mesh
% y = linspace(-yDist, yDist, 100);
% z = linspace(0, zDist, 100);
% x = 0;
%
% [X,Y,Z] = meshgrid(x,y,z);
%
% %remove extra necessary singular dimension
% X = squeeze(X);
% Y = squeeze(Y);
% Z = squeeze(Z);
%
% Yarr = Y(:);
% Zarr = Z(:);
%
% [theta, rho] = cart2pol(Yarr, Zarr);
%
% thetaRhoMeshPred = [theta, rho];
%
% rows = size(X);
%
% radIntMagLine = predict(linearModel, thetaRhoMeshPred);
%
% % a = coeffLine(1);
% % b = coeffLine(2);
% % c = coeffLine(3);
% % radIntMagLine = a + b.*theta + c.*rho;
% radIntMagLine = reshape(radIntMagLine,rows);
%
%
%
% radIntMagQuad = predict(QuadModel, thetaRhoMeshPred);
%
% % a = coeffQuad(1);
% % b = coeffQuad(2);
% % c = coeffQuad(3);
% % d = coeffQuad(4);
% % e = coeffQuad(5);
% % f = coeffQuad(6);
% % radIntMagQuad = a + b.*theta + c.*rho + d.*theta.*rho + e.*theta.^2 + f.*rho.^2;
% radIntMagQuad = reshape(radIntMagQuad,rows);
%
% radIntMagGP = predict(GPmatModel, thetaRhoMeshPred);
% radIntMagGP = reshape(radIntMagGP,rows);
%
%
% %plot 2D view-plane surface with line-scan camera as origin. Only plotting
% %YZ plane.
% figLineReg = figure('Name', 'Radiant intensity on line-scan view-plane using LINEAR regression');
% surf(Y, Z, radIntMagLine, 'EdgeColor', 'none'); hold on;
% xlabel('y');
% ylabel('z');
% title('Linear regression model')
% view(2);
% scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled'); %line-scan origin
% colormap(jet(100));
% colorbar; caxis([0, maxRadiInt]);
% axis equal;
% %plot green light source
% scatter3(lightSrcLS(2), lightSrcLS(3), (maxRadiInt+1), 200, [0,1,0], 'filled');
% %plot frame camera
% scatter3(frameCameraLS(2,4), frameCameraLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled');
% %plot fov of line-scan on view-plane
% plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
%
%
%
%
% %plot 2D view-plane surface with line-scan camera as origin. Only plotting
% %YZ plane.
% figQuadReg = figure('Name', 'Radiant intensity on line-scan view-plane using QUADATRIC regression');
% surf(Y, Z, radIntMagQuad, 'EdgeColor', 'none'); hold on;
% xlabel('y');
% ylabel('z');
% title('Quadratic regression model')
% view(2);
% scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled'); %line-scan origin
% colormap(jet(100));
% colorbar; caxis([0, maxRadiInt]);
% axis equal;
% %plot green light source
% scatter3(lightSrcLS(2), lightSrcLS(3), (maxRadiInt+1), 200, [0,1,0], 'filled');
% %plot frame camera
% scatter3(frameCameraLS(2,4), frameCameraLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled');
% %plot fov of line-scan on view-plane
% plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
%
%
% %plot 2D view-plane surface with line-scan camera as origin. Only plotting
% %YZ plane.
% figGPReg = figure('Name', 'Radiant intensity on line-scan view-plane using GP matern 5/2 regression');
% surf(Y, Z, radIntMagGP, 'EdgeColor', 'none'); hold on;
% xlabel('y');
% ylabel('z');
% title('GP matern 5/2 regression model')
% view(2);
% scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled'); %line-scan origin
% colormap(jet(100));
% colorbar; caxis([0, maxRadiInt]);
% axis equal;
% %plot green light source
% scatter3(lightSrcLS(2), lightSrcLS(3), (maxRadiInt+1), 200, [0,1,0], 'filled');
% %plot frame camera
% scatter3(frameCameraLS(2,4), frameCameraLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled');
% %plot fov of line-scan on view-plane
% plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
%
% %% Differences between groundtruth
%
% diffRadIntLine = abs(radIntMagLine - radIntMagDotTarget);
% diffRadIntQuad = abs(radIntMagQuad - radIntMagDotTarget);
% diffRadIntGP = abs(radIntMagGP - radIntMagDotTarget);
%
% figure('Name', 'abosolute difference of view-plane linear');
% surf(Y, Z, diffRadIntLine, 'EdgeColor', 'none'); hold on;
% xlabel('y');
% ylabel('z');
% title('absolute difference of view-plane linear')
% view(2);
% scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled'); %line-scan origin
% colormap(flipud(hot));
% colorbar; caxis([0, maxRadiInt]);
% axis equal;
% %plot green light source
% scatter3(lightSrcLS(2), lightSrcLS(3), (maxRadiInt+1), 200, [0,1,0], 'filled');
% %plot frame camera
% scatter3(frameCameraLS(2,4), frameCameraLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled');
% %plot fov of line-scan on view-plane
% plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
%
%
% figure('Name', 'abosolute difference of view-plane quadratic');
% surf(Y, Z, diffRadIntQuad, 'EdgeColor', 'none'); hold on;
% xlabel('y');
% ylabel('z');
% title('absolute difference of view-plane quadratic')
% view(2);
% scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled'); %line-scan origin
% colormap(flipud(hot));
% colorbar; caxis([0, maxRadiInt]);
% axis equal;
% %plot green light source
% scatter3(lightSrcLS(2), lightSrcLS(3), (maxRadiInt+1), 200, [0,1,0], 'filled');
% %plot frame camera
% scatter3(frameCameraLS(2,4), frameCameraLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled');
% %plot fov of line-scan on view-plane
% plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
%
%
% figure('Name', 'abosolute difference of view-plane GP');
% surf(Y, Z, diffRadIntGP, 'EdgeColor', 'none'); hold on;
% xlabel('y');
% ylabel('z');
% title('absolute difference of view-plane GP')
% view(2);
% scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled'); %line-scan origin
% colormap(flipud(hot));
% colorbar; caxis([0, maxRadiInt]);
% axis equal;
% %plot green light source
% scatter3(lightSrcLS(2), lightSrcLS(3), (maxRadiInt+1), 200, [0,1,0], 'filled');
% %plot frame camera
% scatter3(frameCameraLS(2,4), frameCameraLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled');
% %plot fov of line-scan on view-plane
% plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);

%% Call_back function

function  distFromSrc_callback(src, ~, lightSrc, fig, S)
lightSrc.PlotRadiantSlice(fig, S, src.Value);
end