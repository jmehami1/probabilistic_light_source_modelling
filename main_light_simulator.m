%Initial script for simulating a non-isotropic point light source.

close all;
clear;

%robotics toolbox
run(['rvctools' filesep 'startup_rvc.m'])

%code for 3D line plane intersection
addpath(genpath('line_plane_intersection'));

%code for this project
addpath('code');

%GPML toolbox
run(['gpml-matlab-master', filesep, 'startup.m']);


%% Parameters to change
getSymmetricPts = true;

%variance of pixel intensity noise in line-scan camera measurements
intNoiseSigma = 0.01;

%% Setting up non-isotropic light source

%Pose of source in world coordinates (frame camera)
locLightSrc = [0.1;0.1;-0.1];
rotLightSrc = eul2rotm(deg2rad([0,-10, 0]), 'ZYX');
poseLightSrc = [rotLightSrc, locLightSrc; 0, 0, 0, 1];

%radiant intensity distribution parameters
maxRadiInt = 10;
mu = 1.5;

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
lightSrc = LightSimulator(locLightSrc, rotLightSrc, maxRadiInt, mu, r_d, distFromSrc, figSim, S);
trplot(poseLightSrc, 'rviz', 'length', 0.1);

%distFromSrc slider
b = uicontrol('Parent',figSim,'Style','slider','Position',[81,120,300,23],...
    'value',distFromSrc, 'min', 0, 'max',1, 'SliderStep', [0.1/1, 0.1/1]);
bgcolor = figSim.Color;
uicontrol('Parent',figSim,'Style','text','Position',[50,90,23,23],...
    'String','0','BackgroundColor',bgcolor);
uicontrol('Parent',figSim,'Style','text','Position',[350,90,23,23],...
    'String','1','BackgroundColor',bgcolor);
uicontrol('Parent',figSim,'Style','text','Position',[150,90,200,23],...
    'String','Distance from the source','BackgroundColor',bgcolor);
disp_dist = uicontrol('Parent',figSim,'Style','text','Position',[200,145,50,20],...
    'String', num2str(distFromSrc),'BackgroundColor', [1,1,1]);
%callback function at the end of the script
b.Callback = @(src, eventData) lightSrc.distFromSrc_callback(src, eventData, disp_dist);


%mu slider
slider_mu = uicontrol('Parent',figSim,'Style','slider','Position',[81,200,300,23],...
    'value', mu, 'min',0, 'max',10, 'SliderStep', [0.5/10, 0.5/10]);
bgcolor = figSim.Color;
uicontrol('Parent',figSim,'Style','text','Position',[50,170,23,23],...
    'String','0','BackgroundColor',bgcolor);
uicontrol('Parent',figSim,'Style','text','Position',[350,170,23,23],...
    'String','10','BackgroundColor',bgcolor);
uicontrol('Parent',figSim,'Style','text','Position',[150,170,200,23],...
    'String','mu','BackgroundColor',bgcolor);
disp_mu = uicontrol('Parent',figSim,'Style','text','Position',[200,225,50,20],...
    'String', num2str(mu),'BackgroundColor', [1,1,1]);
%callback function at the end of the script
slider_mu.Callback = @(src, eventData) lightSrc.mu_callback(src, eventData, disp_mu);

%checkerbox to turn surface plane visibility on/off
surfDirDispCheck = uicontrol('Parent',figSim,'Style','checkbox', 'String', 'Display Direction Plane', 'Position', [20,45,200,20] );
surfDirDispCheck.Callback = @(src, eventData) lightSrc.surfDir_visible_callback(src, eventData);

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
tLS = [-0.1; 0.1; -0.1];
poseLS = [rotLS, tLS; 0, 0, 0, 1];
extLS = poseLS \ eye(4);

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
colormap(hot);
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

%for storing symmetric points (if required)
if getSymmetricPts
    targetSamplePosSymmFrame = cell(1,noSamples);
    targetSampleSymmNormal = zeros(3, noSamples);
end

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
    if getSymmetricPts
        targetSymmPtsLS = zeros(3, length(vtargetPix));
        targetCentSymmPos = SymmetricPoint_LightSrc_Linescan(targetCentPos, locLightSrcLS, rotLightSrcDirLS, vpNormal);
    end
    
    for pixel = 1:length(vtargetPix)
        pnt = normTargetPixHom(:,pixel);
        
        %calculate 3D point of pixel on target relative to line-scan camera
        %using line-plane intersection
        [pntOnTarget, rc] = line_plane_intersection(pnt, pnt, normTargetLS, targetCentPos);
        
        %the point should be on the plane
        if rc ~= 1
            error('line-plane intersection was not on plane');
        end
        
        targetPtsLS(:, pixel) = pntOnTarget;
        
        %find symmetric point about light source direction on view-plane
        if getSymmetricPts
            pntOnTargetSymm = SymmetricPoint_LightSrc_Linescan(pntOnTarget, locLightSrcLS, rotLightSrcDirLS, vpNormal);
            targetSymmPtsLS(:, pixel) = pntOnTargetSymm;
        end
    end
    
    %calculate the symmetric normal, which is actually finding the rotation
    %of the symmetrical target and taking its z-directional axis
    if getSymmetricPts
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
    end
    
    %transform points from line-scan to frame camera
    targetPtsFrame = poseLS*[targetPtsLS; ones(1,length(targetPtsLS(1,:)))];
    targetPtsFrame = targetPtsFrame(1:3,:);
    
    
    %plot line and axes of target in the simulator figure
    figure(figSim);
    line(targetPtsFrame(1,:), targetPtsFrame(2,:), targetPtsFrame(3,:), 'Color', [0,0,0], 'LineWidth', 2); hold on;
    trplot(poseTargetFrame, 'rviz', 'length', (maxRange - minRange)*0.1);
    
    %plot line and axes of symmetric target in the simulator figure
    if getSymmetricPts
        line(targetSymmPtsFrame(1,:), targetSymmPtsFrame(2,:), targetSymmPtsFrame(3,:), 'Color', [0,0,0], 'LineWidth', 2);
        trplot(targetSymmPoseFrame, 'rviz', 'length', (maxRange - minRange)*0.1);
    end
    
    drawnow();
    
    %plot line in the view-plane intensity plot
    figure(figViewPlane);
    line(targetPtsLS(2,:), targetPtsLS(3,:), maxRadiInt.*ones(size(targetPtsLS(2,:))), 'Color', [0,0,0], 'LineWidth', 2);
    
    %plot symmetric line in the view-plane intensity plot
    if getSymmetricPts
        line(targetSymmPtsLS(2,:), targetSymmPtsLS(3,:), maxRadiInt.*ones(size(targetSymmPtsLS(2,:))), 'Color', [0,0,0], 'LineWidth', 2);
    end
    
    drawnow();
    
    %save the target pose in the c.f of the frame camera
    targetSamplePosFrame(sample) = {targetPtsFrame};
    
    %save symmetric target pose in the c.f of the frame camera
    if getSymmetricPts
        targetSamplePosSymmFrame(sample) = {targetSymmPtsFrame};
    end
end

% %% 2D plot of radiant intensity on view-plane dot product with normal
%
% % plot it as a rectangle starting at where the line-scan is positioned and
% % only plot the YZ-plane
% % size of plane in each dimension
% yDist = 0.5;
% zDist = 1;
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
% %These points are in the line-scan coordinate frame, transform them into
% %the frame camera coordinate frame (world coordinates
% %rotation
% XYZrot = [X(:),Y(:),Z(:)]*rotLS;
%
% rows = size(X);
% %reshape to mesh and translate
% XFrame = reshape(XYZrot(:,1),rows) + tLS(1);
% YFrame = reshape(XYZrot(:,2),rows) + tLS(2);
% ZFrame = reshape(XYZrot(:,3),rows) + tLS(3);
%
% % [~, radIntMagVec] = lightSrc.RadiantIntensityMesh(XFrame, YFrame, ZFrame);
%
% radianceOutTarget = zeros(size(XFrame));
%
% targetNormal = targetSampleNormal(:, 1);
%
% for i = 1:size(XFrame,1)
%     for j = 1:size(XFrame, 2)
%         pnt = [XFrame(i,j), YFrame(i,j), ZFrame(i,j)];
%
%         radianceOutTarget(i, j) = lightSrc.RadianceOutMaterialPoint(pnt, targetNormal, target.Reflectance);
%
%     end
% end
%
% radianceOutTargetNoise = radianceOutTarget + normrnd(0, intNoiseSigma, size(radianceOutTarget));
%
% diffRadiance = abs(radianceOutTargetNoise - radianceOutTarget);
%
% maxRadianceOut = maxRadiInt/pi;
%
% %plot 2D view-plane surface with line-scan camera as origin. Only plotting
% %YZ plane.
% figRadianceOut = figure('Name', 'Radiance emitted out of surface with constant surface normal');
% surf(Y, Z, radianceOutTarget, 'EdgeColor', 'none'); hold on;
% xlabel('y');
% ylabel('z');
% title('View-plane RADIANCE emitted from target');
% view(2);
% scatter3(0, 0, (maxRadianceOut+1), 200, [0,0,1], 'filled'); %line-scan origin
% colormap(cool);
% colorbar; caxis([0, maxRadianceOut]);
% axis equal;
%
% %plot green light source
% poseLightSrcLS = T_LS2Frame \ [locLightSrc; 1]; %transform to line-scan coordinate system
% scatter3(poseLightSrcLS(2), poseLightSrcLS(3), (maxRadianceOut+1), 200, [0,1,0], 'filled');
%
% %plot frame camera
% extLightSrc = inv(T_LS2Frame);
% scatter3(extLightSrc(2,4), extLightSrc(3,4), (maxRadianceOut+1), 200, [1,0,0], 'filled');
%
% %plot fov of line-scan on view-plane
% plot3(yPoly, zPoly, (maxRadianceOut+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
%
%
%
% %plot 2D view-plane surface with line-scan camera as origin. Only plotting
% %YZ plane.
% figRadianceOutNoise = figure('Name', 'Radiance emitted out of surface with constant surface normal and added noise');
% surf(Y, Z, radianceOutTargetNoise, 'EdgeColor', 'none'); hold on;rows = size(X);
% %reshape to mesh and translate
% XFrame = reshape(XYZrot(:,1),rows) + tLS(1);
% xlabel('y');
% ylabel('z');
% title(['View-plane RADIANCE emitted from target with added noise of \sigma ', num2str(intNoiseSigma)]);
% view(2);
% scatter3(0, 0, (maxRadianceOut+1), 200, [0,0,1], 'filled'); %line-scan origin
% colormap(cool);
% colorbar; caxis([0, maxRadianceOut]);
% axis equal;
%
% %plot green light source
% scatter3(poseLightSrcLS(2), poseLightSrcLS(3), (maxRadianceOut+1), 200, [0,1,0], 'filled');
%
% %plot frame camera
% scatter3(extLightSrc(2,4), extLightSrc(3,4), (maxRadianceOut+1), 200, [1,0,0], 'filled');
%
% %plot fov of line-scan on view-plane
% plot3(yPoly, zPoly, (maxRadianceOut+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
%
% maxRadianceOut = max(diffRadiance, [], 'all')*1.1;
%
% %plot 2D view-plane surface with line-scan camera as origin. Only plotting
% %YZ plane.
% figRadianceDiff = figure('Name', 'Difference in Radiance from surface with/without noise');
% surf(Y, Z, diffRadiance, 'EdgeColor', 'none'); hold on;
% xlabel('y');
% ylabel('z');
% title('Difference in view-plane RADIANCE emitted from targetwith/without noise ');
% view(2);
% scatter3(0, 0, (maxRadianceOut+1), 200, [0,0,1], 'filled'); %line-scan origin
% colormap(hot);
% colorbar; caxis([0, maxRadianceOut]);
% axis equal;
%
% %plot green light source
% scatter3(poseLightSrcLS(2), poseLightSrcLS(3), (maxRadianceOut+1), 200, [0,1,0], 'filled');
%
% %plot frame camera
% scatter3(extLightSrc(2,4), extLightSrc(3,4), (maxRadianceOut+1), 200, [1,0,0], 'filled');
%
% %plot fov of line-scan on view-plane
% plot3(yPoly, zPoly, (maxRadianceOut+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);

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
        %         angleLightRayLS = targetPolarLS(1,pt);
        %         vignet = betaLS*cos(angleLightRayLS)^4;
        %         targetInten(pt) = lightSrc.RadianceOutMaterialPoint(targetPtsFrame(:,pt), targetNormal, target.Reflectance)
        
        targetInten(pt) = lightSrc.RadianceOutMaterialPoint(targetPtsFrame(:,pt), targetNormal, target.Reflectance);
    end
    
    %add gaussian white noise to target pixel measurements
    targetIntenNoise = targetInten + normrnd(0, intNoiseSigma, size(targetInten));
    targetSampleInten(sample) = {targetIntenNoise};
end

%% Combine the samples such that they can be analysed by regressional approaches

% {noSamples} x {targetL, targetPntNormals, targetPntFrame,, targetPntNormalsSymm, targetPntFrameSymm}
targetTrials = cell(noSamples, 5);

%number of trials is equal to the number of samples. Each consective trial
%will add an extra sample
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
    
    targetTrials(trials,:) = {targetL, targetPntNormals, targetPntFrame, targetPntNormalsSymm, targetPntFrameSymm};
    
%     targetTrials(trials,1) = {targetL};
%     targetTrials(trials,2) = {targetPntNormals};
%     targetTrials(trials,3) = {targetPntFrame};
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

% Phi0 = [maxRadiInt, mu, rAtt];
Phi0 = [5, 1, 1];

for trials = 1:noSamples
    
%     if getSymmetricPts
%         targetL = zeros(1,2*trials*yPix);
%         targetPntNormals = zeros(3,2*trials*yPix);
%         targetPntFrame = zeros(3, 2*trials*yPix);
%     else
%         targetL = zeros(1, trials*yPix); %#ok<*UNRCH>
%         targetPntNormals = zeros(3, trials*yPix);
%         targetPntFrame = zeros(3, trials*yPix);
%     end
%     
%     noPts = 0;
%     
%     %combines all samples into a single array for the current trials
%     for samples = 1:trials
%         noCurPts = length(targetSampleInten{sample});
%         
%         targetL(noPts+1:noPts+noCurPts) = targetSampleInten{sample};
%         targetPntNormals(:, noPts+1:noPts+noCurPts) = repmat(targetSampleNormal(:, sample), [1, noCurPts]);
%         targetPntFrame(:, noPts+1:noPts+noCurPts) = targetSamplePosFrame{sample};
%         
%         noPts = noPts + noCurPts;
%         
%         if getSymmetricPts
%             targetL(noPts+1:noPts+noCurPts) = targetSampleInten{sample};
%             targetPntNormals(:, noPts+1:noPts+noCurPts) = repmat(targetSampleSymmNormal(:, sample), [1, noCurPts]);
%             targetPntFrame(:, noPts+1:noPts+noCurPts) = targetSamplePosSymmFrame{sample};
%             
%             noPts = noPts + noCurPts;
%         end
%     end
    
%     %clip the arrays to correct size
%     targetL = targetL(1:noPts);
%     targetPntNormals = targetPntNormals(:, 1:noPts);
%     targetPntFrame = targetPntFrame(:, 1:noPts);

%     targetTrials(trials,:) = {targetL, targetPntNormals, targetPntFrame, targetPntNormalsSymm, targetPntFrameSymm};

    targetL = targetTrials{trials,1};
    targetL = [targetL, targetL]; %#ok<AGROW>
    
    targetPntFrame = [targetTrials{trials,3}, targetTrials{trials,5}];
    targetPntNormals = [targetTrials{trials,2}, targetTrials{trials,4}];

    
    %The optimisation function
    f = @(H)LightSourceLeastSqrObj(H, targetL', targetPntFrame, ...
        targetPntNormals, locLightSrc, lightSrc.get_SourDirVec(), target.Reflectance);
    
    %Perform optimisation
    [optPhi,resNorm] = lsqnonlin(f,Phi0, [], [], optOptions);
    optPhi = optPhi';
    
    %store optimised parameters and residual
    optPhiTrials(trials, :) = optPhi;
    resNormTrials(trials) = resNorm;
    
    %determine observability of NLS by calculating the Fisher Information
    %Matrix
    %Jacobian of radiance objective function w.r.t to unknown parameters
    jac = JacobianLightRadiance(optPhi, targetPntFrame, ...
        targetPntNormals, locLightSrc, lightSrc.get_SourDirVec(), target.Reflectance);
    
    %create covariance matrix
    cov = diag((1/(intNoiseSigma^2)).*ones(1,size(jac,1)));
    
    FIM = (jac'/cov)*jac;
    
    rankFIM(trials) = rank(FIM);
    eigFIM(trials,:) = log(eig(FIM))';
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
s1 = subplot(1,3,1);
surf(Y, Z, radIntMagGroundtruth, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Ground-truth')
view(2);
scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s1, hot);
caxis([0, maxRadiInt]);
% colorbar; 
% axis equal;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);



s2 = subplot(1,3,2);
surf(Y, Z, radIntMagLeastSqr, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Least Squares')
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

s3 = subplot(1,3,3);
surf(Y, Z, diffRad, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Absolute difference in Least Squares')
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

%% Building model with GP non-symmetric

cab(figSim, figViewPlane, figAbsLS, figLeastSqrViewPlane);


meanfunc = [];                    % empty: don't use a mean function
% covfunc = {@covSE, 'eye', []};              % Squared Exponental covariance function
% covfunc = {@covMaterniso, 10};              % Squared Exponental covariance function
covfunc = @covSEiso;

likfunc = @likGauss;              % Gaussian likelihood
hyp = struct('mean', [], 'cov', [0,0], 'lik', -1);


ptsHomLightSrc = poseLightSrc\ptsHomFrame;
ptsLightSrc = ptsHomLightSrc(1:3, :);

%radius of point from source
% ptsRadius = sqrt(sum(ptsLightSrc.^2, 1));

[ptsTheta, ptsRadius] = cart2pol(ptsLightSrc(1,:), ptsLightSrc(3,:));
ptsTheta = pi/2 - ptsTheta;

%angle from the light source direction vector
% srcDir = [0;0;1];
% ptsTheta = acos(sum(srcDir.*ptsLightSrc,1)./ptsRadius);

testingX = [ptsRadius', ptsTheta'];

figGP_NonSymmetric = figure('Name','GP training Non-Symmetric');
s1 = subplot(1,2,1);
surfGP = surf(Y, Z, zeros(size(Y)), 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Zero-mean NOT-symmetric points GP')
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

s2 = subplot(1,2,2);
surfGP_var = surf(Y, Z, zeros(size(Y)), 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Variance')
view(2);
scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s2, turbo);
caxis([0, 3]);
colorbar; 

%plot light source in the line-scan camera coordinate frame (green)
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);


for trials = noSamples
    
    targetL = targetTrials{trials,1};
    
    targetPntFrame = targetTrials{trials,3};
    targetPntNormals = targetTrials{trials,2};

    %calculate the direction light vector (point to light source)
    pntLightVec = locLightSrc - targetPntFrame;
    dirPntLightVec = pntLightVec./vecnorm(pntLightVec);
            
    %Calculate radiant intensity magnitude used for building model
    radIntMagPnt = (targetL.*pi)./(target.Reflectance.*dot(targetPntNormals, dirPntLightVec,1));
    
    %find the point w.r.t to light source c.f
    targetPntFrameHom = [targetPntFrame; ones(1, length(targetPntFrame(1,:)))];
    targetPntLightSrcHom = poseLightSrc\targetPntFrameHom;
    targetPntLightSrc = targetPntLightSrcHom(1:3, :);
    
    [targetPntTheta, targetPntRadius] = cart2pol(targetPntLightSrc(1,:), targetPntLightSrc(3,:));
    targetPntTheta = pi/2 - targetPntTheta;
    
    GPdata = [targetPntRadius', targetPntTheta', radIntMagPnt'];
    
    %downsample data
    GPDataDown = downsample(GPdata, 50);
    xGP = GPDataDown(:,1:2);
    yGP = GPDataDown(:,3);
    
    %Plotting training points on figure
    [xSrc, zSrc] = pol2cart(pi/2 - xGP(:,2), xGP(:,1));
    ptSrcHom = [xSrc'; zeros(1,length(zSrc)); zSrc'; ones(1,length(zSrc))];   
    ptSrcLSHom = poseLightSrcLS*ptSrcHom;
    
    figure(figGP_NonSymmetric);
    
    %plot light source in the line-scan camera coordinate frame (green)
    scatter3(ptSrcLSHom(2,:), ptSrcLSHom(3,:), (maxRadiInt+1)*ones(size(ptSrcLSHom(3,:))), 20, [0,1,1], 'filled','MarkerEdgeColor', [0,0,0]);
    
    %perform training
    hyp2 = minimize(hyp, @gp, -1000, @infGaussLik, meanfunc, covfunc, likfunc, xGP, yGP);

    [mu, varMu] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, xGP, yGP, testingX);

    radIntGP = reshape(mu,rows);
    radVar = reshape(varMu,rows);
    
    %update figure with GP model
    surfGP.ZData = radIntGP;
    surfGP_var.ZData = radVar;
    drawnow();
end


%plot 2D view-plane surface with line-scan camera as origin. Only plotting
%YZ plane.
figLeastGPViewPlaneNonSymm = figure('Name', 'Radiant intensity view-plane from GP');
s1 = subplot(1,3,1);
surf(Y, Z, radIntMagGroundtruth, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Ground-truth')
view(2);
scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s1, hot);
caxis([0, maxRadiInt]);
% colorbar; 
% axis equal;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);



s2 = subplot(1,3,2);
surf(Y, Z, radIntGP, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Zero-mean GP NOT Symmetric')
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

s3 = subplot(1,3,3);
surf(Y, Z, diffRad, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Absolute difference in GP')
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


%% Building model with GP Symmetric

cab(figSim, figViewPlane, figAbsLS, figLeastSqrViewPlane, figGP_NonSymmetric, figLeastGPViewPlaneNonSymm);


meanfunc = [];                    % empty: don't use a mean function
% covfunc = {@covSE, 'eye', []};              % Squared Exponental covariance function
% covfunc = {@covMaterniso, 10};              % Squared Exponental covariance function
covfunc = @covSEiso;

likfunc = @likGauss;              % Gaussian likelihood
hyp = struct('mean', [], 'cov', [0,0], 'lik', -1);


ptsHomLightSrc = poseLightSrc\ptsHomFrame;
ptsLightSrc = ptsHomLightSrc(1:3, :);

%radius of point from source
% ptsRadius = sqrt(sum(ptsLightSrc.^2, 1));

[ptsTheta, ptsRadius] = cart2pol(ptsLightSrc(1,:), ptsLightSrc(3,:));
ptsTheta = pi/2 - ptsTheta;

%angle from the light source direction vector
% srcDir = [0;0;1];
% ptsTheta = acos(sum(srcDir.*ptsLightSrc,1)./ptsRadius);

testingX = [ptsRadius', ptsTheta'];

figGPSymm = figure('Name','GP training');
s1 = subplot(1,2,1);
surfGP = surf(Y, Z, zeros(size(Y)), 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Zero-mean symmetric points GP')
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

s2 = subplot(1,2,2);
surfGP_var = surf(Y, Z, zeros(size(Y)), 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Variance')
view(2);
scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s2, turbo);
caxis([0, 3]);
colorbar; 

%plot light source in the line-scan camera coordinate frame (green)
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);


for trials = 1:noSamples
    
    targetL = targetTrials{trials,1};
    targetL = [targetL, targetL]; %#ok<AGROW>
    
    targetPntFrame = [targetTrials{trials,3}, targetTrials{trials,5}];
    targetPntNormals = [targetTrials{trials,2}, targetTrials{trials,4}];

    %calculate the direction light vector (point to light source)
    pntLightVec = locLightSrc - targetPntFrame;
    dirPntLightVec = pntLightVec./vecnorm(pntLightVec);
            
    %Calculate radiant intensity magnitude used for building model
    radIntMagPnt = (targetL.*pi)./(target.Reflectance.*dot(targetPntNormals, dirPntLightVec,1));
    
    
    %find the point w.r.t to light source c.f
    targetPntFrameHom = [targetPntFrame; ones(1, length(targetPntFrame(1,:)))];
    targetPntLightSrcHom = poseLightSrc\targetPntFrameHom;
    targetPntLightSrc = targetPntLightSrcHom(1:3, :);
    
    
    [targetPntTheta, targetPntRadius] = cart2pol(targetPntLightSrc(1,:), targetPntLightSrc(3,:));
    targetPntTheta = pi/2 - targetPntTheta;
    
    GPdata = [targetPntRadius', targetPntTheta', radIntMagPnt'];
    
    %downsample data
    GPDataDown = downsample(GPdata, 50);
    xGP = GPDataDown(:,1:2);
    yGP = GPDataDown(:,3);
    
    %Plotting training points on figure
    [xSrc, zSrc] = pol2cart(pi/2 - xGP(:,2), xGP(:,1));
    ptSrcHom = [xSrc'; zeros(1,length(zSrc)); zSrc'; ones(1,length(zSrc))];   
    ptSrcLSHom = poseLightSrcLS*ptSrcHom;
    
    figure(figGPSymm)
    
    %plot light source in the line-scan camera coordinate frame (green)
    scatter3(ptSrcLSHom(2,:), ptSrcLSHom(3,:), (maxRadiInt+1)*ones(size(ptSrcLSHom(3,:))), 20, [0,1,1], 'filled','MarkerEdgeColor', [0,0,0]);
    
    %perform training
    hyp2 = minimize(hyp, @gp, -1000, @infGaussLik, meanfunc, covfunc, likfunc, xGP, yGP);

    [mu, varMu] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, xGP, yGP, testingX);

    radIntGP = reshape(mu,rows);
    radVar = reshape(varMu,rows);
    
    %update figure with GP model
    surfGP.ZData = radIntGP;
    surfGP_var.ZData = radVar;
    drawnow();
end


%plot 2D view-plane surface with line-scan camera as origin. Only plotting
%YZ plane.
figLeastGPViewPlaneSymm = figure('Name', 'Radiant intensity view-plane from GP');
s1 = subplot(1,3,1);
surf(Y, Z, radIntMagGroundtruth, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Ground-truth')
view(2);
scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s1, hot);
caxis([0, maxRadiInt]);
% colorbar; 
% axis equal;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);



s2 = subplot(1,3,2);
surf(Y, Z, radIntGP, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Const Mean GP Symmetric')
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

s3 = subplot(1,3,3);
surf(Y, Z, diffRad, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Absolute difference in GP')
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

%% Building model with GP Symmetric Constant Mean

cab(figSim, figViewPlane, figAbsLS, figLeastSqrViewPlane, figGP_NonSymmetric, ...
    figLeastGPViewPlaneNonSymm, figGPSymm, figLeastGPViewPlaneSymm);


meanfunc = @meanConst;                    % empty: don't use a mean function
% covfunc = {@covSE, 'eye', []};              % Squared Exponental covariance function
% covfunc = {@covMaterniso, 10};              % Squared Exponental covariance function
covfunc = @covSEiso;

likfunc = @likGauss;              % Gaussian likelihood
hyp = struct('mean', [1], 'cov', [0,0], 'lik', -1);


ptsHomLightSrc = poseLightSrc\ptsHomFrame;
ptsLightSrc = ptsHomLightSrc(1:3, :);

%radius of point from source
% ptsRadius = sqrt(sum(ptsLightSrc.^2, 1));

[ptsTheta, ptsRadius] = cart2pol(ptsLightSrc(1,:), ptsLightSrc(3,:));
ptsTheta = pi/2 - ptsTheta;

%angle from the light source direction vector
% srcDir = [0;0;1];
% ptsTheta = acos(sum(srcDir.*ptsLightSrc,1)./ptsRadius);

testingX = [ptsRadius', ptsTheta'];

figGPSymmConst = figure('Name','GP training');
s1 = subplot(1,2,1);
surfGP = surf(Y, Z, zeros(size(Y)), 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Const mean symmetric points GP')
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

s2 = subplot(1,2,2);
surfGP_var = surf(Y, Z, zeros(size(Y)), 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Variance')
view(2);
scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s2, turbo);
caxis([0, 3]);
colorbar; 

%plot light source in the line-scan camera coordinate frame (green)
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);


for trials = 1:noSamples
    
    targetL = targetTrials{trials,1};
    targetL = [targetL, targetL]; %#ok<AGROW>
    
    targetPntFrame = [targetTrials{trials,3}, targetTrials{trials,5}];
    targetPntNormals = [targetTrials{trials,2}, targetTrials{trials,4}];

    %calculate the direction light vector (point to light source)
    pntLightVec = locLightSrc - targetPntFrame;
    dirPntLightVec = pntLightVec./vecnorm(pntLightVec);
            
    %Calculate radiant intensity magnitude used for building model
    radIntMagPnt = (targetL.*pi)./(target.Reflectance.*dot(targetPntNormals, dirPntLightVec,1));
    
    
    %find the point w.r.t to light source c.f
    targetPntFrameHom = [targetPntFrame; ones(1, length(targetPntFrame(1,:)))];
    targetPntLightSrcHom = poseLightSrc\targetPntFrameHom;
    targetPntLightSrc = targetPntLightSrcHom(1:3, :);
    
    
    [targetPntTheta, targetPntRadius] = cart2pol(targetPntLightSrc(1,:), targetPntLightSrc(3,:));
    targetPntTheta = pi/2 - targetPntTheta;
    
    GPdata = [targetPntRadius', targetPntTheta', radIntMagPnt'];
    
    %downsample data
    GPDataDown = downsample(GPdata, 50);
    xGP = GPDataDown(:,1:2);
    yGP = GPDataDown(:,3);
    
    %Plotting training points on figure
    [xSrc, zSrc] = pol2cart(pi/2 - xGP(:,2), xGP(:,1));
    ptSrcHom = [xSrc'; zeros(1,length(zSrc)); zSrc'; ones(1,length(zSrc))];   
    ptSrcLSHom = poseLightSrcLS*ptSrcHom;
    
    figure(figGPSymmConst)
    
    %plot light source in the line-scan camera coordinate frame (green)
    scatter3(ptSrcLSHom(2,:), ptSrcLSHom(3,:), (maxRadiInt+1)*ones(size(ptSrcLSHom(3,:))), 20, [0,1,1], 'filled','MarkerEdgeColor', [0,0,0]);
    
    %perform training
    hyp2 = minimize(hyp, @gp, -1000, @infGaussLik, meanfunc, covfunc, likfunc, xGP, yGP);

    [mu, varMu] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, xGP, yGP, testingX);

    radIntGP = reshape(mu,rows);
    radVar = reshape(varMu,rows);
    
    %update figure with GP model
    surfGP.ZData = radIntGP;
    surfGP_var.ZData = radVar;
    drawnow();
end


%plot 2D view-plane surface with line-scan camera as origin. Only plotting
%YZ plane.
figLeastGPViewPlaneSymmConst = figure('Name', 'Radiant intensity view-plane from GP');
s1 = subplot(1,3,1);
surf(Y, Z, radIntMagGroundtruth, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Ground-truth')
view(2);
scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin
colormap(s1, hot);
caxis([0, maxRadiInt]);
% colorbar; 
% axis equal;

%plot light source in the line-scan camera coordinate frame (green)
scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), (maxRadiInt+1), 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);

%plot frame camera in line-scan coordinate frame
scatter3(extLS(2,4), extLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%plot fov of line-scan on view-plane
plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);



s2 = subplot(1,3,2);
surf(Y, Z, radIntGP, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Const Mean GP Symmetric')
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

s3 = subplot(1,3,3);
surf(Y, Z, diffRad, 'EdgeColor', 'none'); hold on;
xlabel('y');
ylabel('z');
title('Absolute difference in GP')
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


%% OBSOLETE


% targetPntFrameHom = [targetPntFrame; ones(1, length(targetPntFrame(1,:)))];
% 
% targetPntLightSrcHom = poseLightSrc\targetPntFrameHom;
% targetPntLightSrc = targetPntLightSrcHom(1:3, :);
% 
% targetPntRadius = sqrt(sum(targetPntLightSrc.^2, 1));
% 
% srcDir = [0;0;1];
% 
% targetPntTheta = acos(sum(srcDir.*targetPntLightSrc,1)./targetPntRadius);
% 
% inOut = [targetPntRadius', targetPntTheta', targetL'];
% 
% GPDataDown = downsample(inOut, 50);
% 
% x = GPDataDown(:,1:2);
% y = GPDataDown(:,3);
% 
% meanfunc = [];                    % empty: don't use a mean function
% covfunc = @covSEiso;              % Squared Exponental covariance function
% likfunc = @likGauss;              % Gaussian likelihood
% hyp = struct('mean', [], 'cov', [0 0], 'lik', -1);
% 
% hyp2 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, x, y);
% 
% % nlml2 = gp(hyp2, @infGaussLik, [], covfunc, likfunc, x, y);
% 
% % size of plane in each dimension
% zDist = 1;
% 
% %Create mesh
% thetaTest = linspace(-pi/2, pi/2, 100);
% radiusTest = linspace(0, zDist, 100);
% [thetaX, radX] = meshgrid(thetaTest,radiusTest);
% 
% xMesh = [thetaX(:), radX(:)];
% 
% rows = size(thetaX);
% %reshape to mesh and translate
% 
% [mu, s2] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, x, y, xMesh);
% 
% intTest = reshape(mu,rows);
% 
% [xTest, yTest] = pol2cart(thetaX, radX);
% 
% 
% figGP = figure('Name', 'GP TEST');
% surf(xTest, yTest, intTest, 'EdgeColor', 'none'); hold on;
% xlabel('X');
% ylabel('Y');
% colorbar;
% 
% 
% % %%
% %
% %
% % % %% Non-linear optimisation using parameter light source model
% % %
% % % %Calculate incoming radiance recieved through the knowledge of the material
% % % targetL = (pi/target.Reflectance) .* targetInten;
% %
% % %The optimisation function
% % f = @(H)LightSourceLeastSqrObj(H, targetL', targetPntFrame, ...
% %     targetPntNormals, ligSourLoc, lightSrc.get_SourDirVec());
% %
% % Phi0 = [10, 1, 1];
% %
% % optOptions = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'SpecifyObjectiveGradient',true, 'CheckGradients', true, ...
% %     'MaxIterations', 1000000000, 'FunctionTolerance',1e-6, 'MaxFunctionEvaluations',1000000000, 'StepTolerance',1e-6);%,'ScaleProblem', ...
% % %'jacobian', 'InitDamping', 0.01, 'FiniteDifferenceType', 'central');
% % optOptions.Display = 'none';
% %
% % % optOptions = optimoptions('lsqnonlin', 'Algorithm', 'trust-region-reflective', 'SpecifyObjectiveGradient',false, 'CheckGradients', false, ...
% % %     'MaxIterations', 1000000000, 'FunctionTolerance',1e-6, 'MaxFunctionEvaluations',1000000000, 'StepTolerance',1e-15, 'OptimalityTolerance', 1e-6);
% % % optOptions.Display = 'none';
% %
% % %Perform optimisation
% % [optPhi,resNorm] = lsqnonlin(f,Phi0, [] , [], optOptions);
% %
%
%
%
%
% % %% Set up for training
% % xAngle = pi;
% % rotTarget = eul2rotm([0, 0, xAngle], 'ZYX');
% %
% predictorNames = {'theta', 'rho'};
% %
% % % Convert input to table
% inputTable = array2table(targetPolarLS', 'VariableNames', predictorNames);
% %
% predictors = inputTable(:, predictorNames);
% response = targetIntenNoise(:);
% isCategoricalPredictor = [false, false];
%
%
% % %% Train a linear regression model
% % % This code specifies all the model options and trains the model.
% % concatenatedPredictorsAndResponse = predictors;
% % concatenatedPredictorsAndResponse.targetInten = response;
% % linearModel = fitlm(concatenatedPredictorsAndResponse, 'linear', 'RobustOpts', 'off');
% % disp(linearModel);
% % % coeffLine = linearModel.Coefficients.Estimate;
% %
% % %% Train a quadratic regression model
% % % This code specifies all the model options and trains the model.
% % % concatenatedPredictorsAndResponse = predictors;
% % % concatenatedPredictorsAndResponse.targetInten = response;
% % QuadModel = fitlm(concatenatedPredictorsAndResponse, 'quadratic', 'RobustOpts', 'off');
% % disp(QuadModel);
% %
% % % coeffQuad = QuadModel.Coefficients.Estimate;
% %
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
% %
% % %% 2D plot of radiant intensity on view-plane
% %
% % %plot it as a rectangle starting at where the line-scan is positioned and
% % %only plot the YZ-plane
% % %size of plane in each dimension
% yDist = 0.5;
% zDist = 1;
% %
% % %Create mesh
% y = linspace(-yDist, yDist, 100);
% z = linspace(0, zDist, 100);
% x = 0;
% %
% [X,Y,Z] = meshgrid(x,y,z);
% %
% % %remove extra necessary singular dimension
% X = squeeze(X);
% Y = squeeze(Y);
% Z = squeeze(Z);
% %
% Yarr = Y(:);
% Zarr = Z(:);
% %
% [theta, rho] = cart2pol(Yarr, Zarr);
%
% thetaRhoMeshPred = [theta, rho];
% %
% rows = size(X);
% %
% % radIntMagLine = predict(linearModel, thetaRhoMeshPred);
% %
% % % a = coeffLine(1);
% % % b = coeffLine(2);
% % % c = coeffLine(3);
% % % radIntMagLine = a + b.*theta + c.*rho;
% % radIntMagLine = reshape(radIntMagLine,rows);
% %
% %
% %
% % radIntMagQuad = predict(QuadModel, thetaRhoMeshPred);
% %
% % % a = coeffQuad(1);
% % % b = coeffQuad(2);
% % % c = coeffQuad(3);
% % % d = coeffQuad(4);
% % % e = coeffQuad(5);
% % % f = coeffQuad(6);
% % % radIntMagQuad = a + b.*theta + c.*rho + d.*theta.*rho + e.*theta.^2 + f.*rho.^2;
% % radIntMagQuad = reshape(radIntMagQuad,rows);
% %
% radIntMagGP = predict(GPmatModel, thetaRhoMeshPred);
% radIntMagGP = reshape(radIntMagGP,rows);
% %
% %
% % %plot 2D view-plane surface with line-scan camera as origin. Only plotting
% % %YZ plane.
% % figLineReg = figure('Name', 'Radiant intensity on line-scan view-plane using LINEAR regression');
% % surf(Y, Z, radIntMagLine, 'EdgeColor', 'none'); hold on;
% % xlabel('y');
% % ylabel('z');
% % title('Linear regression model')
% % view(2);
% % scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled'); %line-scan origin
% % colormap(jet(100));
% % colorbar; caxis([0, maxRadiInt]);
% % axis equal;
% % %plot green light source
% % scatter3(lightSrcLS(2), lightSrcLS(3), (maxRadiInt+1), 200, [0,1,0], 'filled');
% % %plot frame camera
% % scatter3(frameCameraLS(2,4), frameCameraLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled');
% % %plot fov of line-scan on view-plane
% % plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
% %
% %
% %
% %
% % %plot 2D view-plane surface with line-scan camera as origin. Only plotting
% % %YZ plane.
% % figQuadReg = figure('Name', 'Radiant intensity on line-scan view-plane using QUADATRIC regression');
% % surf(Y, Z, radIntMagQuad, 'EdgeColor', 'none'); hold on;
% % xlabel('y');
% % ylabel('z');
% % title('Quadratic regression model')
% % view(2);
% % scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled'); %line-scan origin
% % colormap(jet(100));
% % colorbar; caxis([0, maxRadiInt]);
% % axis equal;
% % %plot green light source
% % scatter3(lightSrcLS(2), lightSrcLS(3), (maxRadiInt+1), 200, [0,1,0], 'filled');
% % %plot frame camera
% % scatter3(frameCameraLS(2,4), frameCameraLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled');
% % %plot fov of line-scan on view-plane
% % plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
% %
% %
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
% % %% Differences between groundtruth
% %
% % diffRadIntLine = abs(radIntMagLine - radIntMagDotTarget);
% % diffRadIntQuad = abs(radIntMagQuad - radIntMagDotTarget);
% % diffRadIntGP = abs(radIntMagGP - radIntMagDotTarget);
% %
% % figure('Name', 'abosolute difference of view-plane linear');
% % surf(Y, Z, diffRadIntLine, 'EdgeColor', 'none'); hold on;
% % xlabel('y');
% % ylabel('z');
% % title('absolute difference of view-plane linear')
% % view(2);
% % scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled'); %line-scan origin
% % colormap(flipud(hot));
% % colorbar; caxis([0, maxRadiInt]);
% % axis equal;
% % %plot green light source
% % scatter3(lightSrcLS(2), lightSrcLS(3), (maxRadiInt+1), 200, [0,1,0], 'filled');
% % %plot frame camera
% % scatter3(frameCameraLS(2,4), frameCameraLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled');
% % %plot fov of line-scan on view-plane
% % plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
% %
% %
% % figure('Name', 'abosolute difference of view-plane quadratic');
% % surf(Y, Z, diffRadIntQuad, 'EdgeColor', 'none'); hold on;
% % xlabel('y');
% % ylabel('z');
% % title('absolute difference of view-plane quadratic')
% % view(2);
% % scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled'); %line-scan origin
% % colormap(flipud(hot));
% % colorbar; caxis([0, maxRadiInt]);
% % axis equal;
% % %plot green light source
% % scatter3(lightSrcLS(2), lightSrcLS(3), (maxRadiInt+1), 200, [0,1,0], 'filled');
% % %plot frame camera
% % scatter3(frameCameraLS(2,4), frameCameraLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled');
% % %plot fov of line-scan on view-plane
% % plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
% %
% %
% % figure('Name', 'abosolute difference of view-plane GP');
% % surf(Y, Z, diffRadIntGP, 'EdgeColor', 'none'); hold on;
% % xlabel('y');
% % ylabel('z');
% % title('absolute difference of view-plane GP')
% % view(2);
% % scatter3(0, 0, (maxRadiInt+1), 200, [0,0,1], 'filled'); %line-scan origin
% % colormap(flipud(hot));
% % colorbar; caxis([0, maxRadiInt]);
% % axis equal;
% % %plot green light source
% % scatter3(lightSrcLS(2), lightSrcLS(3), (maxRadiInt+1), 200, [0,1,0], 'filled');
% % %plot frame camera
% % scatter3(frameCameraLS(2,4), frameCameraLS(3,4), (maxRadiInt+1), 200, [1,0,0], 'filled');
% % %plot fov of line-scan on view-plane
% % plot3(yPoly, zPoly, (maxRadiInt+1)*ones(size(yPoly)), 'Color', [1,1,1], 'LineWidth', 2);
%
% %% Call_back function
%
% % function  distFromSrc_callback(src, ~, lightSrc, fig, S)
% % lightSrc.set_distFromSrc = src.Value;
% % lightSrc.PlotRadiantSlice(fig, S);
% % end
%
% % function  distFromSrc_callback(src, ~, lightSrc, fig, S, mu)
% % lightSrc.setMu(mu);
% % lightSrc.PlotRadiantSlice(fig, S, src.Value);
% % end


function vp_callback(src, ~, vpPatch)
value = get(src,'Value');

if value
    vpPatch.Visible = true;
else
    vpPatch.Visible = false;
end

end