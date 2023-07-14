% Build models of the light intensity field for a real light source given
% measurements of a planar diffuse target captured by a line-scan frame
% camera system. This will train a light-source-mean GP and a parametric
% least squares models that assume the source is a non-isotropic near-field
% disk. The location and principal direction of the light source needs to
% be determined prior to building the models

%Author: Jasprabhjit Mehami, 13446277
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

%check if mex_ChArUco_Pose has been built
if ~exist(fullfile("ext_lib", "mex_ChArUco_Pose", "bin", "ArucoPixDect.mexa64"), 'file')
    warning("Cannot estimate board pose. Please build mex_ChArUco_Pose submodule")
else
    addpath(genpath(fullfile("ext_lib", "mex_ChArUco_Pose")));
end

%code for this project
addpath('src');

%parameter file
paramFile = fullfile('parameter_files', 'light_modelling.yaml');
if ~exist(paramFile, 'file')
    error("YAML parameter file not found");
end


%% Read YAML file containing the pattern specifications and parameters for code
% All dimensions are in metres

paramYaml = yaml.ReadYaml(paramFile);
displayOn = paramYaml.DisplayOn;
frameCamera = lower(paramYaml.FrameCamera);
usingTestData = paramYaml.UseTestData;
curDownSamplingGP = paramYaml.DownSamplingGP;

%maximum pixel intensity for uint14;
maxIntUint14 = 2^14 - 1;

%STD in the radiance measurements
stdRadianceNoisePer = paramYaml.sigmaRadianceNoisePercent;
cur_stdRadianceNoise = (stdRadianceNoisePer/100);

%% Directories and files

disp('Checking directories...')

commonDir = "/home/jmeh/Experimental_Data/lamb_carcass_une_2022/common_data/";
%frame camera intrinsic parameters file
% frameIntrFile = fullfile('frame_camera_intrinsic', [frameCamera, '.mat']);
frameIntrFile = fullfile(commonDir, "frame_camera_intrinsic.mat");

if ~exist(frameIntrFile, 'file')
    error("Frame camera intrinsic parameters not found");
end

if usingTestData
    sourDir = fullfile('test_data', 'light_map');
else
    %Get source directory where images are located and results will be saved
    sourDir = uigetdir(fullfile('~'),'Provide source directory where light modelling images are located?');
end

%frame image directory
frameDir = fullfile(sourDir, 'Frame'); %Directory containing images
if ~exist(frameDir, 'dir')
    error('Source directory does not contain directory called "Frame" which should contain RGB images');
end

fullPathFrame = fullfile(frameDir, '*.png');

%Need to get number of images in directory
numImagesFrame = numel(dir(fullPathFrame));

if numImagesFrame < 1
    error('no images in provided image directory')
end

%hyperspectral line-scan image directory
hsDir = fullfile(sourDir, 'Line-scan'); %Directory containing images
if ~exist(hsDir, 'dir')
    error('Source directory does not contain directory called "Line-scan" which should contain hyperspectral images');
end

fullPathHS = fullfile(hsDir, '*.png');

%Need to get number of images in directory
numImagesHS = numel(dir(fullPathHS));

if numImagesHS < 1
    error('no images in provided image directory')
end

if numImagesHS ~= numImagesFrame
    error("number of images in folders are not the same");
end

numImages = numImagesHS;

%dark reference image
darkRefImageFile = fullfile(sourDir, "dark_ref.png");
if ~exist(darkRefImageFile, 'file')
    error("Dark reference images not found");
end

%light triangulation results
% ligTrigFile = fullfile(sourDir, 'pt_light.mat');
ligTrigFile = fullfile(commonDir, 'pt_light.mat');
 
if ~exist(ligTrigFile, 'file')
    error("Could not find point light triangulation in the following path:\n %s\n", ligTrigFile);
end

%camera system calibration results
% cameraSysCaliFile = fullfile(sourDir, 'camera_system_optimised_parameters.mat');
cameraSysCaliFile = fullfile(commonDir, 'camera_system_optimised_parameters.mat');
if ~exist(cameraSysCaliFile, 'file')
    error("Could not find camera system calibration in the following path:\n %s\n", cameraSysCaliFile);
end

%white stripe reflectance CSV file
whiteStripeReflFile = fullfile('parameter_files', 'white_stripe_reflectance.csv');
if ~exist(whiteStripeReflFile, 'file')
    error("white stripe reflectance CSV not found");
end

%result directory
resultsDir = fullfile(sourDir, 'Result');
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

%% Load intrinsic parameters and board parameters

disp('Loading frame intrinsic and board parameters...');

%load frame camera intrinsic parameter file
franeIntrStruct = load(frameIntrFile);

cameraParamF = franeIntrStruct.cameraParams;

%extract focal length components in pixels
fx = cameraParamF.FocalLength(1);
fy = cameraParamF.FocalLength(2);

%optical centre
u0 = cameraParamF.PrincipalPoint(1);
v0 = cameraParamF.PrincipalPoint(2);

%array of intrinsic parameters that introduce noise (assume noise in
%intrinsic distortion is zero)
thetaFrameintr = [fx, fy, u0, v0];

%Size of the images from the RGB camera
frameImgSize = cameraParamF.ImageSize;

%distortion parameters
distRad = cameraParamF.RadialDistortion;
distTan = cameraParamF.TangentialDistortion;
% distCoefCV = [distRad(1:2), distTan, distRad(3)]; %array of distortion coefficients in opencv format
distCoefCV = zeros(1,5); %array of distortion coefficients in opencv format

%extract the error in the intrinsic parameters of the frame camera
% std_f = estimationErrors.IntrinsicsErrors.FocalLengthError;
% std_u0v0 = estimationErrors.IntrinsicsErrors.PrincipalPointError;
% stdRGB_Intr = [std_f,std_u0v0];

%ArUco board parameters
xNumMarker = paramYaml.NumCols;
yNumMarker = paramYaml.NumRows;
arucoLen = paramYaml.ArucoSideLength;
sepLen = paramYaml.SeparationLength;
numMarkersExp = paramYaml.NumberExpectedMarker;

%intrinsic object for the RGB camera
frameIntrinsic = cameraIntrinsics(thetaFrameintr(1:2),thetaFrameintr(3:4), frameImgSize);
kFrame = frameIntrinsic.IntrinsicMatrix';

%% Load calibration parameters for the line-scan frame camera system

disp('Loading line-scan calibration parameters ...');

%calibration data
load(cameraSysCaliFile)

%grab the last optimised parameters
caliParam = linescan_Opt_Param(end,:);

%extract individual parameters
K1 = caliParam(9);
K2 = caliParam(10);
K3 = 0;
P1 = 0;
P2 = caliParam(11);
fy = caliParam(7);
v0 = caliParam(8);
t = caliParam(1:3);
rotEul = caliParam(4:6);
rotMat = eul2rotm(rotEul, 'ZYX');

% Transformation of the frame w.r.t to line-scan camera
T_F_2_LS = [rotMat, t'; 0, 0, 0, 1 ];

% Transformation of the line-scan w.r.t to frame camera
T_LS_2_F = T_F_2_LS\eye(4);

%dimensions of line-scan image
rowsLS = 320;
colsLS = 1;

%intrinsic object for the linescan camera
intrLS = cameraIntrinsics([1, fy], [realmin,v0],[rowsLS,colsLS], 'RadialDistortion', [K1,K2,K3], 'TangentialDistortion', [P1,P2]);
cameraParamLS = cameraParameters('IntrinsicMatrix', intrLS.IntrinsicMatrix, 'ImageSize', [rowsLS,colsLS], 'RadialDistortion', [K1,K2,K3], 'TangentialDistortion', [P1,P2]);

%intrinsic matrix of line-scan camera
kLS = intrLS.IntrinsicMatrix';
kLS(1,3) = 0;

%% Load images for both cameras

disp('Loading images...');

%Preallocate space for cell arrays of images
imagesLS = cell(1,numImages);
imagesFrame = cell(1,numImages);

% Load all images
for i = 1:numImages
    imagesFrame{i} = undistortImage(imread(fullfile(frameDir, ['img', num2str(i),'.png'])), cameraParamF);
    imagesLS{i} = imread(fullfile(hsDir, ['hs', num2str(i),'.png']));
end

imageDark = imread(darkRefImageFile);


%% Simulator figure
close all;

%figure for light simulator
figSim = figure('Name', 'Light Simulator with Frame Camera Origin');
%plot frame camera at origin
plotCamera('Location', zeros(1,3), 'Orientation', eye(3), 'Size', 0.05, 'AxesVisible', true); hold on;

xlabel('x');
ylabel('y');
zlabel('z');
grid on;
axis equal;

%plot line-scan camera in simulator figure w.r.t to frame camera
figure(figSim);
plotCamera('Location', T_LS_2_F(1:3,4), 'Orientation', T_LS_2_F(1:3, 1:3)', 'Size', 0.05, 'AxesVisible', true, 'Color', [0,0,1]); hold on;

%% Pose of diffuse ArUco board from each image

disp('Estimating extrinisic pose of each image...');

extMatFile = fullfile(resultsDir, 'imgExtrinsicPose.mat');

%display 3D pose of board w.r.t to frame camera
if displayOn
    %edge pts of pattern on board in homogeneous coordinates
    patEdgePtsHom = [
        0,0,0, 1;
        0, yNumMarker*arucoLen + (yNumMarker - 1)*arucoLen, 0, 1;
        xNumMarker*arucoLen + (xNumMarker - 1)*arucoLen, yNumMarker*arucoLen + (yNumMarker - 1)*arucoLen, 0, 1;
        xNumMarker*arucoLen + (xNumMarker - 1)*arucoLen, 0, 0, 1;
        0,0,0, 1;
        ]';

    global fillHandleCell;
    fillHandleCell = cell(1,numImages);

    %checkerbox to turn surface plane visibility on/off
    boardDispCheck = uicontrol('Parent',figSim,'Style','checkbox', 'String', 'Show Planes', 'Position', [20,45,200,20] );
    boardDispCheck.Callback = @(src, eventData) boardDispCheck_callback(src, eventData, fillHandleCell);
end

if ~exist(extMatFile, 'file')
    figFrameImg = figure('Name','ArUco detected markers');
    %store all the poses of each found pattern
    extPosePattern = zeros(4,4,numImages);
    %used to filter images where the pose can't be found
    goodImages = zeros(1,numImages);
    numGoodImg = 0;

    % 2D marker corner positions in world coordinates (metres)
    markerCornerCell = ArUcoBoardMarkerCornersCell(0, xNumMarker, yNumMarker, arucoLen, sepLen);

    %loop through each image and get pose of board using ArUco pattern
    for imgLoop = 1:numImages

        [rotMat, trans, found, imgDisp] = ArucoPosEst(imagesFrame{imgLoop}, markerCornerCell, cameraParamF, false);

        if ~found
            continue;
        end

        %image is good
        numGoodImg = numGoodImg + 1;
        goodImages(numGoodImg) = imgLoop;

        %store found extrinsic parameter
        extPosePattern(:,:,numGoodImg) = [rotMat,trans'; 0, 0, 0, 1];

        %display the frame camera image with detected markers and 3D pose of
        %the board in simulator figure
        if displayOn
            figure(figFrameImg);
            clf(figFrameImg);
            imshow(imgDisp); hold on;
            drawnow();

            %edge points tr'rg  nsformed into the frame camera c.f
            patEdgePtsFrameHom = extPosePattern(:,:,numGoodImg)*patEdgePtsHom;
            figure(figSim);
            fillboard = fill3(patEdgePtsFrameHom(1,:),patEdgePtsFrameHom(2,:),patEdgePtsFrameHom(3,:),'c', 'FaceAlpha', 0.1);
            fillboard.HandleVisibility = 'callback';
            fillboard.Visible = 'off';
            fillHandleCell(imgLoop) = {fillboard};
            drawnow();
        end
    end

    goodImages = goodImages(1:numGoodImg);
    extPosePattern = extPosePattern(:,:,1:numGoodImg);
    save(extMatFile, 'extPosePattern', 'goodImages', 'numGoodImg')
else
    fprintf("\tLoading existing extrinsic parameters\n");
    load(extMatFile);

    if displayOn

        %loop through each image and get pose of board using ArUco pattern
        for imgLoop = 1:numImages
            %edge points tr'rg  nsformed into the frame camera c.f
            patEdgePtsFrameHom = extPosePattern(:,:,imgLoop)*patEdgePtsHom;
            figure(figSim);
            fillboard = fill3(patEdgePtsFrameHom(1,:),patEdgePtsFrameHom(2,:),patEdgePtsFrameHom(3,:),'c', 'FaceAlpha', 0.1);
            fillboard.HandleVisibility = 'callback';
            fillboard.Visible = 'off';
            fillHandleCell(imgLoop) = {fillboard};
            drawnow();
        end
    end
end

%remove all data from the frame images where we could not find proper
%extrinsic parameters
imagesFrame = imagesFrame(goodImages);
imagesLS = imagesLS(goodImages);
numImages = numGoodImg;

%% Get normalised line-scan image coordinates and their radiance

disp('Calculating normalised image coordinates and extracting radiance...');

% relevant pixel locations in line-scan image
bandStart = 40;
bandEnd = 203;
hypBands = bandStart:bandEnd;
numBands = length(hypBands);

vImgLineDist = 1:rowsLS; %distorted v component pixel positions on image line
imgLineDist = [zeros(size(vImgLineDist));vImgLineDist]'; %distorted pixels on image line (u = 0)
imgLine = undistortPoints(imgLineDist, cameraParamLS); %undistort image line

%determine normalised pixel coordinates
imgLineHom = [imgLine'; ones(1,length(imgLineDist(:,1)))];
normPtImgLS = (kLS\imgLineHom)';

RadiImgCell = cell(1,numImages);

for imgLoop = 1:numImages
    curImg = imagesLS{imgLoop};
    %subtract dark reference
    radiImg = curImg - imageDark;
    radiImg = double(radiImg(hypBands, :))./maxIntUint14;
    RadiImgCell(imgLoop) = {radiImg};
end

%% Determine 3D location of pixels on planar target diffuse board

disp('Calculating 3D location of pixels on planar target board...');

% target is offset from the board by its thickness
T_target_2_board = trvec2tform([0,0,2/1000]);

ptsOnTarget = zeros(50000,3);
ptsSurfNorm = zeros(50000,3);
ptsRadianceTarget = zeros(50000,numBands);

numPtsCount = 0;

if displayOn
    figProjectLineFrame = figure('Name', 'Projected hyperspectral image-line to frame camera');
    bandDisplay = 100; %band used for displaying hyperspectral line
end

%Get the 3D locations of the line-scan pixel points in the frame camera
%coordinate frame
for imgLoop = 1:numImages
    curNormRadiImg = RadiImgCell{imgLoop};

    %pose of target w.r.t to frame camera c.f
    T_tar_2_F = extPosePattern(:,:,imgLoop)*T_target_2_board;

    %pose of target w.r.t to line-scan camera c.f
    T_tar_2_LS = T_F_2_LS*T_tar_2_F;

    %origin taken to be point on target surface
    ptOnTarget = tform2trvec(T_tar_2_LS);

    %surface normal of target which is the z-direction of the board's
    %coordinate frame
    tarSurfNorm = T_tar_2_LS(1:3,3)';

    numPts = size(normPtImgLS, 1);

    %display frame image
    if displayOn
        figure(figProjectLineFrame);
        clf(figProjectLineFrame);
        imshow(imagesFrame{imgLoop}); hold on;
        targetPtsFrameRadiance = zeros(numPts, 4);
        dispPtCount = 0;
    end

    %loop through points for current image and save 3D point location,
    %surface normal at that point, and band reflectances
    for ptLoop = 1:numPts
        normPtLS = normPtImgLS(ptLoop,:);

        %get intersection of pixel ray to target plane
        [pntTargetLS, validIntersection] =  CameraPixelRayIntersectPlane(normPtLS, tarSurfNorm, ptOnTarget);

        %point does not intersect plane
        if ~validIntersection
            continue;
        end

        %transform intersected point to frame camera c.f
        pntTargetLSHom = [pntTargetLS, 1]';
        pntTargetFrameHom = T_LS_2_F*pntTargetLSHom;

        numPtsCount = numPtsCount + 1;

        ptsRadianceTarget(numPtsCount, :) = curNormRadiImg(:,ptLoop)';

        %store surface normal and 3D point in frame camera c.f
        ptsSurfNorm(numPtsCount,:) = T_tar_2_F(1:3,3)';
        ptsOnTarget(numPtsCount,:) = pntTargetFrameHom(1:3);

        if displayOn
            dispPtCount = dispPtCount + 1;
            targetPtsFrameRadiance(dispPtCount, :) = [ptsOnTarget(numPtsCount,:), ptsRadianceTarget(numPtsCount, bandDisplay)];
        end
    end

    %display projected line onto image
    if displayOn
        %project 3D line points and radiance to 2D image pixels
        lineImgPts = projectPoints(targetPtsFrameRadiance, kFrame,eye(4), [], frameImgSize);


        x = lineImgPts(:,1)';
        y = lineImgPts(:,2)';
        z = zeros(size(x));

        %normalise grayscale radiance values and find corresponding rgb
        %from color map
        col = normalize(lineImgPts(:,3)', 'range');
        colRGB = zeros(1,length(col),1);
        colRGB(:,:,1) = col;

        %plot the projected hyperspectral line to the frame camera
        %coordinate frame.
        plNaive = surface([x;x],[y;y], [z;z], [colRGB;colRGB],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',2, 'DisplayName', 'Active', 'FaceColor', [0.9290 0.6940 0.1250]);
        colormap(jet)
    end
end

%clip to size
ptsOnTarget = ptsOnTarget(1:numPtsCount, :);
ptsSurfNorm = ptsSurfNorm(1:numPtsCount, :);
ptsRadianceTarget = ptsRadianceTarget(1:numPtsCount, :);

%% Light source triangulation parameters

%extract location and principal direction.
load(ligTrigFile);
T_S_2_F = [rotLightSrc, locLightSrc; 0,0,0,1];

%light source simulator for visualisation
lightSrc = LightSimulator(locLightSrc, rotLightSrc, figSim, ptsOnTarget, ptsRadianceTarget, numBands, numPtsCount);

%% Diffuse reflectance target stripe

%read the reflectance values for each band of the hyperspectral camera
%stored in file
data = readmatrix(whiteStripeReflFile);
targetReflectances = data(:,2);

%% Prepare Data for building models

ptsFrameHom = [ptsOnTarget'; ones(1,length(ptsOnTarget(:,1)))];
%Transform points to LS C.F
ptsLSHom = T_F_2_LS*ptsFrameHom;

%Create dense mesh of points using the edges of the measured line-scan
%points
x = mean(ptsLSHom(1,:));
rangeExpandRatio = 0.25;

%expand the y and z limits
yzMin = [min(ptsLSHom(2,:)), min(ptsLSHom(3,:))];
yzMin = yzMin - abs(rangeExpandRatio.*yzMin);
yzMax = [max(ptsLSHom(2,:)), max(ptsLSHom(3,:))];
yzMax = yzMax + abs(rangeExpandRatio.*yzMax);
y = linspace(yzMin(1), yzMax(1), 200);
z = linspace(yzMin(2), yzMax(2), 200);
[X,Y,Z] = meshgrid(x,y,z);

%remove extra unnecessary singular dimension
X = squeeze(X);
Y = squeeze(Y);
Z = squeeze(Z);
rows = size(X);

%These points are in the line-scan c.f, transform them into
%the frame camera c.f (world coordinates)
ptsTestLS = [X(:),Y(:),Z(:)]';
ptsTestLSHom = [ptsTestLS; ones(1, size(ptsTestLS, 2))];
ptsTestFrameHom = T_LS_2_F*ptsTestLSHom;

%transform points from frame camera C.F to light source C.F
ptsTestLightSrcHom = T_S_2_F\ptsTestFrameHom;
ptsTestLightSrc = ptsTestLightSrcHom(1:3, :);

%Training points without radiant intensity
[ptsRadius,ptsTheta] = cart2sphZ(ptsTestLightSrc(1,:), ptsTestLightSrc(2,:), ptsTestLightSrc(3,:));

testingX = [ptsRadius', ptsTheta'];

%Reshape the testing points in the frame camera C.F. into a mesh
XTestFr = reshape(ptsTestFrameHom(1,:),rows);
YTestFr = reshape(ptsTestFrameHom(2,:),rows);
ZTestFr = reshape(ptsTestFrameHom(3,:),rows);

%% Least Squares models

disp('Estimating light source model using least squares...');


%flag to run the GP training
runLS = true;
savedOptmFile = fullfile(resultsDir, 'ls_lightsrc_optm.mat');

%check if GP was previously optimised
if exist(savedOptmFile, 'file')
    disp('Loading existing least squares light source mean model');
    load(savedOptmFile);
    runLS = false;
end

radIntVP_LeasSqr = zeros([rows, numBands]);

if runLS
    optOptions = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'SpecifyObjectiveGradient',true, 'CheckGradients', false, ...
        'MaxIterations', 1000000000, 'FunctionTolerance',1e-6, 'MaxFunctionEvaluations',1000000000, 'StepTolerance',1e-8, ...
        'FiniteDifferenceType', 'central', 'ScaleProblem','none');
    optOptions.Display = 'none';

    optPhiBand = zeros(numBands, 3);
    resBand = zeros(numBands, 1);
    Phi0 = [1,1,1];

    for bandLoop = 1:numBands
        targetL = ptsRadianceTarget(:, bandLoop);
        [optPhi, res] = LightSrcOptmLS(lightSrc, Phi0, targetL', ptsOnTarget', ptsSurfNorm', targetReflectances(bandLoop), optOptions, cur_stdRadianceNoise);
        resBand(bandLoop) = res;
        %store optimised parameters and residual
        optPhiBand(bandLoop, :) = optPhi;
        %create temporary light source using current optimised parameters
        tempLigSrcLS = LightSimulator(locLightSrc, rotLightSrc, optPhi(1), optPhi(2), optPhi(3));
        radIntMag = tempLigSrcLS.RadiantIntensityMesh(XTestFr, YTestFr, ZTestFr);
        radIntVP_LeasSqr(:,:,bandLoop) = radIntMag;
    end

    save(savedOptmFile, 'optPhiBand', 'resBand')

else
    %use trained parameters and get irradiance
    for bandLoop = 1:numBands
        optPhi = optPhiBand(bandLoop, :);
        %create temporary light source using current optimised parameters
        tempLigSrcLS = LightSimulator(locLightSrc, rotLightSrc, optPhi(1), optPhi(2), optPhi(3));
        radIntMag = tempLigSrcLS.RadiantIntensityMesh(XTestFr, YTestFr, ZTestFr);
        radIntVP_LeasSqr(:,:,bandLoop) = radIntMag;
    end

end

%% Training GP model

%calculate the direction light vector (point to light source)
pntLightVec = locLightSrc - ptsOnTarget';
dirPntLightVec = pntLightVec./vecnorm(pntLightVec);

%transform points from frame camera C.F to light source C.F
ptsTestLightSrcHom = T_S_2_F\ptsFrameHom;
ptsLightSrc = ptsTestLightSrcHom(1:3, :);

%Training points without radiant intensity
[ptsRadius,ptsTheta] = cart2sphZ(ptsLightSrc(1,:), ptsLightSrc(2,:), ptsLightSrc(3,:));

trainingXY = [ptsRadius', ptsTheta'];

%irradiance which is actually measured by line-scan camera
radIntVP_Meas = zeros(size(ptsRadianceTarget));

dotProductVP = zeros(size(ptsRadianceTarget));

%Calculate radiant intensity magnitude used for building model
for bandLoop = 1:numBands
    curRefl = targetReflectances(bandLoop);
    targetL = ptsRadianceTarget(:, bandLoop);
    radIntMagPnt = (targetL'.*pi)./(curRefl.*dot(ptsSurfNorm', dirPntLightVec,1));
    radIntVP_Meas(:,bandLoop) = radIntMagPnt';
    dotProductVP(:,bandLoop) = dot(ptsSurfNorm', dirPntLightVec,1);
end

%store the radiant intensity from GP model with corresponding variance
radIntVP_GP = zeros([rows, numBands]);
varVP_GP = zeros([rows, numBands]);

%flag to run the GP training
runTraining = true;

savedOptmFile = fullfile(resultsDir, 'gp_lightsrc_optm.mat');

%check if GP was previously optimised
if exist(savedOptmFile, 'file')
    disp('Loading existing trained GP light source mean model');
    load(savedOptmFile);

    %if downsampling and STD radiance noise have not changed, GP training is not necessary.
    if (curDownSamplingGP == downSamplingGP) && (stdRadianceNoise == cur_stdRadianceNoise)
        runTraining = false;
    end
end


if runTraining
    %cell arrays to store hyperparameters and training data for each band
    hypOptCell = cell(numBands, 1);
    downSamplingGP = curDownSamplingGP; %set current downSamplingGP
    stdRadianceNoise = cur_stdRadianceNoise; % set current radiance noise

    for bandLoop = 1:numBands
        trainingXY(:,3) = radIntVP_Meas(:, bandLoop);
        [mu, varMu, curHypOptStruct] = LightSrcOptmGP(testingX, trainingXY, 3, stdRadianceNoise, downSamplingGP, 10000, false);
        radIntVP_GP(:,:,bandLoop) = reshape(mu,rows);
        varVP_GP(:,:,bandLoop) = reshape(varMu,rows);
        hypOptCell(bandLoop) = {curHypOptStruct};
        fprintf('\tGP optimised and queried band %i of %i\n', bandLoop, numBands);
    end

    save(savedOptmFile, 'hypOptCell', 'downSamplingGP', 'stdRadianceNoise');
else
    for bandLoop = 1:numBands
        curHypOptStruct = hypOptCell{bandLoop};
        [mu, varMu] = LightSrcOptmGP(testingX, curHypOptStruct);
        radIntVP_GP(:,:,bandLoop) = reshape(mu,rows);
        varVP_GP(:,:,bandLoop) = reshape(varMu,rows);
        fprintf('\tGP queried band %i of %i\n', bandLoop, numBands);
    end
end

%% Plot the line-scan view-plane to compare results

cab(figSim);

figRadianceVRadiantIntensity_VP = figure('Name', 'Line-scan View-plane');

%%%%%%%%%%%%%%%% Measured Radiance on view-plane
sRadiance = subplot(1,3,1);
scatRadiance = scatter3(ptsLSHom(2,:), ptsLSHom(3,:), ptsRadianceTarget(:,1), 10, ptsRadianceTarget(:,1), 'filled');
hold on;

% maxInt = max(intenTarget(:,1));
% scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), maxInt, 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
% scatter3(extLS(2,4), extLS(3,4), maxInt, 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
% scatter3(0, 0, maxInt, 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin

xlabel('y');
ylabel('z');
% title('Measured Radiance')
view(2);
colormap(sRadiance, hot(10000));
% axis equal;
colorbar();


%%%%%%%%%%%%%%%% Measured Radiance on view-plane
sDot = subplot(1,3,2);
scatDot = scatter3(ptsLSHom(2,:), ptsLSHom(3,:), dotProductVP(:,1), 10, dotProductVP(:,1), 'filled');
hold on;

% maxInt = max(intenTarget(:,1));
% scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), maxInt, 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
% scatter3(extLS(2,4), extLS(3,4), maxInt, 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
% scatter3(0, 0, maxInt, 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin

xlabel('y');
ylabel('z');
% title('Dot Product')
view(2);
colormap(sDot, parula(10000));
% axis equal;
colorbar();



%%%%%%%%%%%%%%%% Estimated Radiant Intensity on view-plane
sMeas = subplot(1,3,3);
scatMeas = scatter3(ptsLSHom(2,:), ptsLSHom(3,:), radIntVP_Meas(:,1), 10, radIntVP_Meas(:,1), 'filled');
hold on;

% maxInt = max(radIntVP_Meas(:,1));
% scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), maxInt, 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
% scatter3(extLS(2,4), extLS(3,4), maxInt, 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
% scatter3(0, 0, maxInt, 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin

xlabel('y');
ylabel('z');
% title('Estimated Radiant Intensity')
view(2);
colormap(sMeas, hot(10000));
% axis equal;
colorbar();




sliderFig1 = figure('Name', 'Band slider', 'Position', [659,832,580,65]);

slider_band = uicontrol('Parent',sliderFig1,'Style','slider','Position',[150,40,400,40],...
    'value',1, 'min', 1, 'max',numBands, 'SliderStep', [1/(numBands+1), 1/(numBands+1)]);
bgcolor = sliderFig1.Color;
uicontrol('Parent',sliderFig1,'Style','text','Position',[150,10,23,23],...
    'String','1','BackgroundColor',bgcolor);
uicontrol('Parent',sliderFig1,'Style','text','Position',[550,10,40,23],...
    'String',num2str(numBands),'BackgroundColor',bgcolor);
uicontrol('Parent',sliderFig1,'Style','text','Position',[300,10,200,23],...
    'String','Band','BackgroundColor',bgcolor);
% callback function at the end of the script
slider_band.Callback = @(src, eventData) band_callback2(src, eventData, scatRadiance, scatDot, scatMeas, ptsRadianceTarget, dotProductVP, radIntVP_Meas);


figLS_VP = figure('Name', 'Line-scan View-plane models');





%%%%%%%%%%%%%%%% measured view-plane
sMeas = subplot(1,3,1);
scatMeas = scatter3(ptsLSHom ...
    (2,:), ptsLSHom(3,:), radIntVP_Meas(:,1), 10, radIntVP_Meas(:,1), 'filled');
hold on;

maxInt = max(radIntVP_Meas(:,1));
% scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), maxInt, 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
% scatter3(extLS(2,4), extLS(3,4), maxInt, 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
% scatter3(0, 0, maxInt, 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin

xlabel('y');
ylabel('z');
% title('Measured')
view(2);
colormap(sMeas, hot(10000));
% axis equal;
xlim([-0.12, 0.12]);
ylim([0.31, 0.57]);
grid off;



%%%%%%%%%%%%%%%% least squares view-plane
sLeasSqr = subplot(1,3,2);
surfLeasSqr = surf(Y, Z, radIntVP_LeasSqr(:,:,1), 'EdgeColor', 'none'); hold on;

maxInt = max([ maxInt, max(radIntVP_LeasSqr(:,:,1), [], 'all')]);


% scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), maxInt, 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
% scatter3(extLS(2,4), extLS(3,4), maxInt, 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
% scatter3(0, 0, maxInt, 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin

xlabel('y');
ylabel('z');
% title('Least Squares')
view(2);
colormap(sLeasSqr, hot(10000));
% axis equal;
xlim([-0.12, 0.12]);
ylim([0.31, 0.57]);
grid off;

%%%%%%%%%%%%%%%% least squares view-plane
sGP = subplot(1,3,3);
surfGP = surf(Y, Z, radIntVP_GP(:,:,1), 'EdgeColor', 'none'); hold on;

maxInt = max([ maxInt, max(radIntVP_GP(:,:,1), [], 'all')]);


% scatter3(poseLightSrcLS(2,4), poseLightSrcLS(3,4), maxInt, 200, [0,1,0], 'filled','MarkerEdgeColor', [0,0,0]);
% scatter3(extLS(2,4), extLS(3,4), maxInt, 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);
% scatter3(0, 0, maxInt, 200, [0,0,1], 'filled', 'MarkerEdgeColor', [0,0,0]); %line-scan origin

xlabel('y');
ylabel('z');
% title('Gaussian Process')
view(2);
colormap(sGP, hot(10000));
% axis equal;
xlim([-0.12, 0.12]);
ylim([0.31, 0.57]);
colorbar();

grid off;

caxis([0, maxInt]);

sliderFig = figure('Name', 'Band slider', 'Position', [659,832,580,65]);

slider_band = uicontrol('Parent',sliderFig,'Style','slider','Position',[150,40,400,40],...
    'value',1, 'min', 1, 'max',numBands, 'SliderStep', [1/(numBands+1), 1/(numBands+1)]);
bgcolor = sliderFig.Color;
uicontrol('Parent',sliderFig,'Style','text','Position',[150,10,23,23],...
    'String','1','BackgroundColor',bgcolor);
uicontrol('Parent',sliderFig,'Style','text','Position',[550,10,40,23],...
    'String',num2str(numBands),'BackgroundColor',bgcolor);
uicontrol('Parent',sliderFig,'Style','text','Position',[300,10,200,23],...
    'String','Band','BackgroundColor',bgcolor);
% callback function at the end of the script
slider_band.Callback = @(src, eventData) band_callback(src, eventData, scatMeas, surfLeasSqr, surfGP, radIntVP_Meas, radIntVP_LeasSqr, radIntVP_GP, sMeas, sLeasSqr, sGP);

%%

%%%%%%%%%%%%%%%% measured view-plane
bandFig = 13;
maxInt = max(radIntVP_Meas(:,bandFig),[], 'all');
% maxInt = max([ maxInt, max(radIntVP_LeasSqr(:,:,bandFig), [], 'all')]);
% maxInt = max([ maxInt, max(radIntVP_GP(:,:,bandFig), [], 'all')]);

hfig = figure();
scatter3(ptsLSHom(2,:), ptsLSHom(3,:), radIntVP_Meas(:,bandFig), 10, radIntVP_Meas(:,bandFig), 'filled');
hold on;

xlabel('y (m)');
ylabel('z (m)');
view(2);
colormap(hot(10000));
xlim([-0.12, 0.12]);
ylim([0.3, 0.57]);
grid off;
colorbar();
caxis([0, maxInt]);

a1 = gca;

imgName = fullfile(resultsDir, 'exp_mes.png');
exportpropergraphic(hfig, imgName,  1.5);


%%%%%%%%%%%%%%%% least squares view-plane

hfig = figure();
surf(Y, Z, radIntVP_LeasSqr(:,:,bandFig), 'EdgeColor', 'none'); hold on;
xlabel('y (m)');
ylabel('z (m)');
view(2);
colormap(hot(10000));
xlim([-0.12, 0.12]);
ylim([0.3, 0.57]);
grid off;
colorbar();
caxis([0, maxInt]);

a2 = gca;

imgName = fullfile(resultsDir, 'exp_lsq.png');
exportpropergraphic(hfig, imgName,  1.5);


%%%%%%%%%%%%%%%% least squares view-plane

hfig = figure();
surf(Y, Z, radIntVP_GP(:,:,bandFig), 'EdgeColor', 'none'); hold on;
xlabel('y (m)');
ylabel('z (m)');
view(2);
colormap(hot(10000));
xlim([-0.12, 0.12]);
ylim([0.3, 0.57]);
grid off;
colorbar();
caxis([0, maxInt]);

a3 = gca;

linkaxes([a1,a2,a3], 'xy');

imgName = fullfile(resultsDir, 'exp_gp.png');
exportpropergraphic(hfig, imgName,  1.5);

%% Plot GP irradiance planes on simulator figure

x = linspace(-0.2, 0.2, 100);
y = 0;
z = linspace(0, 0.5, 100);
[X,Y,Z] = meshgrid(x,y,z);

%remove extra unnecessary singular dimension
X = squeeze(X);
Y = squeeze(Y);
Z = squeeze(Z);
rows = size(X);

%These points are in the line-scan c.f, transform them into
%the frame camera c.f (world coordinates)
ptsTestLightSrc = [X(:),Y(:),Z(:)]';
ptsTestLightSrc = [ptsTestLightSrc; ones(1, size(ptsTestLightSrc, 2))];
ptsTestFrameHom = T_S_2_F*ptsTestLightSrc;


%Training points without radiant intensity
[ptsRadius,ptsTheta] = cart2sphZ(ptsTestLightSrc(1,:), ptsTestLightSrc(2,:), ptsTestLightSrc(3,:));
testingX = [ptsRadius', ptsTheta'];

irrVisSim = zeros([size(X), numBands]);
for bandLoop = 1:numBands
    curHypOptStruct = hypOptCell{bandLoop};
    mu = LightSrcOptmGP(testingX, curHypOptStruct);
    irrVisSim(:,:,bandLoop) = reshape(mu,rows);
end

XTestFr = reshape(ptsTestFrameHom(1,:),rows);
YTestFr = reshape(ptsTestFrameHom(2,:),rows);
ZTestFr = reshape(ptsTestFrameHom(3,:),rows);
lightSrc.PlotIrradianceCube(irrVisSim, XTestFr, YTestFr, ZTestFr,'GP light-src mean XZ');

for bandLoop = 1:numBands
    optPhi = optPhiBand(bandLoop, :);
    %create temporary light source using current optimised parameters
    tempLigSrcLS = LightSimulator(locLightSrc, rotLightSrc, optPhi(1), optPhi(2), optPhi(3));
    radIntMag = tempLigSrcLS.RadiantIntensityMesh(XTestFr, YTestFr, ZTestFr);
    irrVisSim(:,:,bandLoop) = radIntMag;
end

lightSrc.PlotIrradianceCube(irrVisSim, XTestFr, YTestFr, ZTestFr,'Least squares XZ');


%% CALLBACK functions
function band_callback(src, ~, scatMeas, surfLeasSqr, surfGP, ...
    radIntVP_Meas, radIntVP_LeasSqr, radIntVP_GP, sMeas, sLeasSqr, sGP)
i = round(src.Value);

scatMeas.ZData = radIntVP_Meas(:,i);
scatMeas.CData = radIntVP_Meas(:,i);

surfLeasSqr.CData = radIntVP_LeasSqr(:, :, i);
surfLeasSqr.ZData = radIntVP_LeasSqr(:, :, i);

surfGP.CData = radIntVP_GP(:, :, i);
surfGP.ZData = radIntVP_GP(:, :, i);

maxInt = max([ max(radIntVP_Meas(:,i)), max(radIntVP_GP(:,:,i), [], 'all'), max(radIntVP_LeasSqr(:,:,i), [], 'all')]);

caxis(sMeas, [0, maxInt]);
caxis(sLeasSqr, [0, maxInt]);
caxis(sGP, [0, maxInt]);

end

function band_callback2(src, ~, scatRadiance, scatDot, scatMeas, intenTarget, dotProductVP, radIntVP_Meas)
i = round(src.Value);

scatMeas.ZData = radIntVP_Meas(:,i);
scatMeas.CData = radIntVP_Meas(:,i);

scatDot.ZData = dotProductVP(:,i);
scatDot.CData = dotProductVP(:,i);

scatRadiance.ZData = intenTarget(:,i);
scatRadiance.CData = intenTarget(:,i);

end

function boardDispCheck_callback(src, ~, fillHandleCell)
warning('off');
global fillHandleCell; %#ok<REDEFGI>

value = get(src,'Value');

if value
    for i = 1:length(fillHandleCell)
        handlFill = fillHandleCell{i};
        handlFill.Visible = 'on';
    end
else
    for i = 1:length(fillHandleCell)
        handlFill = fillHandleCell{i};
        handlFill.Visible = 'off';
    end
end
warning('on');
end