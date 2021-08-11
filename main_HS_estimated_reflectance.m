% Determines the estimated reflectance of images captured by the line-scan
% hyperspectral and frame camera system with a modelled light source.

%Author: Jasprabhjit Mehami, 13446277

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

%Better pose estimation of a planar board
addpath(genpath(['ext_lib', filesep, 'IPPE']));

%Contains implementation of dichromatic model estimation by Huynh and
%Robles-Kelly
addpath(genpath(['ext_lib', filesep, 'Scyllarus']));

%code for this project
addpath('src_code');

%parameter file
paramFile = ['parameter_files', filesep, 'HS_estimated_reflectance.yaml'];
if ~exist(paramFile, 'file')
    error("YAML parameter file not found");
end

%% Read YAML file containing the pattern specifications and parameters for code
% All dimensions are in metres

paramYaml = yaml.ReadYaml(paramFile);
displayOn = paramYaml.DisplayOn;
usingBlackfly = paramYaml.UsingBlackFly;
usingTestData = paramYaml.UseTestData;

%STD in the radiance measurements
stdRadianceNoise = paramYaml.sigmaRadianceNoise;

%maximum pixel intensity for uint14;
maxIntUint14 = 2^14 - 1;

%% Directories and files

%frame camera intrinsic parameters file
if usingBlackfly
    frameIntrFile = ['frame_camera_intrinsic', filesep, 'blackfly.mat'];
else
    frameIntrFile = ['frame_camera_intrinsic', filesep, 'primesense.mat'];
end

if ~exist(frameIntrFile, 'file')
    error("Frame camera intrinsic parameters not found");
end

if usingTestData
    sourDir = ['test_data', filesep, 'sample_platform_curved'];
else
    %Get source directory where images are located and results will be saved
    sourDir = uigetdir(['~', filesep], 'Provide source directory where light modelling images are located?');
end

%frame image directory
frameDir = [sourDir, filesep, 'Frame']; %Directory containing images
if ~exist(frameDir, 'dir')
    error('Source directory does not contain directory called "Frame" which should contain RGB images');
end

fullPathFrame = fullfile([frameDir, filesep, '*.png']);

%Need to get number of images in directory
numImagesFrame = numel(dir(fullPathFrame));

if numImagesFrame < 1
    error('no images in provided image directory')
end

%hyperspectral line-scan image directory
hsDir = [sourDir, filesep, 'Line-scan_corrected']; %Directory containing images
if ~exist(hsDir, 'dir')
    error('Source directory does not contain directory called "Line-scan" which should contain hyperspectral images');
end

fullPathHS = fullfile([hsDir, filesep, '*.png']);

%Need to get number of images in directory
numImagesHS = numel(dir(fullPathHS));

if numImagesHS < 1
    error('no images in provided image directory')
end

if numImagesHS ~= numImagesFrame
    error("number of images in folders are not the same");
end

numImages = numImagesHS;


% %result directory
% resultDir = [sourDir, filesep, 'Result'];
% if ~exist(resultDir, 'dir')
%     mkdir(resultDir);
% end


if usingTestData
    ligTrigPath = ['test_data', filesep, 'light_trig', filesep, 'Result', filesep];
    ligTrigFile = 'pt_light.yaml';
    
    gpLightMappingFile = ['test_data', filesep, 'light_map', filesep, 'Result', filesep, 'gp_lightsrc_optm.mat'];
    
    load(gpLightMappingFile);
else
    %get light triangulation results parameter file from user
    [ligTrigFile, ligTrigPath] = uigetfile([ligTrigStartDir, '*.yaml'], 'Provide light triangulation results YAML file', 'MultiSelect', 'off');
end

ligSrcYaml = yaml.ReadYaml([ligTrigPath, ligTrigFile]);


%line-scan frame camera calibration parameters file
cameraSysCaliFile = ['camera_system_calibration', filesep, 'camera_system_optimised_parameters.mat'];
if ~exist(cameraSysCaliFile, 'file')
    error('camera system optimised parameters file not found');
end

%% Load images for both cameras

fprintf('Loading images...');

%Preallocate space for cell arrays of images
imagesLS = cell(1,numImages);
imagesFrame = cell(1,numImages);

% Load all images
for i = 1:numImages
    imagesFrame{i} = im2gray(imread([frameDir, filesep, 'img', num2str(i),'.png']));
    imagesLS{i} = imread([hsDir, filesep, 'hs', num2str(i),'.png']);
end

imageDark = imread([sourDir, filesep, 'dark_ref_corrected.png']);

fprintf('Done\n');

%% Load the board YAML file

boardYamlFile = [sourDir, filesep, 'platform_board.yaml'];

if ~exist(boardYamlFile, 'file')
    error('YAML parameter file is not present in directory');
end

paramYaml = yaml.ReadYaml(boardYamlFile);

%% Determine extrinsic of from each image

fprintf('Estimating extrinisic pose of each image...');

load(frameIntrFile); %Load the intrinsic parameters of the camera

%extract focal point components in pixels
fx = cameraParams.FocalLength(1);
fy = cameraParams.FocalLength(2);

%optical centre of camera in pixelstrue
u0 = cameraParams.PrincipalPoint(1);
v0 = cameraParams.PrincipalPoint(2);

%array of intrinsic parameters that introduce noise (assume noise in
%intrinsic distortion is zero)
thetaFrameintr = [fx, fy, u0, v0];

%Size of the images from the RGB camera
frameImgSize = size(imagesFrame{1});
frameImgSize = frameImgSize(1:2);

%distortion parameters
distRad = cameraParams.RadialDistortion;
distTan = cameraParams.TangentialDistortion;
distCoefCV = [distRad(1:2), distTan, distRad(3)]; %array of distortion coefficients in opencv format

%extract the error in the intrinsic parameters of the frame camera
% std_f = estimationErrors.IntrinsicsErrors.FocalLengthError;
% std_u0v0 = estimationErrors.IntrinsicsErrors.PrincipalPointError;
% stdRGB_Intr = [std_f,std_u0v0];


%ChArUco pattern size
xNumMarker = paramYaml.NumCols;
yNumMarker = paramYaml.NumRows;
% checkSize = pattern.CheckerSideLength;
arucoLen = paramYaml.ArucoSideLength;
sepLen = paramYaml.SeparationLength;
numMarkersExp = paramYaml.NumberExpectedMarker;

%intrinsic object for the RGB camera
frameIntrinsic = cameraIntrinsics(thetaFrameintr(1:2),thetaFrameintr(3:4), frameImgSize);
K_frame = frameIntrinsic.IntrinsicMatrix';

%store all the poses of each found pattern
extPosePattern = zeros(4,4,numImages);

if displayOn
    fig = figure('Name','ChArUco pattern pose');
end

% 2D marker corner positions in world coordinates (metres)
markerCornerCell = ArUcoBoardMarkerCornersCell(1, xNumMarker, yNumMarker, arucoLen, sepLen);

%used to filter images where the pose can't be found
goodImages = zeros(1,numImages);
numGoodImg = 0;


for imgLoop = 1:numImages
    
    [rotMat, trans, found, imgDisp] = ArucoPosEst(imagesFrame{imgLoop}, markerCornerCell, cameraParams);
    
    if ~found
        continue;
    end
    
    %image is good
    numGoodImg = numGoodImg + 1;
    goodImages(numGoodImg) = imgLoop;
    
    %store found extrinsic parameter
    extPosePattern(:,:,numGoodImg) = [rotMat,trans'; 0, 0, 0, 1];
    
    % display the frame camera image with the projected axis on the pattern
    if displayOn
        clf(fig);
        imshow(imgDisp); hold on;
        
        Plot3DAxis2Image(extPosePattern(:,:,numGoodImg), arucoLen, K_frame, frameImgSize, []);
        
        drawnow();
    end
end

%remove all data from the frame images where we could not find proper
%extrinsic parameters
goodImages = goodImages(1:numGoodImg);
imagesFrame = imagesFrame(goodImages);
imagesLS = imagesLS(goodImages);
extPosePattern = extPosePattern(:,:,1:numGoodImg);
numImages = numGoodImg;

fprintf('Done\n');

%% Load calibration parameters for the line-scan frame camera system

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
K_LS = intrLS.IntrinsicMatrix';
K_LS(1,3) = 0;

%field of view (radians)
fovLS = 2*atan(rowsLS/(2*fy));

%field of view of pixel (radians)
fovLSPix = fovLS/rowsLS;

%% Get normalised line-scan image coordinates and their normalised radiance

fprintf('Normalising line-scan pixel coordinates and radiance...');

% relevant pixel locations in line-scan image
bandStart = 40;
bandEnd = 203;
hypBands = bandStart:bandEnd;
numBands = length(hypBands);
numBandPix = 256;

normPtsLSCell = cell(1,numImages);
normRadianceImgCell = cell(1,numImages);

for imgLoop = 1:numImages
    curImg = imagesLS{imgLoop};
    
    vImgLineDist = 1:rowsLS; %distorted v component pixel positions on image line
    imgLineDist = [zeros(size(vImgLineDist));vImgLineDist]'; %distorted pixels on image line (u = 0)
    imgLine = undistortPoints(imgLineDist, cameraParamLS); %undistort image line
    
    %determine normalised pixel coordinates
    imgLineHom = [imgLine'; ones(1,length(imgLineDist(:,1)))];
    normPtImgLS = K_LS\imgLineHom;
    normPtsLSCell(imgLoop) = {normPtImgLS'};
    
    %normalise between dark reference and max intensity of uint14 for
    %relevant bands    
    normRadianceImgCell(imgLoop) = {NormaliseRadianceBands(curImg, maxIntUint14.*ones(size(curImg), 'uint16'), imageDark, bandStart, bandEnd)};  
end

fprintf('Done\n');


%% Find intersection with planar pattern and cylinder

% Simulator figure
close all;

%figure for light simulator
figSim = figure('Name', 'Light Simulator with Frame Camera Origin');
%plot frame camera at origin
plotCamera('Location', zeros(1,3), 'Orientation', eye(3), 'Size', 0.05, 'AxesVisible', true); hold on;
plotCamera('Location', T_LS_2_F(1:3,4), 'Orientation', T_LS_2_F(1:3, 1:3)', 'Size', 0.05, 'AxesVisible', true, 'Color', [0,0,1]); hold on;

axis equal;
grid on;
xyzlabel();

% Light source triangulation parameters

%extract location and principal direction.
locLigSrc = cell2mat(ligSrcYaml.locLightSrc);
rotLigSrc = cell2mat(ligSrcYaml.rotLightSrc);
T_S_2_F = [rotLigSrc, locLigSrc; 0,0,0,1];
T_F_2_S = T_S_2_F\eye(4);

%light source simulator for visualisation
lightSrc = LightSimulator(locLigSrc, rotLigSrc, figSim);

%Vertices of board in world coordinates
cornPatPts = [
    0, 0, 0;
    arucoLen*xNumMarker + sepLen*(xNumMarker-1), 0, 0;
    arucoLen*xNumMarker + sepLen*(xNumMarker-1), arucoLen*yNumMarker + sepLen*(yNumMarker-1), 0;
    0, arucoLen*yNumMarker + sepLen*(yNumMarker-1), 0;
    0, 0, 0;
    ];
cornPatPtsHom = [cornPatPts'; ones(1,size(cornPatPts,1))];

%parameters of cylinder in world coordinates
rCyl = 10/100; %radius
lCyl = 26/100; %length
%top left corner position of cylinder
xCorn = 4.15/100;
yCorn = 36.5/100;
%start and end points of cylinder on cylinder axis
xCylStart = [xCorn + rCyl; yCorn - lCyl; 0];
xCylEnd = [xCorn + rCyl; yCorn; 0];

%create mesh of points on half cylinder surface
nCylPts = 100;
[X,Y,Z] = cylinder(rCyl, nCylPts);

%only keep one side of cylinder
X = X(:,1:nCylPts/2 + 1);
Y = Y(:,1:nCylPts/2 + 1);
Z = Z(:,1:nCylPts/2 + 1);

rowCol = size(X);
X = X(:);
Y = Y(:);
Z = Z(:);
Z = lCyl.*Z; %scale by length

%transform to pattern coordinates
pts = [X'; Y'; Z'];
pts = eul2rotm([0, 0, pi/2], 'ZYX')*pts;
pts = pts + [xCorn + rCyl; yCorn; 0];
cylPtsHom = [pts; ones(1,size(pts,2))];

%plot the board, cylinder and axis in the coordinate frame of the frame
%camera
cylSurf = surf('EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.3);
boardFill = fill3(1,1,1, 'c', 'FaceAlpha', 0.1);
boardTR = trplot(eye(4), 'rviz', 'length', 0.05);
ptsCylSCAT = scatter3(0,0,0, 10, 'k');
normCylQUIV = quiver3(0,0,0,0,0,0, 'r-', 'LineWidth', 0.5);
ptsCylSCAT.Visible = false;
normCylQUIV.Visible = false;
view(-25,-40);

% [v, X, Y, Z, nx, ny, nz, band1, band2, ..., band146]
pixRadianceData = zeros(100000, 7 + numBands);
numPixPerLine = zeros(1,numImages);
numData = 0;

for imgLoop = 1:numImages
    curNormRadiImg = normRadianceImgCell{imgLoop};
    
    figure(figSim);
    
    normHSPts = normPtsLSCell{imgLoop}';
    numPts = size(normHSPts,2);
    
    T_Board_2_F = extPosePattern(:,:,imgLoop);
    T_Board_2_LS = T_F_2_LS*T_Board_2_F;
    
    %transform cylinder points to camera coordinate frame
    cylPtsFHom =  T_Board_2_F*cylPtsHom;
    X = reshape(cylPtsFHom(1,:), rowCol);
    Y = reshape(cylPtsFHom(2,:), rowCol);
    Z = reshape(cylPtsFHom(3,:), rowCol);
    
    %update cylinder mesh points
    cylSurf.XData = X;
    cylSurf.YData = Y;
    cylSurf.ZData = Z;
    
    %transform board plane to camera coordinate frame
    cornPatFrameHom = T_Board_2_F*cornPatPtsHom;
    %update board points
    boardFill.XData =  cornPatFrameHom(1,:);
    boardFill.YData =  cornPatFrameHom(2,:);
    boardFill.ZData =  cornPatFrameHom(3,:);
    
    %update 3D axis to current extrinsic transformation
    trplot(boardTR, T_Board_2_F);
    
    transBoard = T_Board_2_LS(1:3, 4);
    surfNormBoard = T_Board_2_LS(1:3,3);
    
    %cylinder start and end points in the coordinate frame of the line-scan
    %camera
    xClyStartLS = T_Board_2_LS*[xCylStart;1];
    xClyStartLS = xClyStartLS(1:3);
    xCylEndLS = T_Board_2_LS*[xCylEnd;1];
    xCylEndLS = xCylEndLS(1:3);
    
    %cylinder start and end points in the coordinate frame of the frame
    %camera
    xClyStartF = T_Board_2_F*[xCylStart;1];
    xClyStartF = xClyStartF(1:3);
    xCylEndF = T_Board_2_F*[xCylEnd;1];
    xCylEndF = xCylEndF(1:3);
    
    
    ptsCylFrame = zeros(numPts,3); %store line-scan pixel points that are on the cylinder
    ptsBoardSurfNorm = zeros(numPts,3); %store the corresponding normal at that point
    
    ptCount = 0;
    
    %trace pixel ray from the line-scan camera to the target frame to
    %determine its 3D location w.r.t LS
    for pixelLoop = 1:numPts
        pnt = normHSPts(:,pixelLoop);
        
        %intersect ray with cylinder
        [pntLS, valid] = CameraPixelRayIntersectCylinder(pnt', rCyl, xClyStartLS', xCylEndLS');
        
        %point is on cylinder
        if valid
            ptCount = ptCount + 1;
            
            %transform point to frame camera
            pntF  = T_LS_2_F*[pntLS';1];
            ptsCylFrame(ptCount, :) = pntF(1:3);
            
            %find the t-parameter of the point on the cylinder axis that will yield the minimum
            %distance to the intersection.
            x0 = pntF(1:3);
            x1 = xClyStartF;
            x2 = xCylEndF;
            t = linePosition3d(x0, [x1', (x2 - x1)']);
            
            %point on cylinder axis with minimum distance to intersection.
            %The line created by this point and the intersection on
            %cylinder surface is normal to the axis
            xl = x1 + (x2-x1).*t;
            
            %calculate surface normal direction vector
            surfNormCyl = x0 - xl;
            surfNormCyl = surfNormCyl./norm(surfNormCyl);
            ptsBoardSurfNorm(ptCount, :) = surfNormCyl';
            
            
            curBandRadi = curNormRadiImg(hypBands, pixelLoop);
            
            numData = numData + 1;
            pixRadianceData(numData, :) = [
                pixelLoop, pntF(1:3)', surfNormCyl', curBandRadi'
                ];
        end
    end
    
    numPixPerLine(imgLoop) = ptCount;
    
    %no points on cylinder
    if ptCount < 1
        %turn off scatter and quiver visibility property
        ptsCylSCAT.Visible = false;
        normCylQUIV.Visible = false;
        continue;
    end
    
    %clip arrays to size
    ptsCylFrame = ptsCylFrame(1:ptCount,:);
    ptsBoardSurfNorm = ptsBoardSurfNorm(1:ptCount,:);
    
    %points on cylinder with surface normals
    s = [ptsCylFrame, ptsBoardSurfNorm];
    
    ptsCylSCAT.Visible = true;
    normCylQUIV.Visible = true;
    
    %update 3D scatter points on cylinder surface
    ptsCylSCAT.XData = s(:,1);
    ptsCylSCAT.YData = s(:,2);
    ptsCylSCAT.ZData = s(:,3);
    
    %downsampling for visualisation of quivers
    if ptCount > 10
        s = downsample(s, 10);
    end
    
    %update quivers to show surface normal on cylinder
    normCylQUIV.XData = s(:,1);
    normCylQUIV.YData = s(:,2);
    normCylQUIV.ZData = s(:,3);
    normCylQUIV.UData = s(:,4);
    normCylQUIV.VData = s(:,5);
    normCylQUIV.WData = s(:,6);
    
    drawnow();
end

pixRadianceData = pixRadianceData(1:numData, :);

%% Use GP light source model to predict radiant intensity

fprintf('Radiant intensity prediction using GP model...');

%data iterator
pData = 1;

%Store radiant intensity for each line-image
radiantIntenCell = cell(1,numImages);

%radiant intensity hypercuber
radiantIntenHyperCube = zeros(rowsLS, numImages, numBands);

%normalised radiance hypercube that contains relevant pixels
normRadianceHyperCube = zeros(rowsLS, numImages, numBands);

for imgLoop = 1:numImages
    numCurLine = numPixPerLine(imgLoop); %number of relevant pixels on line
    
    %need atleast one point to determine radiant intensity
    if numCurLine < 1        
        continue
    end
    
    pixLine = pixRadianceData(pData:pData + numCurLine - 1, 1);
    curNormRadiance = pixRadianceData(pData:pData + numCurLine - 1, 8:end);
    
    %store normalised radiance
    normRadianceHyperCube(pixLine, imgLoop, :) = curNormRadiance;
    
    %store radiant intensity predictions for points
    radiantIntLine = zeros(numBands, numCurLine);
    
    %line-scan 3D points in frame camera coordinates
    ptFrame = pixRadianceData(pData:pData + numCurLine - 1, 2:4);
    %transform pts to light source coordinate frame
    ptsLightSrc = T_F_2_S*[ptFrame'; ones(1,numCurLine)];
    
    %convert points to spherical coordinates
    [ptsRadius,ptsTheta] = cart2sphZ(ptsLightSrc(1,:), ptsLightSrc(2,:), ptsLightSrc(3,:));
    testingX = [ptsRadius', ptsTheta']; %points used to predict radiant intensity
        
    %make prediction for each band
    for bandLoop = 1:numBands
        hypOpt = hypOptCell{bandLoop};
        %prediction of radiant intensity at location from source
        [mu, varMu] = LightSrcOptmGP(3, trainingXY, testingX, downSamplingGP, 100000, stdRadianceNoise, false, hypOpt);
        radiantIntLine(bandLoop,:) = mu';
    end
    
    radiantIntenCell(imgLoop) = {radiantIntLine};
    radiantIntenHyperCube(pixLine, imgLoop, :) = radiantIntLine';
    
    pData = pData + numCurLine; %increment iterator
end

%average radiant intensity per band
radiIntBandAver = squeeze(sum(radiantIntenHyperCube, [1,2])) ./ squeeze(sum(radiantIntenHyperCube > 0, [1,2]));

fprintf('Done\n');


%% Plotting hypercubes

cab('hyperCubeFig', 'figSim');

hyperCubeFig = figure('Name', 'Hypercube visualisation');

hsVis = Hypercube_Visualiser(hyperCubeFig);
hsVis.AddHypercube(radiantIntenHyperCube, "Radiant intensity from GP");
hsVis.AddHypercube(normRadianceHyperCube, "Normalised Radiance");


%% Estimate reflectance using minimisation with no added information

[specularImgOR, shadingFactorImgOR, reflectanceOR] = recover_dichromatic_parameters_LS(normRadianceHyperCube, radiIntBandAver);

hsVis.AddHypercube(reflectanceOR, "Reflectance Original Method")
% clc;
% lineT = 212;
% pixT = 39;
% bandT = 100;
% 
% l = radiantIntenHyperCube(pixT, lineT, bandT);
% G(pixT, lineT)*l*S(pixT, lineT, bandT) + K(pixT, lineT)*l
% normRadianceHyperCube(pixT, lineT, bandT)

%% Calculate estimated reflectance for each line measurement

pData = 1;


%radiant intensity magnitude produced by light source for each line image. There is a
%value at each pixel-band


normRadCube = zeros(rowsLS, numImages);
lightFieldVarCube = zeros(rowsLS, numImages);

reflecHyperCube = zeros(rowsLS, numImages, numBands);

locLS = T_LS_2_F(1:3,4);

% optOptions = optimoptions(@lsqlin, 'Algorithm','trust-region-reflective', 'Display', 'none');
optOptions = optimoptions(@lsqlin, 'Algorithm','interior-point', 'Display', 'none');

% figRes = figure('Name', 'Estimated reflectance result');

maxPixDist = 2;

for imgLoop = 1:numImages
    
    numCurLine = numPixPerLine(imgLoop);

    %Need atleast two point measurements from hyperspectral line
    if numCurLine < 2
        continue;
    end
    
    pixLine = pixRadianceData(pData:pData + numCurLine - 1, 1);
    
    %need atleast one pair of consective pixels that are within maxPixDist
    %apart
    pixDiffAbs = abs(diff(pixLine));
    if ~any(pixDiffAbs <= maxPixDist, 'all')
       continue; 
    end
    
    tic();
    
    radiantIntLine = radiantIntenCell{imgLoop};
    ptFrame = pixRadianceData(pData:pData + numCurLine - 1, 2:4);
    PtFrameHom = [ptFrame'; ones(1,numCurLine)];
    
    %transform pts to light source coordinate frame
    ptsLightSrc = T_F_2_S*PtFrameHom;
    
    %distance from light source to points
    ptsDistLigSrc = cart2sphZ(ptsLightSrc(1,:), ptsLightSrc(2,:), ptsLightSrc(3,:));
    
    %direction light vector
    ptLigVec = locLigSrc' - ptFrame;
    dirLigVec = ptLigVec./vecnorm(ptLigVec, 2, 2); 
    
    
    ptViewVec = locLS' - ptFrame;
    dirViewVec = ptViewVec./vecnorm(ptViewVec, 2, 2); %direction camera viewing vector
    ptsDistLS = vecnorm(ptViewVec,2,2); %distance from line-scan camera to points
    

    
    surfNormal = pixRadianceData(pData:pData + numCurLine - 1, 5:7);
    normRadiance = pixRadianceData(pData:pData + numCurLine - 1, 8:end);
    
    
    %length of sector for each pixel projected onto the distance for each
    %point. Assuming a square pixel, calculate the foreshortened area and
    %then find the projected area onto the surface for each pixel.
    sidePixLS = ptsDistLS.*fovLSPix;
    pixProjArea = sidePixLS.^2 ./ dot(dirViewVec, surfNormal, 2);

    estRefl = SolveEstimatedReflectance(normRadiance', radiantIntLine, surfNormal, ...
        dirLigVec, dirViewVec, ptsDistLigSrc', pixProjArea, pixLine, 5, optOptions, maxPixDist);
    
    reflecHyperCube(pixLine, imgLoop, :) = estRefl';
    
    
    pData = pData + numCurLine;
    
    tElapsed = toc();
    
    fprintf("Elapsed time: %f \t number of pixels: %i\n",  tElapsed, numCurLine);
end
 
% normRadCube = mat2gray(normRadCube);

img = reflecHyperCube(:,:,100);

imshow(img);


lightImg = mat2gray(radiantIntLine);
lightVarImg = mat2gray(lightFieldVarCube);

figure()
subplot(1,3, 1);
imshow(lightImg);
title('normalised light intensity');

subplot(1,3, 2);
imshow(lightVarImg);
colormap(gca, jet(256)); % Ignore pink map and use jet instead.
colorbar(gca);
title('normalised variance');


subplot(1,3, 3);
imshow(normRadCube);
% colormap(gca, jet(256)); % Ignore pink map and use jet instead.
% colorbar(gca);
title('normalised radiance');


