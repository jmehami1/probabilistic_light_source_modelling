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

addpath(genpath(['ext_lib', filesep, 'krebs-shape-reflectance-estimation']));
addpath(genpath(['ext_lib', filesep, 'krebs-reflectance-estimation']));

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
hasOpenCV = paramYaml.HasOpenCV;

%STD in the radiance measurements
stdRadianceNoisePer = paramYaml.sigmaRadianceNoisePercent;
stdRadianceNoise = (stdRadianceNoisePer/100);

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

%white stripe reflectance CSV file
whiteStripeReflFile = ['parameter_files', filesep, 'white_stripe_reflectance.csv'];
if ~exist(paramFile, 'file')
    error("white stripe reflectance CSV not found");
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
imageWhite = imread([sourDir, filesep, 'white_ref_corrected.png']);

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

extMatFile = [sourDir, filesep, 'imgExtrinsicPose.mat'];

if hasOpenCV
    
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
    
    goodImages = goodImages(1:numGoodImg);
    extPosePattern = extPosePattern(:,:,1:numGoodImg);
 
    save(extMatFile, 'extPosePattern', 'goodImages', 'numGoodImg')
else
    load(extMatFile);
end

%remove all data from the frame images where we could not find proper
%extrinsic parameters
imagesFrame = imagesFrame(goodImages);
imagesLS = imagesLS(goodImages);
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

RadiImgCell = cell(1,numImages);


vImgLineDist = 1:rowsLS; %distorted v component pixel positions on image line
imgLineDist = [zeros(size(vImgLineDist));vImgLineDist]'; %distorted pixels on image line (u = 0)
imgLine = undistortPoints(imgLineDist, cameraParamLS); %undistort image line

%determine normalised pixel coordinates
imgLineHom = [imgLine'; ones(1,length(imgLineDist(:,1)))];
normPtImgLS = K_LS\imgLineHom;


for imgLoop = 1:numImages
    curImg = imagesLS{imgLoop};

    radiImg = curImg - imageDark;%subtract dark reference
    radiImg = double(radiImg(hypBands, :))./maxIntUint14; %scale down
    
    RadiImgCell(imgLoop) = {radiImg};
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
    curNormRadiImg = RadiImgCell{imgLoop};
    
    figure(figSim);
    
    numPts = size(normPtImgLS,2);
    
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
    
%     transBoard = T_Board_2_LS(1:3, 4);
%     surfNormBoard = T_Board_2_LS(1:3,3);
    
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
        pnt = normPtImgLS(:,pixelLoop);
        
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
            
            
            curBandRadi = curNormRadiImg(:, pixelLoop);
            
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
    if ptCount > 50
        s = downsample(s, 50);
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
        trainingXY = trainingXYCell{bandLoop};
        
        %prediction of radiant intensity at location from source
        [mu, varMu] = LightSrcOptmGP(3, trainingXY, testingX, stdRadianceNoise, hypOpt);
        radiantIntLine(bandLoop,:) = mu';
    end
    
    radiantIntenHyperCube(pixLine, imgLoop, :) = radiantIntLine';
    
    pData = pData + numCurLine; %increment iterator
end

%average radiant intensity per band
radiIntBandAver = squeeze(sum(radiantIntenHyperCube, [1,2])) ./ squeeze(sum(radiantIntenHyperCube > 0, [1,2]));

fprintf('Done\n');

%% Standardise hyperspectral measurements using white reference

% %read the reflectance values for each band of the hyperspectral camera
% %stored in file
% data = readmatrix(whiteStripeReflFile);
% targetReflectances = data(:,2);
% 
% 
% normWhiteRef = NormaliseRadianceBands(imageWhite, maxIntUint14.*ones(size(imageWhite), 'uint16'), imageDark, bandStart, bandEnd);
% normWhiteRef = normWhiteRef(bandStart:bandEnd,:)';
% whiteRefHyperCube = zeros(numPts,1,numBands);
% whiteRefHyperCube(:,1,:) = normWhiteRef;
% 
% T_Board_2_F = extPosePattern(:,:,1);
% T_Board_2_LS = T_F_2_LS*T_Board_2_F;
% 
% transBoard = T_Board_2_LS(1:3, 4);
% surfNormBoard = T_Board_2_LS(1:3,3);
% 
% normHSPts = normPtImgLS{1}';
% 
% ptsWhiteFrame = zeros(numPts, 3);
% 
% %trace pixel ray from the line-scan camera to the target frame to
% %determine its 3D location w.r.t LS
% for pixelLoop = 1:numPts
%     pnt = normHSPts(:,pixelLoop);
%     
%     pntLS = CameraPixelRayIntersectPlane(pnt', surfNormBoard', transBoard');
%     
%     pntF  = T_LS_2_F*[pntLS';1];
%     ptsWhiteFrame(pixelLoop, :) = pntF(1:3);
% end
% 
% ptsLightSrc = T_F_2_S*[ptsWhiteFrame'; ones(1,numPts)];
% 
% %convert points to spherical coordinates
% [ptsRadius,ptsTheta] = cart2sphZ(ptsLightSrc(1,:), ptsLightSrc(2,:), ptsLightSrc(3,:));
% testingX = [ptsRadius', ptsTheta']; %points used to predict radiant intensity
% 
% %store radiant intensity predictions for points
% radiantIntLine = zeros(numBands, numPts);
% 
% %make prediction for each band
% for bandLoop = 1:numBands
%     hypOpt = hypOptCell{bandLoop};
%     %prediction of radiant intensity at location from source
%     [mu, varMu] = LightSrcOptmGP(3, trainingXY, testingX, downSamplingGP, 1000, stdRadianceNoise, false, hypOpt);
%     radiantIntLine(bandLoop,:) = mu';
% end
% 
% radiantIntenWhiteHyperCube = zeros(numPts,1,numBands);
% 
% radiantIntenWhiteHyperCube(:,1,:) = radiantIntLine';
% 
% %direction light vector
% ptLigVec = locLigSrc' - ptsWhiteFrame;
% dirLigVec = ptLigVec./vecnorm(ptLigVec, 2, 2);
% 
% surfNormal = T_Board_2_F(1:3,3)'.*ones(size(dirLigVec));
% 
% dirL_dot_n = dot(surfNormal, dirLigVec, 2);
% 
% addpath(genpath(['ext_lib', filesep, 'krebs-reflectance-estimation']));
% 
% 
% alphaBalFactor = EstimateBalancingFactor(whiteRefHyperCube, radiantIntenWhiteHyperCube, dirL_dot_n, targetReflectances);


%% Estimate reflectance using Krebs method

cab('figSim');

surfNormal = pixRadianceData(:, 5:7);
ptFrame = pixRadianceData(:, 2:4);

%direction light vector
ptLigVec = locLigSrc' - ptFrame;
dirLigVec = ptLigVec./vecnorm(ptLigVec, 2, 2);
dirL_dot_n = dot(surfNormal, dirLigVec, 2);

%image of dot product between light source direction and surface normal
lDotnImg = zeros(rowsLS, numImages);

%iterator
pData = 1;

%filling lDotnImg for each relevant pixel-band
for imgLoop = 1:numImages
    numCurLine = numPixPerLine(imgLoop);
    
    %skip if no pixels
    if numCurLine < 1
        continue;
    end
    
    pixLine = pixRadianceData(pData:pData + numCurLine - 1, 1);
    pixLdotN = dirL_dot_n(pData:pData + numCurLine - 1);
    
    lDotnImg(pixLine, imgLoop) = pixLdotN;
    pData = pData + numCurLine;
end

lDotnImg(lDotnImg < 0) = 0; %ignore any negative dotproduct

figure;imagesc(lDotnImg);

%image mask using the dot product image
binMaskG = lDotnImg > 0;


%% Estimate reflectance using Krebs method with shape

fprintf('**********************************\n');
fprintf('Estimate reflectance using Krebs with shape information method...');


[S_krebsShape, k_krebsShape, g_krebsShape] = ReflectanceEstimationKrebs(normRadianceHyperCube,radiantIntenHyperCube, binMaskG, lDotnImg, 3, pi*pi/3600);

fprintf('Done\n');

%% Estimate reflectance krebs no shape information

fprintf('**********************************\n');
fprintf('Estimate reflectance using Krebs method...');

%radiance with radiant intensity removed
R_hyperCube = normRadianceHyperCube./radiantIntenHyperCube;
R_hyperCube(isnan(R_hyperCube)) = 0; %ignore any nan
R_hyperCube(isinf(R_hyperCube)) = 0; %ignore any inf

[S_krebs, g_Krebs, k_krebs] = method_krebs_logversion(R_hyperCube, binMaskG);

fprintf('Done\n');

%% Robles least squares patch method

fprintf('**********************************\n');
fprintf('Estimate reflectance using Robles method...');

normRadianceHyperCubeMasked = bsxfun(@times, normRadianceHyperCube, binMaskG);
radiantIntenHyperCubeMasked = bsxfun(@times, radiantIntenHyperCube, binMaskG);

[k_robles, g_robles, S_robles] = recover_dichromatic_parameters_LS(normRadianceHyperCubeMasked, radiantIntenHyperCubeMasked);

fprintf('Done\n');

%% Save hypercube and images

hyperDir = [sourDir, filesep, 'Hypercubes'];

hyperCube = S_krebs;
save([hyperDir, filesep, 'estimated_reflectance_krebs.mat'], 'hyperCube');

hyperCube = S_krebsShape;
save([hyperDir, filesep, 'estimated_reflectance_krebs_with_shape.mat'], 'hyperCube');

hyperCube = S_robles;
save([hyperDir, filesep, 'estimated_reflectance_robles.mat'], 'hyperCube');

%%

cab('figSim');


figure;imagesc(k_krebs); title('k Krebs'); colorbar;
figure;imagesc(g_Krebs); title('g Krebs'); colorbar;

figure;imagesc(lDotnImg); title('n dot l');


figure;imagesc(k_krebsShape); title('k Krebs with Shape'); colorbar;
figure;imagesc(g_krebsShape); title('g Krebs with Shape'); colorbar;


figure;imagesc(k_robles); title('k Robles'); colorbar;
figure;imagesc(g_robles); title('g Robles'); colorbar;
% hyperspectralViewer(S_krebs.*15);
% hyperspectralViewer(S_new./alphaBalFactor);


% %% Estimate reflectance using minimisation with no added information
% 
% % cab('figSim');
% 
% % radiIntBandAver = estimate_lightsource(normRadianceHyperCube)';
% % 
% % 
% % [kAverRadi, GAverRadi, phoAverRadi] = recover_dichromatic_parameters_LS(normRadianceHyperCube, radiIntBandAver);
% % 
% % hyperspectralViewer(phoAverRadi);
% 
% [k_Robles, GFullRadi, phoFullRadi] = recover_dichromatic_parameters_LS(normRadianceHyperCube, radiantIntenHyperCube);
% 
% GFullRadi(isnan(GFullRadi)) = 0;
% 
% hyperspectralViewer(phoFullRadi);
% 
% 
% hyperCubeFig = figure('Name', 'Hypercube visualisation');
% 
% hsVis = Hypercube_Visualiser(hyperCubeFig);
% hsVis.AddHypercube(phoAverRadi, "Reflectance Constant band radiant");
% hsVis.AddHypercube(phoFullRadi, "Reflectance Full Radiant");
% 
% %%
% cab('figSim', 'hyperCubeFig');
% 
% figKG = figure('Name', 'Shading factor and specular images');
% 
% subplot(2,2,1);
% imshow(mat2gray(GAverRadi));
% title('Constant band radiant')
% 
% subplot(2,2,2);
% imshow(mat2gray(kAverRadi));
% 
% subplot(2,2,3);
% imshow(mat2gray(GFullRadi));
% title('Full radiant')
% 
% subplot(2,2,4);
% imshow(mat2gray(k_Robles));
% 
% %%
% 
% 
% %% Calculate Shading factor of pixel 3D locations
% 
% cab('figSim');
% 
% surfNormal = pixRadianceData(:, 5:7);
% 
% ptFrame = pixRadianceData(:, 2:4);
% 
% %direction light vector
% ptLigVec = locLigSrc' - ptFrame;
% dirLigVec = ptLigVec./vecnorm(ptLigVec, 2, 2);
% 
% %transform pts to light source coordinate frame
% PtFrameHom = [ptFrame'; ones(1,size(ptFrame,1))];
% ptsLightSrc = T_F_2_S*PtFrameHom;
% 
% %distance from light source to points
% ptsDistLigSrc = cart2sphZ(ptsLightSrc(1,:), ptsLightSrc(2,:), ptsLightSrc(3,:))';
%     
% %distance from line-scan camera to points
% locLS = T_LS_2_F(1:3,4);
% ptViewVec = locLS' - ptFrame;
% ptsDistLS = vecnorm(ptViewVec,2,2);
% 
% %direction camera viewing vector
% dirViewVec = ptViewVec./ptsDistLS; 
%     
% %length of sector for each pixel projected onto the distance for each
% %point. Assuming a square pixel, calculate the foreshortened area and
% %then find the projected area onto the surface for each pixel.
% sidePixLS = ptsDistLS.*fovLSPix;
% pixProjArea = sidePixLS.^2 ./ dot(dirViewVec, surfNormal, 2);
% pixProjArea(pixProjArea < 0) = 0;
% 
% %shading factor of all pixels
% gVec = ShadingFactorPixels(surfNormal, dirLigVec, ptsDistLigSrc, pixProjArea);
% 
% % %minimum shading factor (these pixels won't provide enough constraints)
% % Gmin = 0.01;
% % 
% % %mask to remove any pixels with low G
% % Gmask = gVec > Gmin;
% 
% %shading factor image
% imageG = zeros(rowsLS, numImages);
% 
% %iterator
% pData = 1;
% 
% for imgLoop = 1:numImages
%     numCurLine = numPixPerLine(imgLoop);
%     
%     %skip if no pixels
%     if numCurLine < 1
%         continue;
%     end
%     
%     pixLine = pixRadianceData(pData:pData + numCurLine - 1, 1);   
%     pixG = gVec(pData:pData + numCurLine - 1);
%     imageG(pixLine, imgLoop) = pixG;
%     
%     pData = pData + numCurLine;
% end
% 
% figure;imagesc(imageG);
% 
% binMaskG = imageG > 0;
% rpG = regionprops(binMaskG, 'Centroid', 'BoundingBox');
% 
% figure;imshow(binMaskG); hold on;
% 
% BB = floor(rpG(1).BoundingBox);
% rectangle('Position', [BB(1),BB(2),BB(3),BB(4)],'EdgeColor','r','LineWidth',2) ;
% 
% rpG_crop = imcrop(binMaskG, BB);
% figure;imshow(rpG_crop);
% 
% patchSize = [5, 5];
% 
% 
% rpG_tileCell = mat2tiles(rpG_crop, patchSize);
% [tRows, tCols] = size(rpG_tileCell);
% 
% 
% %apply mask to radiant intensity hypercube
% radiantIntenCubeMasked = radiantIntenHyperCube.*binMaskG;
% radiantIntenCubeMasked_crop = imcrop3(radiantIntenCubeMasked, [BB(1:2), 1, BB(3:4), numBands - 1]);
% radiantIntenCubeMasked_crop_tileCells = mat2tiles(radiantIntenCubeMasked_crop, [patchSize, Inf]);
% 
% normRadianceCubeMasked = normRadianceHyperCube.*binMaskG;
% normRadianceCubeMasked_crop = imcrop3(normRadianceHyperCube, [BB(1:2), 1, BB(3:4), numBands - 1]);
% normRadianceCubeMasked_crop_tileCells = mat2tiles(normRadianceCubeMasked_crop, [patchSize, Inf]);
% 
% imgG_crop = imcrop(imageG, BB);
% imgG_crop_tileCells = mat2tiles(imgG_crop, patchSize);
% 
% rhoImg_crop_tileCells = cell(tRows, tCols);
% 
% for i = 1:tRows
%    for j = 1:tCols
%        curMask = rpG_tileCell{i,j};
%        numValidPix = sum(curMask, 'all');
%        
%        if numValidPix < 2
%            rhoPatch = zeros(size(curMask));
%            rhoImg_crop_tileCells(i,j) = {rhoPatch};
%            continue
%        end
%        
%        SolveEstimatedReflectancePatch(curMask)
%    end
% end
% 
% 
%     estRefl = SolveEstimatedReflectance(normRadiance', radiantIntLine, surfNormal, ...
%         dirLigVec, dirViewVec, ptsDistLigSrc', pixProjArea, pixLine, 5, optOptions, maxPixDist);
%     
%     reflecHyperCube(pixLine, imgLoop, :) = estRefl';
%     
%     
%     pData = pData + numCurLine;
% 
% 
% 
% normRadCube = zeros(rowsLS, numImages);
% lightFieldVarCube = zeros(rowsLS, numImages);
% 
% reflecHyperCube = zeros(rowsLS, numImages, numBands);
% 
% 
% % optOptions = optimoptions(@lsqlin, 'Algorithm','trust-region-reflective', 'Display', 'none');
% optOptions = optimoptions(@lsqlin, 'Algorithm','interior-point', 'Display', 'none');
% 
% % figRes = figure('Name', 'Estimated reflectance result');
% 
% maxPixDist = 2;
% 
% for imgLoop = 1:numImages
%     
%     numCurLine = numPixPerLine(imgLoop);
% 
%     %Need atleast two point measurements from hyperspectral line
%     if numCurLine < 2
%         continue;
%     end
%     
%     pixLine = pixRadianceData(pData:pData + numCurLine - 1, 1);
%     
%     %need atleast one pair of consective pixels that are within maxPixDist
%     %apart
%     pixDiffAbs = abs(diff(pixLine));
%     if ~any(pixDiffAbs <= maxPixDist, 'all')
%        continue; 
%     end
%     
% %     tic();
%     
%     radiantIntLine = radiantIntenCell{imgLoop};
%     ptFrame = pixRadianceData(pData:pData + numCurLine - 1, 2:4);
%     PtFrameHom = [ptFrame'; ones(1,numCurLine)];
%     
%     %transform pts to light source coordinate frame
%     ptsLightSrc = T_F_2_S*PtFrameHom;
%     
%     %distance from light source to points
%     ptsDistLigSrc = cart2sphZ(ptsLightSrc(1,:), ptsLightSrc(2,:), ptsLightSrc(3,:));
%     
%     %direction light vector
%     ptLigVec = locLigSrc' - ptFrame;
%     dirLigVec = ptLigVec./vecnorm(ptLigVec, 2, 2); 
%     
%     
% 
%     
% 
%     
%     surfNormal = pixRadianceData(pData:pData + numCurLine - 1, 5:7);
%     normRadiance = pixRadianceData(pData:pData + numCurLine - 1, 8:end);
%     
%     
%     %length of sector for each pixel projected onto the distance for each
%     %point. Assuming a square pixel, calculate the foreshortened area and
%     %then find the projected area onto the surface for each pixel.
%     sidePixLS = ptsDistLS.*fovLSPix;
%     pixProjArea = sidePixLS.^2 ./ dot(dirViewVec, surfNormal, 2);
% 
%     estRefl = SolveEstimatedReflectance(normRadiance', radiantIntLine, surfNormal, ...
%         dirLigVec, dirViewVec, ptsDistLigSrc', pixProjArea, pixLine, 5, optOptions, maxPixDist);
%     
%     reflecHyperCube(pixLine, imgLoop, :) = estRefl';
%     
%     
%     pData = pData + numCurLine;
%     
% %     tElapsed = toc();
%     
% %     fprintf("Elapsed time: %f \t number of pixels: %i\n",  tElapsed, numCurLine);
% end
% 
% 
% 
% 
% 
% %% Plotting hypercubes
% 
% % cab('hyperCubeFig', 'figSim');
% % 
% % hyperCubeFig = figure('Name', 'Hypercube visualisation');
% % 
% % hsVis = Hypercube_Visualiser(hyperCubeFig);
% % hsVis.AddHypercube(radiantIntenHyperCube, "Radiant intensity from GP");
% % hsVis.AddHypercube(normRadianceHyperCube, "Normalised Radiance");
% % 
% % 
% % %% Estimate reflectance using minimisation with no added information
% % % 
% % % [specularImgOR, shadingFactorImgOR, reflectanceOR] = recover_dichromatic_parameters_LS(normRadianceHyperCube, radiIntBandAver);
% % % 
% % % hyperspectralViewer(reflectanceOR);
% % 
% % [specularImgOR, shadingFactorImgOR, reflectanceNEW] = recover_dichromatic_parameters_LS(normRadianceHyperCube, radiantIntenHyperCube, 5, 0.00);
% % 
% % hyperspectralViewer(reflectanceNEW);
% 
% 
% % hsVis.AddHypercube(reflectanceOR, "Reflectance Original Method")
% % clc;
% % lineT = 212;
% % pixT = 39;
% % bandT = 100;
% % 
% % l = radiantIntenHyperCube(pixT, lineT, bandT);
% % G(pixT, lineT)*l*S(pixT, lineT, bandT) + K(pixT, lineT)*l
% % normRadianceHyperCube(pixT, lineT, bandT)
% 
% %% 
% 
% 
% %% Calculate estimated reflectance for each line measurement
% 
% pData = 1;
% 
% 
% %radiant intensity magnitude produced by light source for each line image. There is a
% %value at each pixel-band
% 
% 
% normRadCube = zeros(rowsLS, numImages);
% lightFieldVarCube = zeros(rowsLS, numImages);
% 
% reflecHyperCube = zeros(rowsLS, numImages, numBands);
% 
% locLS = T_LS_2_F(1:3,4);
% 
% % optOptions = optimoptions(@lsqlin, 'Algorithm','trust-region-reflective', 'Display', 'none');
% optOptions = optimoptions(@lsqlin, 'Algorithm','interior-point', 'Display', 'none');
% 
% % figRes = figure('Name', 'Estimated reflectance result');
% 
% maxPixDist = 2;
% 
% for imgLoop = 1:numImages
%     
%     numCurLine = numPixPerLine(imgLoop);
% 
%     %Need atleast two point measurements from hyperspectral line
%     if numCurLine < 2
%         continue;
%     end
%     
%     pixLine = pixRadianceData(pData:pData + numCurLine - 1, 1);
%     
%     %need atleast one pair of consective pixels that are within maxPixDist
%     %apart
%     pixDiffAbs = abs(diff(pixLine));
%     if ~any(pixDiffAbs <= maxPixDist, 'all')
%        continue; 
%     end
%     
% %     tic();
%     
%     radiantIntLine = radiantIntenCell{imgLoop};
%     ptFrame = pixRadianceData(pData:pData + numCurLine - 1, 2:4);
%     PtFrameHom = [ptFrame'; ones(1,numCurLine)];
%     
%     %transform pts to light source coordinate frame
%     ptsLightSrc = T_F_2_S*PtFrameHom;
%     
%     %distance from light source to points
%     ptsDistLigSrc = cart2sphZ(ptsLightSrc(1,:), ptsLightSrc(2,:), ptsLightSrc(3,:));
%     
%     %direction light vector
%     ptLigVec = locLigSrc' - ptFrame;
%     dirLigVec = ptLigVec./vecnorm(ptLigVec, 2, 2); 
%     
%     
%     ptViewVec = locLS' - ptFrame;
%     dirViewVec = ptViewVec./vecnorm(ptViewVec, 2, 2); %direction camera viewing vector
%     ptsDistLS = vecnorm(ptViewVec,2,2); %distance from line-scan camera to points
%     
% 
%     
%     surfNormal = pixRadianceData(pData:pData + numCurLine - 1, 5:7);
%     normRadiance = pixRadianceData(pData:pData + numCurLine - 1, 8:end);
%     
%     
%     %length of sector for each pixel projected onto the distance for each
%     %point. Assuming a square pixel, calculate the foreshortened area and
%     %then find the projected area onto the surface for each pixel.
%     sidePixLS = ptsDistLS.*fovLSPix;
%     pixProjArea = sidePixLS.^2 ./ dot(dirViewVec, surfNormal, 2);
% 
%     estRefl = SolveEstimatedReflectance(normRadiance', radiantIntLine, surfNormal, ...
%         dirLigVec, dirViewVec, ptsDistLigSrc', pixProjArea, pixLine, 5, optOptions, maxPixDist);
%     
%     reflecHyperCube(pixLine, imgLoop, :) = estRefl';
%     
%     
%     pData = pData + numCurLine;
%     
% %     tElapsed = toc();
%     
% %     fprintf("Elapsed time: %f \t number of pixels: %i\n",  tElapsed, numCurLine);
% end
%  
% % normRadCube = mat2gray(normRadCube);
% 
% img = reflecHyperCube(:,:,100);
% 
% imshow(img);
% 
% 
% lightImg = mat2gray(radiantIntLine);
% lightVarImg = mat2gray(lightFieldVarCube);
% 
% figure()
% subplot(1,3, 1);
% imshow(lightImg);
% title('normalised light intensity');
% 
% subplot(1,3, 2);
% imshow(lightVarImg);
% colormap(gca, jet(256)); % Ignore pink map and use jet instead.
% colorbar(gca);
% title('normalised variance');
% 
% 
% subplot(1,3, 3);
% imshow(normRadCube);
% % colormap(gca, jet(256)); % Ignore pink map and use jet instead.
% % colorbar(gca);
% title('normalised radiance');
% 
% 
