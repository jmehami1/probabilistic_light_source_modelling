% Determines the estimated reflectance of hyperspectral images captured by 
% the line-scan hyperspectral and frame camera system. It will incorporate
% prior light source modelling and shape information (surface normals). The
% estimated reflectance methods are as follows:
% 
%   robles - solves Dichromatic model with equal reflectance patches 
%       Huynh, C. P., & Robles-Kelly, A. (2010). A solution of the dichromatic
%       model for multispectral photometric invariance. International Journal of
%       Computer Vision, 90(1), 1–27.
% 
%   krebs - separates the specular and diffuse components to find
%       reflectance
%       Krebs, A., Benezeth, Y., & Marzani, F. (2020). Intrinsic RGB and
%       multispectral images recovery by independent quadratic programming. PeerJ
%       Computer Science, 6, 1–15. https://doi.org/10.7717/PEERJ-CS.256
% 
%   krebsShape - krebs method with surface normal information (Our method)
% 
% Author: Jasprabhjit Mehami, 13446277

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

%MEX functions for ArUco pose estimation
addpath(genpath(['ext_lib', filesep, 'mex_aruco']));


%Contains implementation of dichromatic model estimation by Huynh and
%Robles-Kelly
addpath(genpath(['ext_lib', filesep, 'Scyllarus']));

%code for this project
addpath('src_code');

addpath(genpath('krebs-shape-reflectance-estimation'));
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
frameCamera = paramYaml.FrameCamera;
usingTestData = paramYaml.UseTestData;
testDataDir = paramYaml.TestDataDir;
hasOpenCV = paramYaml.HasOpenCV;

%STD in the radiance measurements
stdRadianceNoisePer = paramYaml.sigmaRadianceNoisePercent;
stdRadianceNoise = (stdRadianceNoisePer/100);

runRobles = paramYaml.RunRobles;
runKrebs = paramYaml.RunKrebs;
runKrebsShape = paramYaml.RunKrebsShape;


%maximum pixel intensity for uint14;
maxIntUint14 = 2^14 - 1;

%% Directories and files

%frame camera intrinsic parameters file
frameIntrFile = ['frame_camera_intrinsic', filesep, frameCamera, '.mat'];

% if frameCamera
%     frameIntrFile = ['frame_camera_intrinsic', filesep, 'blackfly.mat'];
% else
%     frameIntrFile = ['frame_camera_intrinsic', filesep, 'primesense.mat'];
% end

if ~exist(frameIntrFile, 'file')
    error("Frame camera intrinsic parameters not found");
end

if usingTestData
    sourDir = ['test_data', filesep, testDataDir];
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
    
    if ~exist(gpLightMappingFile, 'file')
        error('GP lighting model has not been built for test_data or missing MAT file');
    end
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

paramBoardYaml = yaml.ReadYaml(boardYamlFile);

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
xNumMarker = paramBoardYaml.NumCols;
yNumMarker = paramBoardYaml.NumRows;
% checkSize = pattern.CheckerSideLength;
arucoLen = paramBoardYaml.ArucoSideLength;
sepLen = paramBoardYaml.SeparationLength;
numMarkersExp = paramBoardYaml.NumberExpectedMarker;

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

% Light source triangulation parameters
%extract location and principal direction.
locLigSrc = cell2mat(ligSrcYaml.locLightSrc);
rotLigSrc = cell2mat(ligSrcYaml.rotLightSrc);
T_S_2_F = [rotLigSrc, locLigSrc; 0,0,0,1];
T_F_2_S = T_S_2_F\eye(4);

if displayOn
    %figure for light simulator
    figSim = figure('Name', 'Light Simulator with Frame Camera Origin');
    %plot frame camera at origin
    plotCamera('Location', zeros(1,3), 'Orientation', eye(3), 'Size', 0.05, 'AxesVisible', true); hold on;
    plotCamera('Location', T_LS_2_F(1:3,4), 'Orientation', T_LS_2_F(1:3, 1:3)', 'Size', 0.05, 'AxesVisible', true, 'Color', [0,0,1]); hold on;
    
    axis equal;
    grid on;
    xyzlabel();
    
    %light source simulator for visualisation
    lightSrc = LightSimulator(locLigSrc, rotLigSrc, figSim);
end


%Vertices of board in world coordinates
cornPatPts = [
    0, 0, 0;
    arucoLen*xNumMarker + sepLen*(xNumMarker-1), 0, 0;
    arucoLen*xNumMarker + sepLen*(xNumMarker-1), arucoLen*yNumMarker + sepLen*(yNumMarker-1), 0;
    0, arucoLen*yNumMarker + sepLen*(yNumMarker-1), 0;
    0, 0, 0;
    ];
cornPatPtsHom = [cornPatPts'; ones(1,size(cornPatPts,1))];

%size of sample region on board in world coordinates. This relates to the
%size of the cylinder used for curved surface
rSample = 10/100; %radius
lSample = 26/100; %length
%top left corner position of cylinder
xCorn = 4.15/100;
yCorn = 36.5/100;

if strcmp(testDataDir, 'deli_meats_curved')
    
    %start and end points of cylinder on cylinder axis
    xCylStart = [xCorn + rSample; yCorn - lSample; 0];
    xCylEnd = [xCorn + rSample; yCorn; 0];
    
    %create mesh of points on half cylinder surface
    nCylPts = 100;
    [X,Y,Z] = cylinder(rSample, nCylPts);
    
    %only keep one side of cylinder
    X = X(:,1:nCylPts/2 + 1);
    Y = Y(:,1:nCylPts/2 + 1);
    Z = Z(:,1:nCylPts/2 + 1);
    
    rowCol = size(X);
    X = X(:);
    Y = Y(:);
    Z = Z(:);
    Z = lSample.*Z; %scale by length
    
    %transform to pattern coordinates
    pts = [X'; Y'; Z'];
    pts = eul2rotm([0, 0, pi/2], 'ZYX')*pts;
    pts = pts + [xCorn + rSample; yCorn; 0];
    cylPtsHom = [pts; ones(1,size(pts,2))];
    
    if displayOn
        %plot the board, cylinder and axis in the coordinate frame of the frame
        %camera
        cylSurf = surf('EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', 0.3);
        
        ptsSCAT = scatter3(0,0,0, 10, 'k');
        normQUIV = quiver3(0,0,0,0,0,0, 'r-', 'LineWidth', 0.5);
        ptsSCAT.Visible = false;
        normQUIV.Visible = false;
        boardFill = fill3(1,1,1, 'c', 'FaceAlpha', 0.1);
        boardTR = trplot(eye(4), 'rviz', 'length', 0.05);
        view(-25,-40);
    end
    
    % [v, X, Y, Z, nx, ny, nz, band1, band2, ..., band146]
    pixRadianceData = zeros(100000, 7 + numBands);
    numPixPerLine = zeros(1,numImages);
    numData = 0;
    
    for imgLoop = 1:numImages
        curNormRadiImg = RadiImgCell{imgLoop};
        
        numPts = size(normPtImgLS,2);
        
        T_Board_2_F = extPosePattern(:,:,imgLoop);
        T_Board_2_LS = T_F_2_LS*T_Board_2_F;
        surfNormBoardLS = T_Board_2_LS(1:3,3);
        surfNormBoardFrame= T_Board_2_F(1:3,3);
        
        tBoardLS = T_Board_2_LS(1:3,4);
        
        
        %transform cylinder points to camera coordinate frame
        cylPtsFHom =  T_Board_2_F*cylPtsHom;
        X = reshape(cylPtsFHom(1,:), rowCol);
        Y = reshape(cylPtsFHom(2,:), rowCol);
        Z = reshape(cylPtsFHom(3,:), rowCol);
        
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
        
        
        %transform board plane to camera coordinate frame
        cornPatFrameHom = T_Board_2_F*cornPatPtsHom;
        
        if displayOn
            %update cylinder mesh points
            cylSurf.XData = X;
            cylSurf.YData = Y;
            cylSurf.ZData = Z;
            
            %update board points
            boardFill.XData =  cornPatFrameHom(1,:);
            boardFill.YData =  cornPatFrameHom(2,:);
            boardFill.ZData =  cornPatFrameHom(3,:);
            
            %update 3D axis to current extrinsic transformation
            trplot(boardTR, T_Board_2_F);
        end
        
        ptsCylFrame = zeros(numPts,3); %store line-scan pixel points that are on the cylinder
        ptsBoardSurfNorm = zeros(numPts,3); %store the corresponding normal at that point
        
        ptCount = 0;
        
        %trace pixel ray from the line-scan camera to the target frame to
        %determine its 3D location w.r.t LS
        for pixLoop = 1:numPts
            pnt = normPtImgLS(:,pixLoop);
            
            %intersect ray with cylinder
            [pntLS, valid] = CameraPixelRayIntersectCylinder(pnt', rSample, xClyStartLS', xCylEndLS');
            
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
                surfNormPt = x0 - xl;
                surfNormPt = surfNormPt./norm(surfNormPt);
                ptsBoardSurfNorm(ptCount, :) = surfNormPt';
                
                
                curBandRadi = curNormRadiImg(:, pixLoop);
                
                numData = numData + 1;
                pixRadianceData(numData, :) = [
                    pixLoop, pntF(1:3)', surfNormPt', curBandRadi'
                    ];
            end
        end
        
        numPixPerLine(imgLoop) = ptCount;
        
        %no points on cylinder
        if ptCount < 1
            if displayOn
                %turn off scatter and quiver visibility property
                ptsSCAT.Visible = false;
                normQUIV.Visible = false;
            end
            continue;
        end
        
        %clip arrays to size
        ptsCylFrame = ptsCylFrame(1:ptCount,:);
        ptsBoardSurfNorm = ptsBoardSurfNorm(1:ptCount,:);
        
        if displayOn
            %points on cylinder with surface normals
            s = [ptsCylFrame, ptsBoardSurfNorm];
            
            ptsSCAT.Visible = true;
            normQUIV.Visible = true;
            
            %update 3D scatter points on cylinder surface
            ptsSCAT.XData = s(:,1);
            ptsSCAT.YData = s(:,2);
            ptsSCAT.ZData = s(:,3);
            
            %downsampling for visualisation of quivers
            if ptCount > 50
                s = downsample(s, 50);
            end
            
            %update quivers to show surface normal on cylinder
            normQUIV.XData = s(:,1);
            normQUIV.YData = s(:,2);
            normQUIV.ZData = s(:,3);
            normQUIV.UData = s(:,4);
            normQUIV.VData = s(:,5);
            normQUIV.WData = s(:,6);
            
            drawnow();
        end
    end
    
elseif strcmp(testDataDir, 'deli_meats_flat')
    
    %vertices of rectangle region containing sample
    sampRecPts = [
        xCorn, yCorn, 0.001;
        xCorn + 2*rSample, yCorn, 0.001;
        xCorn + 2*rSample, yCorn - lSample, 0.001;
        xCorn, yCorn - lSample, 0.001;
        xCorn, yCorn, 0.001;
        ]';
    
    sampRecPtsHom = [sampRecPts; ones(1,5)];
    
    if displayOn
        sampFill = fill3(1,1,1, 'r', 'FaceAlpha', 0.1);
        ptsSCAT = scatter3(0,0,0, 10, 'k');
        normQUIV = quiver3(0,0,0,0,0,0, 'r-', 'LineWidth', 0.5);
        ptsSCAT.Visible = false;
        normQUIV.Visible = false;
        boardFill = fill3(1,1,1, 'c', 'FaceAlpha', 0.1);
        boardTR = trplot(eye(4), 'rviz', 'length', 0.05);
        view(-25,-40);
    end
    % [v, X, Y, Z, nx, ny, nz, band1, band2, ..., band146]
    pixRadianceData = zeros(100000, 7 + numBands);
    numPixPerLine = zeros(1,numImages);
    numData = 0;
    
    for imgLoop = 1:numImages
        curNormRadiImg = RadiImgCell{imgLoop};
        
        numPts = size(normPtImgLS,2);
        
        T_Board_2_F = extPosePattern(:,:,imgLoop);
        T_Board_2_LS = T_F_2_LS*T_Board_2_F;
        surfNormBoardLS = T_Board_2_LS(1:3,3);
        surfNormBoardFrame= T_Board_2_F(1:3,3);
        
        tBoardLS = T_Board_2_LS(1:3,4);
        
        sampRecPtsFrame = T_Board_2_F*sampRecPtsHom;
        
        %transform board plane to camera coordinate frame
        cornPatFrameHom = T_Board_2_F*cornPatPtsHom;
        
        if displayOn
            sampFill.XData = sampRecPtsFrame(1,:);
            sampFill.YData = sampRecPtsFrame(2,:);
            sampFill.ZData = sampRecPtsFrame(3,:);
            
            %update board points
            boardFill.XData =  cornPatFrameHom(1,:);
            boardFill.YData =  cornPatFrameHom(2,:);
            boardFill.ZData =  cornPatFrameHom(3,:);
            
            %update 3D axis to current extrinsic transformation
            trplot(boardTR, T_Board_2_F);
        end
        
        ptsCylFrame = zeros(numPts,3); %store line-scan pixel points that are on the cylinder
        ptsBoardSurfNorm = zeros(numPts,3); %store the corresponding normal at that point
        
        ptCount = 0;
        
        %trace pixel ray from the line-scan camera to the target frame to
        %determine its 3D location w.r.t LS
        for pixLoop = 1:numPts
            pnt = normPtImgLS(:,pixLoop);
            
            [pntLS, valid] = CameraPixelRayIntersectPlane(pnt', surfNormBoardLS', tBoardLS');
            
            %point is on flat infinite plane of the board
            if valid
                pntBoard = T_Board_2_LS\[pntLS';1];
                
                %check if point is within sample area
                if inpolygon(pntBoard(1), pntBoard(2), sampRecPts(1,:), sampRecPts(2,:))
                    
                    %transform point to frame camera
                    pntF  = T_LS_2_F*[pntLS';1];
                    surfNormPt = surfNormBoardFrame;
                    
                    %record point
                    ptCount = ptCount + 1;
                    ptsCylFrame(ptCount, :) = pntF(1:3);
                    ptsBoardSurfNorm(ptCount, :) = surfNormPt';
                    
                    curBandRadi = curNormRadiImg(:, pixLoop);
                    
                    numData = numData + 1;
                    pixRadianceData(numData, :) = [
                        pixLoop, pntF(1:3)', surfNormPt', curBandRadi'
                        ];
                    
                end
            end
        end
        
        numPixPerLine(imgLoop) = ptCount;
        
        %no points on cylinder
        if ptCount < 1
            if displayOn
                %turn off scatter and quiver visibility property
                ptsSCAT.Visible = false;
                normQUIV.Visible = false;
            end
            
            continue;
        end
        
        %clip arrays to size
        ptsCylFrame = ptsCylFrame(1:ptCount,:);
        ptsBoardSurfNorm = ptsBoardSurfNorm(1:ptCount,:);
        
        if displayOn
            %points on cylinder with surface normals
            s = [ptsCylFrame, ptsBoardSurfNorm];
            
            ptsSCAT.Visible = true;
            normQUIV.Visible = true;
            
            %update 3D scatter points on cylinder surface
            ptsSCAT.XData = s(:,1);
            ptsSCAT.YData = s(:,2);
            ptsSCAT.ZData = s(:,3);
            
            %downsampling for visualisation of quivers
            if ptCount > 50
                s = downsample(s, 50);
            end
            
            %update quivers to show surface normal on cylinder
            normQUIV.XData = s(:,1);
            normQUIV.YData = s(:,2);
            normQUIV.ZData = s(:,3);
            normQUIV.UData = s(:,4);
            normQUIV.VData = s(:,5);
            normQUIV.WData = s(:,6);
            
            drawnow();
        end
    end
else
    error('Unknown dataset');
end

pixRadianceData = pixRadianceData(1:numData, :);

%% Use GP light source model to predict radiant intensity

fprintf('Radiant intensity prediction using GP model...');

load(gpLightMappingFile);

%data iterator
pData = 1;

%radiant intensity hypercuber
radiantIntenHyperCube = zeros(rowsLS, numImages, numBands);

%normalised radiance hypercube that contains relevant pixels
normRadianceHyperCube = zeros(rowsLS, numImages, numBands);

numPts = sum(numPixPerLine, 'all');
testingX = zeros(numPts,2);


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
    
    %line-scan 3D points in frame camera coordinates
    ptFrame = pixRadianceData(pData:pData + numCurLine - 1, 2:4);
    %transform pts to light source coordinate frame
    ptsLightSrc = T_F_2_S*[ptFrame'; ones(1,numCurLine)];
    
    %convert points to spherical coordinates
    [ptsRadius,ptsTheta] = cart2sphZ(ptsLightSrc(1,:), ptsLightSrc(2,:), ptsLightSrc(3,:));
    
    testingX(pData:pData + numCurLine - 1, :) = [ptsRadius', ptsTheta']; %points used to predict radiant intensity
    
    pData = pData + numCurLine; %increment iterator
end

%store radiant intensity predictions for points
radiantIntLine = zeros(numPts, numBands);

%make prediction for each band
for bandLoop = 1:numBands
    hypOpt = hypOptCell{bandLoop};
    %     trainingXY = trainingXYCell{bandLoop};
    
    %prediction of radiant intensity at location from source
    [mu, varMu] = LightSrcOptmGP(testingX, hypOpt);
    radiantIntLine(:,bandLoop) = mu;
end

%data iterator
pData = 1;

for imgLoop = 1:numImages
    numCurLine = numPixPerLine(imgLoop); %number of relevant pixels on line
    
    %need atleast one point to determine radiant intensity
    if numCurLine < 1
        continue
    end
    
    pixLine = pixRadianceData(pData:pData + numCurLine - 1, 1);
    
    radiantIntenHyperCube(pixLine, imgLoop, :) = radiantIntLine(pData:pData + numCurLine - 1,:);
    pData = pData + numCurLine; %increment iterator
    
end

fprintf('Done\n');


%% Setup for estimated reflectance methods

cab('figSim');

figDir = [sourDir, filesep, 'Figures'];

if ~exist(figDir, 'dir')
    mkdir(figDir);
end

surfNormal = pixRadianceData(:, 5:7);
ptFrame = pixRadianceData(:, 2:4);

%direction light vector
ptLigVec = locLigSrc' - ptFrame;
dirLigVec = ptLigVec./vecnorm(ptLigVec, 2, 2);
dotLN = dot(surfNormal, dirLigVec, 2);

%normal map
surfImg = zeros(rowsLS, numImages, 3);

%light source direction map
ligImg = zeros(rowsLS, numImages, 3);

%image of dot product between light source direction and surface normal
lDotnImg = zeros(rowsLS, numImages);

%data iterator
pData = 1;

%filling images/maps for each relevant pixel-band
for imgLoop = 1:numImages
    numCurLine = numPixPerLine(imgLoop);
    
    %skip if no pixels
    if numCurLine < 1
        continue;
    end
    
    pixLine = pixRadianceData(pData:pData + numCurLine - 1, 1);
    pixLdotN = dotLN(pData:pData + numCurLine - 1);
    pixSurf = surfNormal(pData:pData + numCurLine - 1, :);
    pixLig = dirLigVec(pData:pData + numCurLine - 1, :);
    
    
    lDotnImg(pixLine, imgLoop) = pixLdotN;
    surfImg(pixLine, imgLoop,:) = pixSurf;
    ligImg(pixLine, imgLoop,:) = pixLig;
    
    pData = pData + numCurLine;
end


figure('Name', 'Surface normal map');
imshow(uint8(((surfImg+1)./2).*255));
title('Surface Normal Map');
drawnow();
savefig([figDir, filesep, 'surface_map.fig']);


figure('Name', 'Light source direction map on surface');
imshow(uint8(((ligImg+1)./2).*255));
title('Light Source Direction Map');
drawnow();
savefig([figDir, filesep, 'light_source_map.fig']);

lDotnImg(lDotnImg < 0) = 0; %ignore any negative dotproduct
figure('Name', 'Dot product between light source ray and surface normal');
imagesc(lDotnImg);
title('Dot Product between Light Ray and Surface Normals'); colormap(turbo);colorbar();
drawnow();
savefig([figDir, filesep, 'dot_lightsource_surfnormal.fig']);

%image mask using the dot product image
binMaskG = lDotnImg > 0;

%estimated radiant intensity per band from radiance hypercube using grey
%world method
radiantIntenEstimatedBand = estimate_lightsource(normRadianceHyperCube, 'grey_world')';

hyperDir = [sourDir, filesep, 'Hypercubes'];

%% Robles least squares patch method

cab('figSim');

if runRobles
    
    method = 'robles';
    
    fprintf('**********************************\n');
    fprintf(['Estimated reflectance using ', method, ' method:\n']);
    
    normRadianceHyperCubeMasked = bsxfun(@times, normRadianceHyperCube, binMaskG);
    radiantIntenHyperCubeMasked = bsxfun(@times, radiantIntenHyperCube, binMaskG);
    
    
    fprintf('\tRunning hypercube with estimated radiant intensity...');
    
    %estimated reflectance using estimated radiant intensity per band
    [k, g, S] = recover_dichromatic_parameters_LS(normRadianceHyperCubeMasked, radiantIntenEstimatedBand, 0, 5, 0.01);
    hyperCube = S;
    
    save([hyperDir, filesep, 'estimated_reflectance_', method, '.mat'], 'hyperCube');
    figure();
    subplot(1,2,1)
    imagesc(k); title('Specular Coefficient'); colormap(gray);colorbar;
    subplot(1,2,2)
    imagesc(g); title('Shading Factor'); colormap(gray);colorbar;
    sgtitle(method);
    savefig([figDir, filesep, 'gk_', method, '.fig']);
    
    fprintf('DONE\n');
    
    fprintf('\tRunning hypercube with GP radiant intensity...');
    
    %estimated reflectance using radiant intensity hypercube from GP
    [k, g, S] = recover_dichromatic_parameters_LS(normRadianceHyperCubeMasked, radiantIntenHyperCubeMasked, 0, 5, 0.01);
    hyperCube = S;
    
    save([hyperDir, filesep, 'estimated_reflectance_', method, '_gp.mat'], 'hyperCube');
    figure();
    subplot(1,2,1)
    imagesc(k); title('Specular Coefficient'); colormap(gray);colorbar;
    subplot(1,2,2)
    imagesc(g); title('Shading Factor'); colormap(gray);colorbar;
    sgtitle([method, ' with GP Light Source']);
    savefig([figDir, filesep, 'gk_', method, '_gp.fig']);
    
    fprintf('DONE\n');
end

%% Estimate reflectance using Krebs method

cab('figSim');

if runKrebs
    
    method = 'krebs';
    
    fprintf('**********************************\n');
    fprintf(['Estimated reflectance using ', method, ' method:\n']);
    
    %radiance divided by radiant intensity from grey world method(effects of illumination removed)
    tempRadi(1,:,:) = radiantIntenEstimatedBand';
    R_hyperCube = bsxfun(@rdivide, normRadianceHyperCube, tempRadi);
    R_hyperCube(isnan(R_hyperCube)) = 0; %ignore any nan
    R_hyperCube(isinf(R_hyperCube)) = 0; %ignore any inf
    
    
    fprintf('\tRunning hypercube with estimated radiant intensity...');
    
    [S, g, k] = method_krebs_logversion(R_hyperCube, binMaskG);
    hyperCube = S;
    
    save([hyperDir, filesep, 'estimated_reflectance_', method, '.mat'], 'hyperCube');
    figure();
    subplot(1,2,1)
    imagesc(k); title('Specular Coefficient'); colormap(gray);colorbar;
    subplot(1,2,2)
    imagesc(g); title('Shading Factor'); colormap(gray);colorbar;
    sgtitle(method);
    savefig([figDir, filesep, 'gk_', method, '.fig']);
    
    fprintf('DONE\n');
    
    %radiance divided by radiant intensity from GP(effects of illumination removed)
    R_hyperCube = normRadianceHyperCube./radiantIntenHyperCube;
    R_hyperCube(isnan(R_hyperCube)) = 0; %ignore any nan
    R_hyperCube(isinf(R_hyperCube)) = 0; %ignore any inf
    
    fprintf('\tRunning hypercube with GP radiant intensity...');
    
    [S, g, k] = method_krebs_logversion(R_hyperCube, binMaskG);
    hyperCube = S;
    
    save([hyperDir, filesep, 'estimated_reflectance_', method, '_gp.mat'], 'hyperCube');
    figure();
    subplot(1,2,1)
    imagesc(k); title('Specular Coefficient'); colormap(gray);colorbar;
    subplot(1,2,2)
    imagesc(g); title('Shading Factor'); colormap(gray);colorbar;
    sgtitle([method, ' with GP Light Source']);
    savefig([figDir, filesep, 'gk_', method, '_gp.fig']);
    
    fprintf('DONE\n');
    
end


%% Estimate reflectance using Krebs method with shape

cab('figSim');

if runKrebsShape
    
    method = 'krebsShape';
    
    fprintf('**********************************\n');
    fprintf(['Estimated reflectance using ', method, ' method:\n']);
    fprintf('\tRunning hypercube with GP radiant intensity...');
    
    
    [S, k, g] = ReflectanceEstimationKrebs(normRadianceHyperCube,radiantIntenHyperCube, binMaskG, lDotnImg, 3, 0.01);
    hyperCube = S;
    
    save([hyperDir, filesep, 'estimated_reflectance_', method, '.mat'], 'hyperCube');
    figure();
    subplot(1,2,1)
    imagesc(k); title('Specular Coefficient'); colormap(gray);colorbar;
    subplot(1,2,2)
    imagesc(g); title('Shading Factor'); colormap(gray);colorbar;
    sgtitle(method);
    savefig([figDir, filesep, 'gk_', method, '.fig']);
    
    fprintf('DONE\n');
    
end
