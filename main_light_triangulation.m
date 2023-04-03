%Triangulation of real point source's principal direction and location
%relative to a frame camera origin. This uses a reflective metal sphere on
%a ChArUco board. By fitting an ellipse to the brightest spot on the sphere
%we can find a constraint to determine the unknown direction and location.

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

%Better export of figures to images
addpath(fullfile('ext_lib', 'export_fig'))

%code for this project
addpath('src_code');

%check if mex_ChArUco_Pose has been built
if ~exist(fullfile("ext_lib", "mex_ChArUco_Pose", "bin", "CharucoPosEst.mexa64"), 'file')
    error("Please build mex_ChArUco_Pose submodule")
else
    addpath(genpath(fullfile("ext_lib", "mex_ChArUco_Pose")));
end

%source code for this project
addpath('src_code');

%parameter file
paramFile = fullfile('parameter_files', 'light_triangulation.yaml');
if ~exist(paramFile, 'file')
    error("YAML parameter file not found");
end


%% Read YAML file containing the pattern specifications and parameters for code
% All dimensions are in metres

paramYaml = yaml.ReadYaml(paramFile);
displayOn = paramYaml.DisplayOn;
frameCamera = lower(paramYaml.FrameCamera);
usingTestData = paramYaml.UseTestData;

%% Directories and files

disp('Checking directories...')

%frame camera intrinsic parameters file
frameIntrFile = fullfile('frame_camera_intrinsic', [frameCamera, '.mat']);

if ~exist(frameIntrFile, 'file')
    error('Frame camera intrinsic parameters not found for %s', frameCamera);
end

%Get source directory where images are located and results will be saved
if usingTestData
    sourDir = fullfile('test_data', 'light_trig');
else
    sourDir = uigetdir('Provide source directory where light triangulation images are located?');
end

%frame image directory
frameDir = fullfile(sourDir, 'Frame'); %Directory containing images
if ~exist(frameDir, 'dir')
    error('Source directory does not contain directory called "Frame" which should contain RGB images');
end

fullPathFrame = fullfile(frameDir, '*.png');

%Need to get number of images in directory
numImages = numel(dir(fullPathFrame));

if numImages < 1
    error('no images in provided image directory')
end

fprintf("Found directory with %i images\n", numImages);

%result directory
resultDir = fullfile(sourDir, 'Result');
if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end

%% Load intrinsic parameters and board parameters

disp('Loading intrinsic and board parameters...');

% load(frameIntrFile); %Load the intrinsic parameters of the camera

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
frameImgSize = cameraParams.ImageSize;
% frameImgSize = frameImgSize(1:2);

%distortion parameters
distRad = cameraParams.RadialDistortion;
distTan = cameraParams.TangentialDistortion;
distCoefCV = zeros(1,5); %array of distortion coefficients in opencv format

%extract the error in the intrinsic parameters of the frame camera
if exist('estimationErrors', 'var')
    %extract the error in the intrinsic parameters of the frame camera
    std_f = estimationErrors.IntrinsicsErrors.FocalLengthError;
    std_u0v0 = estimationErrors.IntrinsicsErrors.PrincipalPointError;
    stdRGB_Intr = [std_f,std_u0v0];
else
    stdRGB_Intr = zeros(1,4);
end

%ChArUco pattern size
xNumCheck = paramYaml.NumCols;
yNumCheck = paramYaml.NumRows;
checkSize = paramYaml.CheckerSideLength;
arucoSize = paramYaml.ArucoSideLength;

%intrinsic object for the RGB camera
frameIntrinsic = cameraIntrinsics(thetaFrameintr(1:2),thetaFrameintr(3:4), frameImgSize);
K_frame = frameIntrinsic.IntrinsicMatrix;

%% Load images

disp('Loading images...');

%Preallocate space for cell arrays of images
imagesFrame = cell(1,numImages);

% Load all images
for i = 1:numImages
    imagesFrame{i} = undistortImage(imread(fullfile(frameDir, ['img', num2str(i),'.png'])), cameraParams);
end


%% Get pose of the plane using the ChArUco pattern
disp('Getting pose of ChArUco pattern...');

%store all the poses of each found pattern
extPosePattern = zeros(4,4,numImages);

%again, used to filter images where the pose can't be found
goodImages = zeros(1,numImages);
numGoodImg = 0;

if displayOn
    fig = figure('Name','ChArUco pattern pose');
end

for i = 1:numImages
    %MEX file for getting the pose using the ChArUco pattern
    [rotMat, trans, found, img] = CharucoPosEst(imagesFrame{i}, K_frame, distCoefCV, ...
        xNumCheck, yNumCheck, checkSize, arucoSize);

    %if pose not found, ignore this image
    if ~found
        continue;
    end

    %image is good
    numGoodImg = numGoodImg + 1;
    goodImages(numGoodImg) = i;

    %store found extrinsic parameter
    extPosePattern(:,:,numGoodImg) = [rotMat,trans'; 0, 0, 0, 1];

    %display the frame camera image with the projected axis on the pattern
    if displayOn
        clf(fig);
        imshow(img); hold on;
    end
end

%remove all data from the frame images where we could not find proper
%extrinsic parameters
goodImages = goodImages(1:numGoodImg);
imagesFrame = imagesFrame(goodImages);
extPosePattern = extPosePattern(:,:,1:numGoodImg);
numImages = numGoodImg;

%% Find brighest pixel location on reflective hemisphere by fitting an ellipse

close all;

disp('Finding brightest pixel on reflective hemisphere images...');


%relevant dimensions of reflective hemisphere on board
xCent = paramYaml.Xhemi;
yCent = paramYaml.Yhemi;
hemiRad = paramYaml.Dhemi/2;

%points on the edge of the reflective hemisphere w.r.t to the board
noPts = 100;
theta = linspace(0, 2*pi, noPts);
xHemiPts = xCent + hemiRad.*cos(theta);
yhemiPts = yCent + hemiRad.*sin(theta);

patPts = [
    xHemiPts;
    yhemiPts;
    zeros(1,noPts);
    ];
patPtsHom = [patPts; ones(1,noPts)];

if displayOn
    figG = figure('Name','gray');
end

reflCentImg = zeros(numImages, 2);

SE = strel('disk',1);

%again, used to filter images where the brightest spot can't be found
goodImages = zeros(1,numImages);
numGoodImg = 0;

for imgLoop = 1:numImages
    if displayOn
        clf(figG);
    end

    %Get the pixel coordinates of the edge of the reflective sphere on the
    %board
    ext = extPosePattern(:,:,imgLoop);
    patPtsHomFrame = ext*patPtsHom;
    patPtsFrame = patPtsHomFrame(1:3,:)';
    patImgPts = projectPoints(patPtsFrame, K_frame', eye(4), [], frameImgSize);

    %clip pixels to remain within the size of the frame image
    for row = 1:noPts
        pixU = round(patImgPts(row,1));

        if pixU < 1
            pixU = 1;
        elseif pixU > frameImgSize(2)
            pixU = frameImgSize(2);
        end

        patImgPts(row,1) = pixU;

        pixV = round(patImgPts(row,2));

        if pixV < 1
            pixV = 1;
        elseif pixV > frameImgSize(1)
            pixV = frameImgSize(1);
        end

        patImgPts(row,2) = pixV;
    end

    patImgPts = patImgPts';

    %convert color to grayscale
    imgGray = im2gray(imagesFrame{imgLoop});


    imgGraySphere = zeros(size(imgGray), 'uint8');

    %extract on the pixels that contain the reflective sphere
    for v = min(patImgPts(2, :), [], 'all'):max(patImgPts(2, :), [], 'all')
        for u = min(patImgPts(1, :), [], 'all'):max(patImgPts(1, :), [], 'all')
            if (inpolygon(u, v, patImgPts(1,:), patImgPts(2,:)))
                imgGraySphere(v, u) = imgGray(v, u);
            end
        end
    end

    %convert grayscale to normalised grayscale
    imgGraySphereNorm = mat2gray(imgGraySphere);

    %thresold all pixels on the sphere
    imgBrSpot = imgGraySphereNorm > 0.95;

    %erode
    imgBrSpot = imerode(imgBrSpot, SE);

    %find regions
    s = regionprops(imgBrSpot, {'Centroid','Orientation','MajorAxisLength','MinorAxisLength', 'Area'});

    %no region found
    if length(s) < 1
        if displayOn
            set(0, 'CurrentFigure', figG);
            imshow(imgGray);hold on;
            title('FAILED bright spot detection');
            hold off;
            drawnow();
        end
        
        continue;
    end
    %find the region with the largest area
    ind = 1;
    maxArea = s(1).Area;

    for i = 2:length(s)
        if s(i).Area > maxArea
            ind = i;
            maxArea = s(i).Area;
        end
    end

    numGoodImg = numGoodImg + 1;
    goodImages(numGoodImg) = imgLoop;

    %get the centre of the region
    regCent = s(ind).Centroid;
    reflCentImg(numGoodImg, :) = regCent;

    % Calculate the ellipse line
    theta = linspace(0,2*pi, noPts);
    col = (s(ind).MajorAxisLength/2)*cos(theta);
    row = (s(ind).MinorAxisLength/2)*sin(theta);
    M = makehgtform('translate',[s(ind).Centroid, 0],'zrotate',deg2rad(-1*s(ind).Orientation));
    D = M*[col;row;zeros(1,numel(row));ones(1,numel(row))];

    if displayOn
        set(0, 'CurrentFigure', figG);
        imshow(imgGray);hold on;
        plot(D(1,:),D(2,:),'r','LineWidth',1);
        plot(regCent(1),regCent(2),'b+','LineWidth',1);
        hold off
        drawnow();
    end
end

%remove all data from the frame images where we could not find proper
%bright spot
goodImages = goodImages(1:numGoodImg);
imagesFrame = imagesFrame(goodImages);
extPosePattern = extPosePattern(:,:,goodImages);
reflCentImg = reflCentImg(1:numGoodImg, :);
numImages = numGoodImg;

%% Find location of brighest spot on reflective sphere in world coordinates

disp('Finding location of brightest spot in world coordinates...');

%centre of the hemisphere relative to the pattern
hemiCent = [xCent; yCent; 0];

ptReflHemiFrame = zeros(numImages, 3);
normReflHemi = zeros(numImages, 3);

for imgLoop = 1:numImages
    centPtImg =  reflCentImg(imgLoop, :)';
    centPtImgHom = [centPtImg; 1];

    %normalised image coordinates
    normPt = K_frame'\centPtImgHom;

    %centre of hemisphere relative to the frame camera
    ext = extPosePattern(:,:,imgLoop);
    locPatFrame = tform2trvec(ext)';
    rotPatFrame = tform2rotm(ext);
    locHemiCentFrame = locPatFrame + rotPatFrame*hemiCent;

    [pt, rc] = intersect_Line_Sphere([[0,0,0], normPt'], [locHemiCentFrame', hemiRad]);

    if ~rc
        error("line-sphere intersection was not on sphere")
    end

    ptReflHemiFrame(imgLoop, :) = pt;

    %surface normal at the intersection point on the reflective hemisphere
    normHemi = pt' - locHemiCentFrame;
    normHemi = normHemi./norm(normHemi);
    normReflHemi(imgLoop, :) = normHemi;
end

%% Find the light source direction and location

disp('Finding point source direction and location...');

[locLightSrc, dirLightSrc] = SolveLightSrcDirLoc(ptReflHemiFrame, normReflHemi, numImages);

%principal direction is along the z-axis, find how much the z-axis has
%rotated from its original direction. This rotation would be about some
%perpindular axes
%Go from direction vector to rotation matrix
zAxis = [0;0;1];
perpenDir = cross(zAxis, dirLightSrc);
DirTheta = atan2(norm(cross(dirLightSrc,zAxis)),dot(dirLightSrc,zAxis));
rotLightSrc = axang2rotm([perpenDir', DirTheta]);

disp('Saving results...');

save(fullfile(resultDir, 'pt_light.mat'), 'locLightSrc', 'rotLightSrc');

display(locLightSrc);
display(dirLightSrc);
display(rotLightSrc);

%% Setting up light simulator figure

%radiant intensity distribution parameters
maxRadiInt = 10;
mu = 1.5;
%attentuation
r_d = 1; %effective radius of disk source

%Starting distance from the source in metres for visualisation
distFromSrc = 0.2;



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
lightSrc = LightSimulator(locLightSrc, rotLightSrc, figSim);
