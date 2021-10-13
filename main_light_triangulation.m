%Triangulation of real point source's principal direction and location
%relative to a frame camera origin. This uses a reflective metal sphere on
%a ChArUco board. By fitting an ellipse to the brightest spot on the sphere
%we can find a constraint to determine the unknown direction and location.

%Author: Jasprabhjit Mehami, 13446277

close all;
clear;

%robotics toolbox for visualising axes
run(['rvctools' filesep 'startup_rvc.m']);

%yaml reader package
addpath(genpath('yamlmatlab-master'));

%parameter file location
addpath(genpath('parameter_files'));

%source code for this project
addpath('src_code');

%MEX functions for ArUco pose estimation
addpath(genpath(['ext_lib', filesep, 'mex_aruco']));


%parameter file
paramFile = ['parameter_files', filesep, 'light_triangulation.yaml'];
if ~exist(paramFile, 'file')
    error("YAML parameter file not found");
end


%% Read YAML file containing the pattern specifications and parameters for code
% All dimensions are in metres

paramYaml = yaml.ReadYaml(paramFile);
displayOn = paramYaml.DisplayOn;
usingBlackfly = paramYaml.UsingBlackFly;


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

%Get source directory where images are located and results will be saved
sourDir = uigetdir(['~', filesep], 'Provide source directory where light triangulation images are located?');

%frame image directory
frameDir = [sourDir, filesep, 'Frame']; %Directory containing images
if ~exist(frameDir, 'dir')
    error('Source directory does not contain directory called "Frame" which should contain RGB images');
end

fullPathFrame = fullfile([frameDir, filesep, '*.png']);

%Need to get number of images in directory
numImages = numel(dir(fullPathFrame));

if numImages < 1
   error('no images in provided image directory') 
end

%result directory
resultDir = [sourDir, filesep, 'Result'];
if ~exist(resultDir, 'dir')
   mkdir(resultDir);
end

%% Load images

%Preallocate space for cell arrays of images
imagesFrame = cell(1,numImages);

% Load all images
for i = 1:numImages
    imagesFrame{i} = imread([frameDir, filesep, 'img', num2str(i),'.png']);
end


%% Get pose of the plane using the ChArUco pattern
disp('Getting pose of ChArUco pattern...');

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
std_f = estimationErrors.IntrinsicsErrors.FocalLengthError;
std_u0v0 = estimationErrors.IntrinsicsErrors.PrincipalPointError;
stdRGB_Intr = [std_f,std_u0v0];

%ChArUco pattern size
xNumCheck = paramYaml.NumCols;
yNumCheck = paramYaml.NumRows;
checkSize = paramYaml.CheckerSideLength;
arucoSize = paramYaml.ArucoSideLength;

%intrinsic object for the RGB camera
frameIntrinsic = cameraIntrinsics(thetaFrameintr(1:2),thetaFrameintr(3:4), frameImgSize);
K_frame = frameIntrinsic.IntrinsicMatrix;

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
    imagesFrame{i} = undistortImage(imagesFrame{i}, cameraParams);
    
    
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


for curImg = 1:numImages
    if displayOn
        clf(figG);
    end
    
    %Get the pixel coordinates of the edge of the reflective sphere on the
    %board
    ext = extPosePattern(:,:,curImg);
    patPtsHomFrame = ext*patPtsHom;
    patPtsFrame = patPtsHomFrame(1:3,:)';
    patImgPts = projectPoints(patPtsFrame, K_frame', eye(4), distCoefCV, frameImgSize);
    
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
    
    img = imagesFrame{curImg};
    imgGray = im2gray(img);
    imgPatGr = uint8(zeros(size(imgGray)));
    
    %extract on the pixels that contain the reflective sphere
    for v = min(patImgPts(2, :), [], 'all'):max(patImgPts(2, :), [], 'all')
        for u = min(patImgPts(1, :), [], 'all'):max(patImgPts(1, :), [], 'all')
            if (inpolygon(u, v, patImgPts(1,:), patImgPts(2,:)))
                imgPatGr(v, u) = imgGray(v, u);
            end
        end
    end
    
    
    %thresold
    imgBrSpot = imgPatGr > 250;
    
    %erode
    SE = strel('disk',1);
    imgBrSpot = imerode(imgBrSpot, SE);
    
    %find regions
    s = regionprops(imgBrSpot, {'Centroid','Orientation','MajorAxisLength','MinorAxisLength', 'Area'});
    
    %find the region with the largest area
    ind = 1;
    maxArea = s(1).Area;
    
    for i = 2:length(s)
        if s(i).Area > maxArea
            ind = i;
            maxArea = s(i).Area;
        end
    end
    
    %get the centre of the region
    regCent = s(ind).Centroid;
    reflCentImg(curImg, :) = regCent;
    
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

%% Find location of brighest spot on reflective sphere in world coordinates

disp('Finding location of brightest spot in world coordinates...');

%centre of the hemisphere relative to the pattern
hemiCent = [xCent; yCent; 0];

ptReflHemiFrame = zeros(numImages, 3);
normReflHemi = zeros(numImages, 3);

for curImg = 1:numImages
    centPtImg =  reflCentImg(curImg, :)';
    centPtImgHom = [centPtImg; 1];
    
    %normalised image coordinates
    normPt = K_frame'\centPtImgHom;
    normPt = normPt./normPt(3,:);
    
    %centre of hemisphere relative to the frame camera
    ext = extPosePattern(:,:,curImg);
    locPatFrame = tform2trvec(ext)';
    rotPatFrame = tform2rotm(ext);
    locHemiCentFrame = locPatFrame + rotPatFrame*hemiCent;
    
    [pt, rc] = intersect_Line_Sphere([normPt', normPt'], [locHemiCentFrame', hemiRad]);
    
    if ~rc
        error("line-sphere intersection was not on sphere")
    end
    
    ptReflHemiFrame(curImg, :) = pt;
    
    %surface normal at the intersection point on the reflective hemisphere
    normHemi = pt' - locHemiCentFrame;
    normHemi = normHemi./norm(normHemi);
    normReflHemi(curImg, :) = normHemi;
end

%% Find the light source direction and location

disp('Finding point source direction and location...');

[locLightSrc, dirLightSrc] = SolveLightSrcDirLoc(ptReflHemiFrame, normReflHemi, numImages);

%principal direction is along the z-axis, find how much the z-axis has
%rotated from its original direction. This rotation would be about some
%perpindular axes
zAxis = [0;0;1];
perpenDir = cross(zAxis, dirLightSrc);
DirTheta = atan2(norm(cross(dirLightSrc,zAxis)),dot(dirLightSrc,zAxis));
rotLightSrc = axang2rotm([perpenDir', DirTheta]);

resultYaml.locLightSrc = locLightSrc;
resultYaml.rotLightSrc = rotLightSrc;
yaml.WriteYaml([resultDir, filesep, 'pt_light.yaml'], resultYaml);

display(locLightSrc);
display(dirLightSrc);


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
lightSrc = LightSimulator(locLightSrc, rotLightSrc, maxRadiInt, mu, r_d, distFromSrc, figSim, S);
