%Test the ArucoPixDect mex function.

%Author: Jasprabhjit Mehami, 13446277

clear;
close all;

%mex function should be in the bin folder
if exist("bin", "dir")
    addpath("bin");
else
    error("mex file not built");
end

%robotics toolbox for visualising
run(['rvctools' filesep 'startup_rvc.m']);

%test image
img = imread(['Images', filesep,'aruco_example.png']);

%camera parameters
intrMat = [532.568131996427,0,0;0,531.905416600879,0;327.499527166381,231.227840418968,1]; %intrinsic matrix for opencv format
distRad = [0.0346875042867809,-0.0917743770901257,-0.0897944587524139];
distTan = [-0.00415109739624088,0.00571543700759848];
distCoefCV = [distRad(1:2), distTan, distRad(3)]; %array of distortion coefficients in opencv format


camParam = cameraParameters('IntrinsicMatrix', intrMat, 'ImageSize', size(img, [1,2]), ...
    'RadialDistortion', distRad, 'TangentialDistortion', distTan);
%% Get pattern points of Aruco board

%ChArUco pattern size
xNumMarker = 8;
yNumMarker = 5;
arucoLen = 0.04;
sepLen = 0.01;

markerCornerCell = cell(yNumMarker*xNumMarker, 1);
ind = 0;

markerIDs = zeros(yNumMarker*xNumMarker, 1);

for i = 1:yNumMarker
    for j = 1:xNumMarker
        x_tl = (j-1)*arucoLen + (j - 1)*sepLen;
        y_tl = (yNumMarker - i + 1)*arucoLen + (yNumMarker - i)*sepLen;
        
        x_tr = j*arucoLen + (j - 1)*sepLen;
        y_tr = (yNumMarker - i + 1)*arucoLen + (yNumMarker - i)*sepLen;
        
        x_br = j*arucoLen + (j - 1)*sepLen;
        y_br = (yNumMarker - i)*arucoLen + (yNumMarker - i)*sepLen;
        
        x_bl = (j-1)*arucoLen + (j - 1)*sepLen;
        y_bl = (yNumMarker - i)*arucoLen + (yNumMarker - i)*sepLen;
        
        markerCorner = [
            x_tl, y_tl, 0;
            x_tr, y_tr, 0;
            x_br, y_br, 0;
            x_bl, y_bl, 0;
            ];
        
        markerIDs(ind + 1) = ind;
        ind = ind + 1;
        markerCornerCell(ind) = {markerCorner};
    end
end

% remove the markers that do not exist on the board
for i = 18:2:23
    markerCornerCell(i:i+5) = [];
    markerIDs(i:i+5) = [];
end

numMarkers = length(markerIDs);
markerPatPts = zeros(numMarkers*4, 3);

figure(1);


for i = 1:numMarkers
    curCorner =  markerCornerCell{i};
    
    plot([curCorner(:,1); curCorner(1,1)], [curCorner(:,2); curCorner(1,2)], 'b-'); hold on;
    
    markerPatPts((i-1)*4 + 1: 4*i, :) = curCorner;
end

grid on;
axis equal;




figB = figure('Name', 'Before');
imshow(img);

imgGray = im2gray(img);


[ids, markerCorner, imgOut] = ArucoPixDect(imgGray);


figure('Name', 'After');
imshow(imgOut);

if length(markerIDs) ~= length(ids)
    error('did not find all markers');
end

[ids, indx] = sort(ids, 'ascend');

markCornerSort = markerCorner(indx, :);
numMarkers = length(ids);

markerImgPts = zeros(4*numMarkers, 2);


for i = 1:length(ids)
    corners = markCornerSort(i,:);
    
    uPts = corners(1:2:end);
    vPts = corners(2:2:end);
    
    
    markerImgPts((i-1)*4 + 1 : 4*i, :) = [uPts', vPts'];
end

figure(figB);
hold on;
for i = 1:size(markerImgPts, 1)
    plot(markerImgPts(i,1), markerImgPts(i,2), 'r*');
end


[rotMat, trans] = extrinsics(markerImgPts, markerPatPts(:, 1:2), camParam);


patFig = figure('Name', 'Pose of pattern');

%plot frame camera
plotCamera('Location', zeros(1,3), 'Orientation', eye(3), 'Size', 0.05, 'AxesVisible', true); hold on;

ext = [rotMat,trans'; 0, 0, 0, 1];

trplot(ext, 'frame', 'Pat', 'rviz', 'length', 0.1); hold on;

% bl = [0, 0, 0, 1];
% br = [xNumCheck*checkSize, 0, 0, 1];
% tl = [0, yNumCheck*checkSize, 0, 1];
% tr = [xNumCheck*checkSize, yNumCheck*checkSize, 0, 1];

coor = ext*[markerPatPts'; ones(1, numMarkers*4)];

plot3(coor(1,:),coor(2,:),coor(3,:),'r*')

axis equal; grid on;
xlabel("tx (m)");
ylabel("ty (m)");
zlabel("tz (m)");

