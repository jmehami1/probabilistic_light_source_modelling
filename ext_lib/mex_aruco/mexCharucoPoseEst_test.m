%Test the CharucoPoseEst mex function.

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
img = imread(['Images', filesep,'patexample.png']);

%camera parameters
intrMat = [532.568131996427,0,0;0,531.905416600879,0;327.499527166381,231.227840418968,1]; %intrinsic matrix for opencv format
distRad = [0.0346875042867809,-0.0917743770901257,-0.0897944587524139];
distTan = [-0.00415109739624088,0.00571543700759848];
distCoefCV = [distRad(1:2), distTan, distRad(3)]; %array of distortion coefficients in opencv format

%ChArUco pattern size
xNumCheck = 8;
yNumCheck = 6;
checkSize = 0.04;
arucoSize = 0.03;

figure('Name', 'Before');
imshow(img);

%mex function for getting pose of the pattern
[rotMat, trans, found, imgOut] = CharucoPosEst(img, intrMat, distCoefCV, xNumCheck, yNumCheck, checkSize, arucoSize);

figure('Name', 'After');
imshow(imgOut);

%if found valid pose, plot the pattern with the camera using
if found
    patFig = figure('Name', 'Pose of pattern');
    
    %plot frame camera
    plotCamera('Location', zeros(1,3), 'Orientation', eye(3), 'Size', 0.05, 'AxesVisible', true); hold on;
    
    ext = [rotMat,trans'; 0, 0, 0, 1];
    
    trplot(ext, 'frame', 'Pat', 'rviz', 'length', 0.1); hold on;
    
    bl = [0, 0, 0, 1];
    br = [xNumCheck*checkSize, 0, 0, 1];
    tl = [0, yNumCheck*checkSize, 0, 1];
    tr = [xNumCheck*checkSize, yNumCheck*checkSize, 0, 1];
    
    coor = ext*([bl ; br; tr; tl ]') ;
    
    fill3(coor(1,:),coor(2,:),coor(3,:),'r', 'FaceAlpha', 0.5)
    
    axis equal; grid on;
    xlabel("tx (m)");
    ylabel("ty (m)");
    zlabel("tz (m)");
    
end
