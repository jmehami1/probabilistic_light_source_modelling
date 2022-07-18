% Calculated estimated reflectance of simulated data from Mitsuba render.

close all;
clear;

%external library directory
addpath('ext_lib');

%robotics toolbox
run(fullfile('ext_lib', 'rvctools', 'startup_rvc.m'));

%Contains implementation of dichromatic model estimation by Huynh and
%Robles-Kelly
addpath(genpath(fullfile('ext_lib', 'Scyllarus')));

addpath(genpath('krebs-shape-reflectance-estimation'));
addpath(genpath(fullfile('ext_lib', 'krebs-reflectance-estimation')));

addpath("openexr-matlab-master");

%code for this project
addpath('src_code');

numChannels = 100;

%% Load data

datasetName = 'bunny_mixture';

mainDir = fullfile('/home/jmeh/PAPER_SCENES', datasetName);
resultDir = fullfile(mainDir, 'results');

if ~exist(resultDir, 'dir')
    mkdir(resultDir);
end


%file paths
rgbPathEXR = fullfile(mainDir, 'rgb.exr');
specReflEXR = fullfile(mainDir, 'spectral_reflectance.exr');
specRadiEXR = fullfile(mainDir, 'spectral_radiant-intensity.exr');

dataMap = exrreadchannels(rgbPathEXR);
% infoExr = exrinfo(rgbPathEXR);

figure();
tiledlayout(2,3, 'Padding', 'tight');



%read color image
curLayerEXR = "color";

imgChannelsCell = values(dataMap,{curLayerEXR+".R", curLayerEXR+".G", curLayerEXR+".B"});

imgRGB = zeros([size(imgChannelsCell{1}),3]);
for i=1:length(imgChannelsCell)
    imgRGB(:,:,i) = imgChannelsCell{i};
end
imgRGB = rescale(imgRGB, 0, 1);

nexttile;
imshow(imgRGB); hold on;
title(curLayerEXR);
imwrite(imgRGB, fullfile(resultDir, datasetName + "_" + curLayerEXR + ".png"));


%read ray distance image
curLayerEXR = "distance";
imgChannelsCell = values(dataMap,{curLayerEXR+".Y"});
imgRayDist = cell2mat(imgChannelsCell);
imgRayDist(imgRayDist < 0.1) = nan;

nexttile;
imshow(imgRayDist, 'DisplayRange', [min(imgRayDist, [], 'all'), max(imgRayDist, [], 'all')])
title(curLayerEXR);
colorbar;

%read normal image
curLayerEXR = "normal";
imgChannelsCell = values(dataMap,{curLayerEXR+".X", curLayerEXR+".Y", curLayerEXR+".Z"});
imgSurfNorm = zeros([size(imgChannelsCell{1}),3]);
for i=1:length(imgChannelsCell)
    imgSurfNorm(:,:,i) = imgChannelsCell{i};
end

nexttile;
imshow((imgSurfNorm + 1)./2);
title(curLayerEXR);
% imwrite(imgSurfNorm, fullfile(resultDir, datasetName + "_" + curLayerEXR + ".png"));


%read point position relative to camera image
curLayerEXR = "position_camera";
imgChannelsCell = values(dataMap,{curLayerEXR+".X", curLayerEXR+".Y", curLayerEXR+".Z"});
imgPtCam = zeros([size(imgChannelsCell{1}),3]);
for i=1:length(imgChannelsCell)
    imgPtCam(:,:,i) = imgChannelsCell{i};
end

nexttile;
imshow(rescale(imgPtCam, 0, 1));
title(curLayerEXR);
imwrite(imgSurfNorm, fullfile(resultDir, datasetName + "_" + curLayerEXR + ".png"));

%read point position relative to world image
% curLayerEXR = "position_world";
% imgChannelsCell = values(dataMap,{curLayerEXR+".X", curLayerEXR+".Y", curLayerEXR+".Z"});
% imgPtWorld = zeros([size(imgChannelsCell{1}),3]);
% for i=1:length(imgChannelsCell)
%     imgPtWorld(:,:,i) = imgChannelsCell{i};
% end
% 
% nexttile;
% imshow(rescale(imgPtWorld, 0, 1));
% title(curLayerEXR);


dataMap = exrreadchannels(specReflEXR);
infoExr = exrinfo(specReflEXR);


%read radiance hypercube
% curLayerEXR = "radiance";
keySet = infoExr.channels(1:numChannels);
imgChannelsCell = values(dataMap,keySet);

hyperRadiance = zeros([size(imgChannelsCell{i}), numChannels]);

for i=1:length(imgChannelsCell)
    hyperRadiance(:,:,i) = imgChannelsCell{i};
end


%read reflectance hypercube
% curLayerEXR = "reflectance";
keySet = infoExr.channels(numChannels+1:end);
imgChannelsCell = values(dataMap,keySet);

hyperReflectance = zeros([size(imgChannelsCell{i}), numChannels]);

for i=1:length(imgChannelsCell)
    hyperReflectance(:,:,i) = imgChannelsCell{i};
end

dataMap = exrreadchannels(specRadiEXR);
infoExr = exrinfo(specRadiEXR);

%read radiant intensity hypercube
curLayerEXR = "radiant_Intensity";
keySet = infoExr.channels(numChannels+1:end);
imgChannelsCell = values(dataMap,keySet);

hyperRadiIntensity = zeros([size(imgChannelsCell{i}), numChannels]);

for i=1:length(imgChannelsCell)
    hyperRadiIntensity(:,:,i) = imgChannelsCell{i};
end

clear dataMap imgChannelsCell keySet infoExr

%% Setup for estimated reflectance methods

%position of camera and light in world coordinate frame.
posLightSrc = [-0.1, 0.05, -0.05];
imgdirLightVec = zeros(1,1,3);

for i = 1:length(posLightSrc)
    imgdirLightVec(:,:,i) = posLightSrc(i);
end

imgdirLightVec = minus(imgdirLightVec, imgPtCam);
imgdirLightVec = imgdirLightVec./vecnorm(imgdirLightVec,2,3);
lDotnImg = dot(imgSurfNorm, imgdirLightVec, 3);
binMaskG = (lDotnImg > 0) & (all(hyperRadiance,3)) & (imgRayDist > 0.1);
binMaskG = bwmorph(binMaskG, 'clean', Inf);
binMaskG = bwmorph(binMaskG, 'spur', Inf);
binMaskG = bwmorph(binMaskG, 'clean', Inf);


lDotnImg(~binMaskG) = 0;

imgSurfNorm = bsxfun(@times, imgSurfNorm, binMaskG);
imgdirLightVec = bsxfun(@times, imgdirLightVec, binMaskG);

hyperRadiIntensity = bsxfun(@times, hyperRadiIntensity, lDotnImg);


figure();
tiledlayout(1,3, 'Padding', 'tight');

nexttile;
imshow((imgSurfNorm + 1)./2);
title("Surface normal");
imwrite((imgSurfNorm + 1)./2, fullfile(resultDir, datasetName + "_normal.png"));

nexttile;

imshow((imgdirLightVec + 1)./2);

title("Normalised light vector");
imwrite((imgdirLightVec + 1)./2, fullfile(resultDir, datasetName + "_light_source_direction.png"));


nexttile;
imshow(lDotnImg);
title("dot product");
imwrite(lDotnImg, fullfile(resultDir, datasetName + "_l_dot_n.png"));


%% estimated reflectance methods
close all;

%false colour channels
channels = [20,40,80];

methods = ["Robles"; "Krebs"; "KrebsShape"];

figure();
tiledlayout(1,4, 'Padding', 'tight');


%***Robles***
method = methods(1);
hyperRadianceMasked = bsxfun(@times, hyperRadiance, binMaskG);
hyperReflectanceMasked = bsxfun(@times, hyperReflectance, binMaskG);
[k, g, S] = recover_dichromatic_parameters_LS(hyperRadianceMasked, hyperRadiIntensity, 0, 5, 1);
resTable = SAM_Hypercube(S, hyperReflectanceMasked);

g(isnan(g)) = 0;

falseImg = zeros([size(S, [1,2]), 3]);

for i = 1:3
    falseChanGray = S(:,:,channels(i));
    falseImg(:,:,i) = falseChanGray;
end

% falseRow = reshape(falseImg, 1, []);

falseImg = rescale(falseImg, 'InputMax', 1);

nexttile;
imshow(falseImg);
title(method);

imwrite(falseImg, fullfile(resultDir, datasetName + "_" + lower(method) + "_refl.png"));
imwrite(mat2gray(k), fullfile(resultDir, datasetName + "_" + lower(method) + "_k.png"));
imwrite(mat2gray(g), fullfile(resultDir, datasetName + "_" + lower(method) + "_g.png"));

%***Krebs***
method = methods(2);
R_hyperCube = hyperRadianceMasked./hyperRadiIntensity;
R_hyperCube(isnan(R_hyperCube)) = 0; %ignore any nan
R_hyperCube(isinf(R_hyperCube)) = 0; %ignore any inf

[S, g, k] = method_krebs_logversion(R_hyperCube, binMaskG);
resTableTemp = SAM_Hypercube(S, hyperReflectanceMasked);
resTable = [resTable; resTableTemp];

falseImg = rescale(falseImg);

for i = 1:3
    falseChanGray = S(:,:,channels(i));
    falseImg(:,:,i) = falseChanGray;
end

falseImg = rescale(falseImg);

nexttile;
imshow(falseImg);
title(method);

imwrite(falseImg, fullfile(resultDir, datasetName + "_" + lower(method) + "_refl.png"));
imwrite(mat2gray(k), fullfile(resultDir, datasetName + "_" + lower(method) + "_k.png"));
imwrite(mat2gray(g), fullfile(resultDir, datasetName + "_" + lower(method) + "_g.png"));

%***Krebs with Shape***
method = methods(3);
[S, k, g] = ReflectanceEstimationKrebs(hyperRadianceMasked,hyperRadiIntensity, binMaskG, lDotnImg, 3, 0.01);
resTableTemp = SAM_Hypercube(S, hyperReflectanceMasked);
resTable = [resTable; resTableTemp];


for i = 1:3
    falseChanGray = S(:,:,channels(i));
    falseImg(:,:,i) = falseChanGray;
end

falseImg = rescale(falseImg);

nexttile;
imshow(falseImg);
title(method);

imwrite(falseImg, fullfile(resultDir, datasetName + "_" + lower(method) + "_refl.png"));
imwrite(mat2gray(k), fullfile(resultDir, datasetName + "_" + lower(method) + "_k.png"));
imwrite(mat2gray(g), fullfile(resultDir, datasetName + "_" + lower(method) + "_g.png"));


%***Reflectance Groundtruth***
method = "gt";
S = hyperReflectanceMasked;
for i = 1:3
    falseChanGray = S(:,:,channels(i));
    falseImg(:,:,i) = falseChanGray;
end

falseImg = rescale(falseImg);

nexttile;
imshow(falseImg);
title(method);

imwrite(falseImg, fullfile(resultDir, datasetName + "_" + lower(method) + "_refl.png"));

resTableTemp = table(methods, 'VariableNames', "Estimation Method");
resTable = [resTableTemp,resTable];


%% Results output

disp(resTable);

fprintf("Number of actual pixels %i\n", sum(binMaskG, 'all'));

writetable(resTable, fullfile(resultDir, "summary.csv"));
