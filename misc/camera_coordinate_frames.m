close all;
clear;

%robotics toolbox
run(['rvctools' filesep 'startup_rvc.m'])




%% Frame camera setup (world coordinate frame)

figFrame = figure('Name', 'Frame world coordinate frame');


rotW = eye(3);
tW = [0;0;0];

plotCamera('Location', tW, 'Orientation', rotW, 'Size', 0.05, 'AxesVisible', true, 'Color', [1,0,0]); hold on;
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');


%% Line-scan camera setup

rotLS = eul2rotm(deg2rad([90,0,0]), 'ZYX');
tLS = 0.5.*[1;0.8;0.4];

poseLS = [rotLS, tLS; 0,0,0,1];

plotCamera('Location', tLS, 'Orientation', rotLS', 'Size', 0.05, 'AxesVisible', true, 'Color', [0,0,1]); hold on;


%% Point in space

pt = [0.2; 0.2; 0.4];
scatter3(pt(1), pt(2), pt(3), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

%% line-scan camera is world coordinate frame

figLine = figure('Name', 'Line-scan world coordinate frame');

plotCamera('Location', tW, 'Orientation', rotW, 'Size', 0.05, 'AxesVisible', true, 'Color', [0,0,1]); hold on;
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');
% Line-scan pose


%pose of frame camera w.r.t to line-scan, is equal to inverse of pose of
%line-scan w.r.t to frame, which is also equal to the extrinsic of the
%line-scan camera
poseFR = inv(poseLS);
rotFR = tform2rotm(poseFR);
tFR = tform2trvec(poseFR);

plotCamera('Location', tFR, 'Orientation', rotFR', 'Size', 0.05, 'AxesVisible', true, 'Color', [1,0,0]); hold on;


ptLS_Hom = poseFR*[pt;1];
ptLS = ptLS_Hom(1:3)

ptLS = rotFR*pt + tFR'
scatter3(ptLS(1), ptLS(2), ptLS(3), 200, [1,0,0], 'filled', 'MarkerEdgeColor', [0,0,0]);

