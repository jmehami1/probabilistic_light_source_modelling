clear all, close all, clc

I = imread('apple/original.png');
mask = imread('apple/mask.png');
[S,g,k] = method_krebs_logversion(double(I), mask);
figure, imagesc(g/prctile(g(:),99));
figure, imagesc(S/prctile(S(:),99));
figure, imagesc(k/prctile(k(:),99));