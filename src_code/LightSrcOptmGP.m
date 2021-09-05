% function [mu, varMu, hypOpt] = LightSrcOptmGP(meanType, trainingData, testingX, sigmaNoise, downSamp, iter, verboseFlag, varargin)

function [mu, varMu, hypOpt, GPDataDown] = LightSrcOptmGP(meanType, trainingData, testingX, sigmaNoise, varargin)
% Builds the intensity field of a non-isotropic disk light source model
% given data using Gaussian Processes

% INPUTS:
%       meanType - Type of GP mean function. 0 - zero mean, 1 - constant mean, 2 - light source mean function (DEFAULT)
%       trainingData - training data [radius, theta, radiant intensity magnitude]
%       testingX - Testing data [radius, theta]
%       downSamp - Amount of downsampling for the given training data
%       iter - maximum number of iterations during optimisation
%       verboseFlag - verbose flag to show output to terminal (true/false)

% OUTPUTS:
%       mu - Queried radiant intensity mean of testing data from optimised
%           GP model
%       varMu - Queried radiant intensity variance of testing data from
%           optimised GP model
%       hypOpt - optimised hyperparameter struct

%Author: Jasprabhjit Mehami, 13446277

covfunc = @covSEiso; %covariance function
likfunc = @likGauss; %gaussian likelihood

switch (meanType)
    %zero-mean
    case 0
        meanfunc = [];
        meanHyper = [];
        %constant mean
    case 1
        meanfunc = @meanConst;
        meanHyper = 1;
        %light source mean function
    case 2
        meanfunc = @meanLightSrc;
        meanHyper = [1,1,1];
    case 3
        meanfunc = @meanLightSrcExp;
        meanHyper = [1,1,1];
    otherwise
        meanfunc = @meanLightSrc;
        meanHyper = [1,1,1];
end

%hyperparameter struct
hyp = struct('mean', meanHyper, 'cov', [0,0], 'lik', log(sigmaNoise));

%only querying optimised model
if nargin == 5
    xGP = trainingData(:,1:2);
    yGP = trainingData(:,3);
    
    hypOpt = varargin{1};
%training of GP and querying optimised model
elseif nargin == 7
    downSamp = varargin{1};
    iter = varargin{2};
    verboseFlag = varargin{3};
    
    %downsample data
    GPDataDown = downsample(trainingData, downSamp);
    xGP = GPDataDown(:,1:2);
    yGP = GPDataDown(:,3);
    
    hypOpt = minimize(hyp, @gp, -iter, verboseFlag, @infGaussLik, meanfunc, covfunc, likfunc, xGP, yGP);
else
    error('There should be either 5 or 7 input arguments');
end

%use the optimised hyperparameter struct to get the mean and variance of
%given testing points
[mu, varMu] = gp(hypOpt, @infGaussLik, meanfunc, covfunc, likfunc, xGP, yGP, testingX);

%no imaginery elements
mu = real(mu);

end

